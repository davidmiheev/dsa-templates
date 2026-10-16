//! Convolutions by the number-theoretic transform (NTT).
//!
//! * [`convolution`]: `c[k] = sum a[i] * b[k - i] mod 998244353`.
//! * [`convolution_mod`]: the same for any supported NTT prime (a const generic), when the result
//!   fits the transform size `2^v2(p - 1)`.
//! * [`convolution_arbitrary_mod`]: any modulus below `2^31`, by three primes and Garner's
//!   algorithm (results up to about `7.9e25 / m^2` terms stay exact).
//! * [`convolution_large`]: modulo 998244353 for results longer than the `2^23` the prime itself
//!   supports (block decomposition with `2^23` transforms).
//! * [`convolution_u64`]: the *exact* result when every coefficient of the true convolution is below
//!   `2^64` (computed modulo three primes and recombined), e.g. for counting problems.
//! * [`convolution_mod_2_64`]: the convolution of `u64` sequences modulo `2^64` (wrapping), by
//!   splitting into four 16-bit limbs.
//! * [`ntt`] / [`intt`]: the raw transforms (`u32` slices modulo a prime), for crates that need
//!   to combine spectra themselves, e.g. polynomial algorithms.
//!
//! The transforms are radix-4, without bit reversal (decimation in frequency forward, in time
//! inverse). Small inputs fall back to the schoolbook product.
//!
//! # Example
//!
//! ```
//! use convolution::{convolution, convolution_arbitrary_mod, convolution_u64};
//!
//! assert_eq!(convolution(&[1, 2, 3], &[4, 5]), vec![4, 13, 22, 15]);
//! assert_eq!(convolution_u64(&[1 << 40, 1], &[1 << 20, 1]), vec![1 << 60, (1 << 40) + (1 << 20), 1]);
//! assert_eq!(convolution_arbitrary_mod(&[1_000_000_006; 3], &[1_000_000_006; 3], 1_000_000_007), vec![1, 2, 3, 2, 1]);
//! ```

mod ntt;

pub use ntt::{intt, intt_unscaled, ntt, primitive_root};
use ntt::pow_mod;

/// 998244353 = 119 * 2^23 + 1, primitive root 3.
pub const P1: u32 = 998_244_353;
/// 167772161 = 5 * 2^25 + 1, primitive root 3.
pub const P2: u32 = 167_772_161;
/// 469762049 = 7 * 2^26 + 1, primitive root 3.
pub const P3: u32 = 469_762_049;
/// 2013265921 = 15 * 2^27 + 1, primitive root 31.
pub const P4: u32 = 2_013_265_921;

fn naive<const M: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut c = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            c[i + j] = (c[i + j] + x as u64 * y as u64) % M as u64;
        }
    }
    c.into_iter().map(|x| x as u32).collect()
}

/// Convolution modulo the NTT-friendly prime `M` ([`P1`]..[`P4`], or `754974721`, `924844033`,
/// `1004535809`). The inputs must be reduced (`< M`).
///
/// # Panics
/// Panics for an unsupported modulus or if the result length exceeds the largest power of two
/// dividing `M - 1`.
pub fn convolution_mod<const M: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let len = a.len() + b.len() - 1;
    if a.len().min(b.len()) <= 24 {
        return naive::<M>(a, b);
    }
    let n = len.next_power_of_two();
    let mut fa = vec![0u32; n];
    let mut fb = vec![0u32; n];
    fa[..a.len()].copy_from_slice(a);
    fb[..b.len()].copy_from_slice(b);
    ntt::ntt::<M>(&mut fa);
    ntt::ntt::<M>(&mut fb);
    for (x, &y) in fa.iter_mut().zip(&fb) {
        *x = (*x as u64 * y as u64 % M as u64) as u32;
    }
    ntt::intt::<M>(&mut fa);
    fa.truncate(len);
    fa
}

/// `c[k] = sum_{i + j = k} a[i] * b[j] mod 998244353`. Inputs must be `< 998244353`.
pub fn convolution(a: &[u32], b: &[u32]) -> Vec<u32> {
    convolution_mod::<P1>(a, b)
}

/// Garner's algorithm for three coprime primes `p`, with the inverses precomputed.
struct Garner3 {
    p: [u64; 3],
    inv12: u64,
    inv123: u64,
    p01_mod_m: u64,
    m: u64,
}

impl Garner3 {
    fn new(p: [u64; 3], m: u64) -> Self {
        Self {
            p,
            inv12: pow_mod(p[0], p[1] - 2, p[1]),
            inv123: pow_mod(p[0] * p[1] % p[2], p[2] - 2, p[2]),
            p01_mod_m: p[0] * p[1] % m,
            m,
        }
    }

    /// The value with residues `r` modulo the three primes, reduced modulo `m`.
    fn combine(&self, r: [u32; 3]) -> u32 {
        let p = &self.p;
        let x1 = r[0] as u64;
        let x2 = (r[1] as u64 + p[1] - x1 % p[1]) % p[1] * self.inv12 % p[1];
        let x12 = x1 + p[0] * x2; // < p0 * p1 < 2^64
        let x3 = (r[2] as u64 + p[2] - x12 % p[2]) % p[2] * self.inv123 % p[2];
        (((x12 % self.m) + self.p01_mod_m * x3) % self.m) as u32
    }
}

/// Convolution modulo `m` computed modulo the three NTT primes `A`, `B`, `C` and recombined by
/// Garner's algorithm. Exact when the true coefficients (below `len * m^2`) are less than
/// `A * B * C`; the inputs must be below `m`.
pub fn convolution_three_primes<const A: u32, const B: u32, const C: u32>(a: &[u32], b: &[u32], m: u32) -> Vec<u32> {
    let reduce = |v: &[u32], p: u32| -> Vec<u32> { v.iter().map(|&x| x % p).collect() };
    let c1 = convolution_mod::<A>(&reduce(a, A), &reduce(b, A));
    let c2 = convolution_mod::<B>(&reduce(a, B), &reduce(b, B));
    let c3 = convolution_mod::<C>(&reduce(a, C), &reduce(b, C));
    let g = Garner3::new([A as u64, B as u64, C as u64], m as u64);
    (0..c1.len()).map(|k| g.combine([c1[k], c2[k], c3[k]])).collect()
}

/// Convolution modulo any `m < 2^31` (inputs must be `< m`), by three NTT primes. Exact as long
/// as `min(len a, len b) * (m - 1)^2 < 7.9e25`, which holds for `len <= 2^23` and `m <= 10^9 + 7`.
///
/// # Panics
/// Panics if the result is longer than `2^23`.
pub fn convolution_arbitrary_mod(a: &[u32], b: &[u32], m: u32) -> Vec<u32> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    assert!(m >= 2 && m < 1 << 31, "modulus {m} out of range");
    if a.len().min(b.len()) <= 24 {
        let mut c = vec![0u64; a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                c[i + j] = (c[i + j] + x as u64 * y as u64) % m as u64;
            }
        }
        return c.into_iter().map(|x| x as u32).collect();
    }
    convolution_three_primes::<P1, P2, P3>(a, b, m)
}

/// Convolution modulo 998244353 of sequences of **any** length. The prime only supports transforms
/// of size `2^23`, so long inputs are cut into blocks of `2^22`: every block is transformed once
/// (size `2^23`), the spectra of block pairs are multiplied and accumulated per output block
/// `i + j`, and the inverse transforms are overlap-added. Memory is about
/// `(blocks_a + blocks_b + blocks_a + blocks_b - 1) * 32 MB`. Inputs must be `< 998244353`.
pub fn convolution_large(a: &[u32], b: &[u32]) -> Vec<u32> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    const SIZE: usize = 1 << 23;
    const BLOCK: usize = SIZE / 2;
    if a.len() + b.len() - 1 <= SIZE {
        return convolution(a, b);
    }
    let spectra = |v: &[u32]| -> Vec<Vec<u32>> {
        v.chunks(BLOCK)
            .map(|chunk| {
                let mut f = vec![0u32; SIZE];
                f[..chunk.len()].copy_from_slice(chunk);
                ntt::ntt::<P1>(&mut f);
                f
            })
            .collect()
    };
    let (fa, fb) = (spectra(a), spectra(b));
    let mut c = vec![0u32; a.len() + b.len() - 1];
    for k in 0..fa.len() + fb.len() - 1 {
        // pairs (i, k - i) of blocks whose product lands in output block k; at most 16 terms of
        // < 2^60 fit a u64, so a single reduction per point suffices
        let pairs: Vec<(&[u32], &[u32])> = (k.saturating_sub(fb.len() - 1)..=k.min(fa.len() - 1))
            .map(|i| (&fa[i][..], &fb[k - i][..]))
            .collect();
        assert!(pairs.len() <= 16, "too many blocks per output block");
        let mut block = vec![0u32; SIZE];
        for (t, out) in block.iter_mut().enumerate() {
            let mut acc = 0u64;
            for &(x, y) in &pairs {
                acc += x[t] as u64 * y[t] as u64;
            }
            *out = (acc % P1 as u64) as u32;
        }
        ntt::intt::<P1>(&mut block);
        let offset = k * BLOCK;
        for (t, &v) in block.iter().enumerate().take(c.len().saturating_sub(offset)) {
            let slot = &mut c[offset + t];
            *slot = {
                let s = *slot + v;
                if s >= P1 {
                    s - P1
                } else {
                    s
                }
            };
        }
    }
    c
}

/// The exact convolution of two `u64` sequences, provided every coefficient of the true result
/// is below `2^64` (checked modulo three primes; the result is otherwise garbage). The inputs
/// themselves must be below `2^64` but are reduced internally.
///
/// # Panics
/// Panics if the result is longer than `2^23`.
pub fn convolution_u64(a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let reduce = |v: &[u64], m: u32| -> Vec<u32> { v.iter().map(|&x| (x % m as u64) as u32).collect() };
    let c1 = convolution_mod::<P1>(&reduce(a, P1), &reduce(b, P1));
    let c2 = convolution_mod::<P2>(&reduce(a, P2), &reduce(b, P2));
    let c3 = convolution_mod::<P3>(&reduce(a, P3), &reduce(b, P3));
    // Garner's algorithm for x = c1 + P1 * t1 + P1 * P2 * t2, wrapping modulo 2^64.
    let (m1, m2, m3) = (P1 as u64, P2 as u64, P3 as u64);
    let inv_m1_mod_m2 = pow_mod(m1, m2 - 2, m2);
    let inv_m1m2_mod_m3 = pow_mod(m1 * m2 % m3, m3 - 2, m3);
    let m1m2 = m1.wrapping_mul(m2);
    (0..c1.len())
        .map(|k| {
            let (r1, r2, r3) = (c1[k] as u64, c2[k] as u64, c3[k] as u64);
            let t1 = (r2 + m2 - r1 % m2) % m2 * inv_m1_mod_m2 % m2;
            let x12 = r1 + m1 * t1; // < m1 * m2 < 2^64
            let t2 = (r3 + m3 - x12 % m3) % m3 * inv_m1m2_mod_m3 % m3;
            x12.wrapping_add(m1m2.wrapping_mul(t2))
        })
        .collect()
}

/// `c[k] = sum_{i + j = k} a[i] * b[j] mod 2^64` (wrapping arithmetic).
///
/// Each number is cut into four 16-bit limbs. For every output limb position `s = 0..3` the
/// products `A_i * B_j` with `i + j = s` are summed in the transform domain; the sums stay below
/// `2^54` for lengths up to `2^21`, so two NTT primes recover them exactly. Twelve transforms
/// per prime (eight forward, four inverse).
///
/// # Panics
/// Panics if the result is longer than `2^21` elements (the exactness bound) .
pub fn convolution_mod_2_64(a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let len = a.len() + b.len() - 1;
    if a.len().min(b.len()) <= 16 {
        let mut c = vec![0u64; len];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                c[i + j] = c[i + j].wrapping_add(x.wrapping_mul(y));
            }
        }
        return c;
    }
    assert!(a.len().min(b.len()) <= 1 << 21, "sequences too long for the 2^54 exactness bound");
    let n = len.next_power_of_two();
    // residues of the four sums for each prime
    fn per_prime<const P: u32>(a: &[u64], b: &[u64], n: usize) -> [Vec<u32>; 4] {
        let limbs = |v: &[u64]| -> [Vec<u32>; 4] {
            let mut out = [vec![0u32; n], vec![0u32; n], vec![0u32; n], vec![0u32; n]];
            for (i, &x) in v.iter().enumerate() {
                for (k, o) in out.iter_mut().enumerate() {
                    o[i] = ((x >> (16 * k)) & 0xFFFF) as u32;
                }
            }
            for o in out.iter_mut() {
                ntt::ntt::<P>(o);
            }
            out
        };
        let (fa, fb) = (limbs(a), limbs(b));
        let mut res = [vec![0u32; n], vec![0u32; n], vec![0u32; n], vec![0u32; n]];
        for (s, r) in res.iter_mut().enumerate() {
            for t in 0..n {
                let mut acc = 0u64;
                for i in 0..=s {
                    acc += fa[i][t] as u64 * fb[s - i][t] as u64 % P as u64;
                }
                r[t] = (acc % P as u64) as u32;
            }
            ntt::intt::<P>(r);
        }
        res
    }
    let r1 = per_prime::<P1>(a, b, n);
    let r2 = per_prime::<P2>(a, b, n);
    let inv = pow_mod(P1 as u64, P2 as u64 - 2, P2 as u64);
    (0..len)
        .map(|k| {
            let mut total = 0u64;
            for s in 0..4 {
                let (x1, x2) = (r1[s][k] as u64, r2[s][k] as u64);
                let t = (x2 + P2 as u64 - x1 % P2 as u64) % P2 as u64 * inv % P2 as u64;
                total = total.wrapping_add((x1 + P1 as u64 * t) << (16 * s));
            }
            total
        })
        .collect()
}
