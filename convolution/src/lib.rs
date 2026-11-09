//! Convolutions by the number-theoretic transform (NTT).
//!
//! * [`convolution`]: `c[k] = sum a[i] * b[k - i] mod 998244353`.
//! * [`convolution_mod`]: the same for the primes [`P1`], [`P2`], [`P3`] (selected by a const
//!   generic), when the lengths allow (`len(a) + len(b) - 1 <= 2^23`, `2^25`, `2^24`).
//! * [`convolution_u64`]: the *exact* result when every coefficient of the true convolution is
//!   below `2^64` (computed modulo three primes and recombined), e.g. for counting problems.
//!
//! Small inputs fall back to the schoolbook product. Transforms are iterative radix-2 without
//! bit reversal (decimation in frequency forward, in time inverse).
//!
//! # Example
//!
//! ```
//! use convolution::{convolution, convolution_u64};
//!
//! assert_eq!(convolution(&[1, 2, 3], &[4, 5]), vec![4, 13, 22, 15]);
//! assert_eq!(convolution_u64(&[1 << 40, 1], &[1 << 20, 1]), vec![1 << 60, (1 << 40) + (1 << 20), 1]);
//! ```

/// 998244353 = 119 * 2^23 + 1, primitive root 3.
pub const P1: u32 = 998_244_353;
/// 167772161 = 5 * 2^25 + 1, primitive root 3.
pub const P2: u32 = 167_772_161;
/// 469762049 = 7 * 2^26 + 1, primitive root 3.
pub const P3: u32 = 469_762_049;

const fn pow_mod(mut b: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1;
    b %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % m;
        }
        b = b * b % m;
        e >>= 1;
    }
    r
}

/// Twiddle tables for transforms of size `n` (a power of two): `fw[half + j] = w_{2 half}^j`.
struct Tables {
    fw: Vec<u32>,
    inv: Vec<u32>,
}

impl Tables {
    fn new<const M: u32>(n: usize) -> Self {
        let m = M as u64;
        let (mut fw, mut inv) = (vec![1u32; n.max(2)], vec![1u32; n.max(2)]);
        let mut half = 1;
        while half < n {
            let w = pow_mod(3, (m - 1) / (2 * half as u64), m);
            let wi = pow_mod(w, m - 2, m);
            let (mut c, mut ci) = (1u64, 1u64);
            for j in 0..half {
                fw[half + j] = c as u32;
                inv[half + j] = ci as u32;
                c = c * w % m;
                ci = ci * wi % m;
            }
            half *= 2;
        }
        Tables { fw, inv }
    }
}

/// In-place forward NTT (decimation in frequency): natural order in, bit-reversed out.
fn ntt<const M: u32>(a: &mut [u32], t: &Tables) {
    let n = a.len();
    let mut half = n / 2;
    while half >= 1 {
        for block in a.chunks_exact_mut(2 * half) {
            let (lo, hi) = block.split_at_mut(half);
            for j in 0..half {
                let (u, v) = (lo[j], hi[j]);
                let s = u + v;
                lo[j] = if s >= M { s - M } else { s };
                let d = if u >= v { u - v } else { u + M - v };
                hi[j] = (d as u64 * t.fw[half + j] as u64 % M as u64) as u32;
            }
        }
        half /= 2;
    }
}

/// In-place inverse NTT (decimation in time): bit-reversed in, natural order out, scaled.
fn intt<const M: u32>(a: &mut [u32], t: &Tables) {
    let n = a.len();
    let mut half = 1;
    while half < n {
        for block in a.chunks_exact_mut(2 * half) {
            let (lo, hi) = block.split_at_mut(half);
            for j in 0..half {
                let u = lo[j];
                let v = (hi[j] as u64 * t.inv[half + j] as u64 % M as u64) as u32;
                let s = u + v;
                lo[j] = if s >= M { s - M } else { s };
                hi[j] = if u >= v { u - v } else { u + M - v };
            }
        }
        half *= 2;
    }
    let n_inv = pow_mod(n as u64, M as u64 - 2, M as u64);
    for x in a.iter_mut() {
        *x = (*x as u64 * n_inv % M as u64) as u32;
    }
}

fn naive<const M: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut c = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            c[i + j] = (c[i + j] + x as u64 * y as u64) % M as u64;
        }
    }
    c.into_iter().map(|x| x as u32).collect()
}

/// Convolution modulo the NTT-friendly prime `M` (one of [`P1`], [`P2`], [`P3`]). The inputs
/// must be reduced (`< M`).
///
/// # Panics
/// Panics if `M` is not one of the three supported primes, or the result length exceeds the
/// largest power of two dividing `M - 1`.
pub fn convolution_mod<const M: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert!(M == P1 || M == P2 || M == P3, "unsupported modulus {M}");
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let len = a.len() + b.len() - 1;
    if a.len().min(b.len()) <= 24 {
        return naive::<M>(a, b);
    }
    let n = len.next_power_of_two();
    assert!(n as u64 <= 1 << (M - 1).trailing_zeros(), "result too long for this modulus");
    let t = Tables::new::<M>(n);
    let (mut fa, mut fb) = (vec![0u32; n], vec![0u32; n]);
    fa[..a.len()].copy_from_slice(a);
    fb[..b.len()].copy_from_slice(b);
    ntt::<M>(&mut fa, &t);
    ntt::<M>(&mut fb, &t);
    for (x, &y) in fa.iter_mut().zip(&fb) {
        *x = (*x as u64 * y as u64 % M as u64) as u32;
    }
    intt::<M>(&mut fa, &t);
    fa.truncate(len);
    fa
}

/// `c[k] = sum_{i + j = k} a[i] * b[j] mod 998244353`. Inputs must be `< 998244353`.
pub fn convolution(a: &[u32], b: &[u32]) -> Vec<u32> {
    convolution_mod::<P1>(a, b)
}

/// The exact convolution of two `u64` sequences, provided every coefficient of the true result
/// is below `2^64` (checked modulo three primes; the result is otherwise garbage). The inputs
/// themselves must be below `2^64` but are reduced internally.
///
/// # Panics
/// Panics if the result is longer than `2^24`.
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
