//! Radix-4 number-theoretic transforms (no bit reversal): forward = decimation in frequency,
//! inverse = decimation in time, so `intt(ntt(a)) == n * a` and pointwise products in between
//! are convolutions. The twiddle factors are generated on the fly from `O(log n)` constants.

#[inline(always)]
fn mul<const P: u32>(a: u32, b: u32) -> u32 {
    (a as u64 * b as u64 % P as u64) as u32
}

pub(crate) const fn pow_mod(mut b: u64, mut e: u64, m: u64) -> u64 {
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

/// A primitive root of the supported NTT-friendly primes.
///
/// # Panics
/// Panics for any other modulus.
pub const fn primitive_root(p: u32) -> u32 {
    match p {
        998_244_353 | 167_772_161 | 469_762_049 | 1_004_535_809 => 3,
        2_013_265_921 => 31,
        754_974_721 => 11,
        924_844_033 => 5,
        _ => panic!("unsupported NTT modulus"),
    }
}

/// The 2-adic order of `P - 1`: the largest transform size is `2^rank2`.
pub(crate) const fn rank2<const P: u32>() -> u32 {
    (P - 1).trailing_zeros()
}

/// Constants of one prime.
struct Info {
    root: [u32; 32],
    iroot: [u32; 32],
    rate2: [u32; 32],
    irate2: [u32; 32],
    rate3: [u32; 32],
    irate3: [u32; 32],
}

fn info<const P: u32>() -> Info {
    let r2 = rank2::<P>() as usize;
    let m = P as u64;
    let mut root = [0u32; 32];
    let mut iroot = [0u32; 32];
    root[r2] = pow_mod(primitive_root(P) as u64, (m - 1) >> r2, m) as u32;
    iroot[r2] = pow_mod(root[r2] as u64, m - 2, m) as u32;
    for i in (0..r2).rev() {
        root[i] = mul::<P>(root[i + 1], root[i + 1]);
        iroot[i] = mul::<P>(iroot[i + 1], iroot[i + 1]);
    }
    let (mut rate2, mut irate2, mut rate3, mut irate3) = ([0u32; 32], [0u32; 32], [0u32; 32], [0u32; 32]);
    let (mut prod, mut iprod) = (1u32, 1u32);
    for i in 0..r2.saturating_sub(1) {
        rate2[i] = mul::<P>(root[i + 2], prod);
        irate2[i] = mul::<P>(iroot[i + 2], iprod);
        prod = mul::<P>(prod, iroot[i + 2]);
        iprod = mul::<P>(iprod, root[i + 2]);
    }
    let (mut prod, mut iprod) = (1u32, 1u32);
    for i in 0..r2.saturating_sub(2) {
        rate3[i] = mul::<P>(root[i + 3], prod);
        irate3[i] = mul::<P>(iroot[i + 3], iprod);
        prod = mul::<P>(prod, iroot[i + 3]);
        iprod = mul::<P>(iprod, root[i + 3]);
    }
    Info { root, iroot, rate2, irate2, rate3, irate3 }
}

/// In-place forward transform of a power-of-two length slice (output in bit-reversed-like order).
pub fn ntt<const P: u32>(a: &mut [u32]) {
    let n = a.len();
    assert!(n.is_power_of_two(), "length {n} is not a power of two");
    let h = n.trailing_zeros() as usize;
    assert!(h <= rank2::<P>() as usize, "length 2^{h} exceeds the 2-adic order of {P} - 1");
    let t = info::<P>();
    let imag = t.root[2];
    let mut len = 0;
    while len < h {
        if h - len == 1 {
            let p = 1 << (h - len - 1);
            let mut rot = 1u32;
            for s in 0..1usize << len {
                let offset = s << (h - len);
                let (lo, hi) = a[offset..offset + 2 * p].split_at_mut(p);
                for (l, r) in lo.iter_mut().zip(hi.iter_mut()) {
                    let left = *l;
                    let right = mul::<P>(*r, rot);
                    *l = if left + right >= P { left + right - P } else { left + right };
                    *r = if left >= right { left - right } else { left + P - right };
                }
                if s + 1 != 1 << len {
                    rot = mul::<P>(rot, t.rate2[(!s).trailing_zeros() as usize]);
                }
            }
            len += 1;
        } else {
            let p = 1 << (h - len - 2);
            let mut rot = 1u32;
            for s in 0..1usize << len {
                let rot2 = mul::<P>(rot, rot);
                let rot3 = mul::<P>(rot2, rot);
                let offset = s << (h - len);
                let block = &mut a[offset..offset + 4 * p];
                let (b01, b23) = block.split_at_mut(2 * p);
                let (b0, b1) = b01.split_at_mut(p);
                let (b2, b3) = b23.split_at_mut(p);
                for i in 0..p {
                    let a0 = b0[i] as u64;
                    let a1 = mul::<P>(b1[i], rot) as u64;
                    let a2 = mul::<P>(b2[i], rot2) as u64;
                    let a3 = mul::<P>(b3[i], rot3) as u64;
                    let a1na3imag = (a1 + P as u64 - a3) * imag as u64 % P as u64;
                    let a0pa2 = a0 + a2;
                    let a0na2 = a0 + P as u64 - a2;
                    b0[i] = ((a0pa2 + a1 + a3) % P as u64) as u32;
                    b1[i] = ((a0pa2 + 2 * P as u64 - (a1 + a3)) % P as u64) as u32;
                    b2[i] = ((a0na2 + a1na3imag) % P as u64) as u32;
                    b3[i] = ((a0na2 + P as u64 - a1na3imag) % P as u64) as u32;
                }
                if s + 1 != 1 << len {
                    rot = mul::<P>(rot, t.rate3[(!s).trailing_zeros() as usize]);
                }
            }
            len += 2;
        }
    }
}

/// In-place inverse transform (input as produced by [`ntt`]); the result is multiplied by `n`
/// compared to the true inverse, so divide by `n` (see [`intt`]).
pub fn intt_unscaled<const P: u32>(a: &mut [u32]) {
    let n = a.len();
    assert!(n.is_power_of_two(), "length {n} is not a power of two");
    let h = n.trailing_zeros() as usize;
    assert!(h <= rank2::<P>() as usize, "length 2^{h} exceeds the 2-adic order of {P} - 1");
    let t = info::<P>();
    let iimag = t.iroot[2];
    let mut len = h;
    while len > 0 {
        if len == 1 {
            let p = 1 << (h - len);
            let mut irot = 1u32;
            for s in 0..1usize << (len - 1) {
                let offset = s << (h - len + 1);
                let (lo, hi) = a[offset..offset + 2 * p].split_at_mut(p);
                for (l, r) in lo.iter_mut().zip(hi.iter_mut()) {
                    let left = *l;
                    let right = *r;
                    *l = if left + right >= P { left + right - P } else { left + right };
                    let d = if left >= right { left - right } else { left + P - right };
                    *r = mul::<P>(d, irot);
                }
                if s + 1 != 1 << (len - 1) {
                    irot = mul::<P>(irot, t.irate2[(!s).trailing_zeros() as usize]);
                }
            }
            len -= 1;
        } else {
            let p = 1 << (h - len);
            let mut irot = 1u32;
            for s in 0..1usize << (len - 2) {
                let irot2 = mul::<P>(irot, irot);
                let irot3 = mul::<P>(irot2, irot);
                let offset = s << (h - len + 2);
                let block = &mut a[offset..offset + 4 * p];
                let (b01, b23) = block.split_at_mut(2 * p);
                let (b0, b1) = b01.split_at_mut(p);
                let (b2, b3) = b23.split_at_mut(p);
                for i in 0..p {
                    let a0 = b0[i] as u64;
                    let a1 = b1[i] as u64;
                    let a2 = b2[i] as u64;
                    let a3 = b3[i] as u64;
                    let a2na3iimag = (a2 + P as u64 - a3) * iimag as u64 % P as u64;
                    b0[i] = ((a0 + a1 + a2 + a3) % P as u64) as u32;
                    b1[i] = mul::<P>(((a0 + P as u64 - a1 + a2na3iimag) % P as u64) as u32, irot);
                    b2[i] = mul::<P>(((a0 + a1 + 2 * P as u64 - a2 - a3) % P as u64) as u32, irot2);
                    b3[i] = mul::<P>(((a0 + P as u64 - a1 + P as u64 - a2na3iimag) % P as u64) as u32, irot3);
                }
                if s + 1 != 1 << (len - 2) {
                    irot = mul::<P>(irot, t.irate3[(!s).trailing_zeros() as usize]);
                }
            }
            len -= 2;
        }
    }
}

/// In-place inverse transform including the `1 / n` scaling.
pub fn intt<const P: u32>(a: &mut [u32]) {
    intt_unscaled::<P>(a);
    let n_inv = pow_mod(a.len() as u64, P as u64 - 2, P as u64) as u32;
    for x in a.iter_mut() {
        *x = mul::<P>(*x, n_inv);
    }
}
