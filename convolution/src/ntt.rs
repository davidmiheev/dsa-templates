//! Radix-4 number-theoretic transforms (no bit reversal): forward = decimation in frequency,
//! inverse = decimation in time, so `intt(ntt(a)) == n * a` and pointwise products in between
//! are convolutions. The twiddle factors are generated on the fly from `O(log n)` constants.

#[inline(always)]
pub(crate) fn mul<const P: u32>(a: u32, b: u32) -> u32 {
    (a as u64 * b as u64 % P as u64) as u32
}

/// `-P^{-1} mod 2^32`, for Montgomery reduction.
pub(crate) const fn neg_inv(p: u32) -> u32 {
    let mut x = p; // correct to 3 bits, each step doubles the number of correct bits
    let mut i = 0;
    while i < 5 {
        x = x.wrapping_mul(2u32.wrapping_sub(p.wrapping_mul(x)));
        i += 1;
    }
    x.wrapping_neg()
}

/// Montgomery product `a * b / 2^32 mod P` of residues `< P`, result in `[0, P)`.
#[inline(always)]
pub(crate) fn mmul<const P: u32>(a: u32, b: u32) -> u32 {
    let t = a as u64 * b as u64;
    let m = (t as u32).wrapping_mul(neg_inv(P));
    let u = ((t + m as u64 * P as u64) >> 32) as u32;
    if u >= P {
        u - P
    } else {
        u
    }
}

/// `x * 2^32 mod P`: the Montgomery form of `x`.
pub(crate) fn to_mont<const P: u32>(x: u32) -> u32 {
    (((x as u64) << 32) % P as u64) as u32
}

#[inline(always)]
pub(crate) fn add<const P: u32>(a: u32, b: u32) -> u32 {
    let s = a + b;
    if s >= P {
        s - P
    } else {
        s
    }
}

#[inline(always)]
pub(crate) fn sub<const P: u32>(a: u32, b: u32) -> u32 {
    if a >= b {
        a - b
    } else {
        a + P - b
    }
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

/// Constants of one prime; the rate tables and `imag` are in Montgomery form.
pub(crate) struct Info {
    pub root: [u32; 32],
    pub iroot: [u32; 32],
    pub rate2: [u32; 32],
    pub irate2: [u32; 32],
    pub rate3: [u32; 32],
    pub irate3: [u32; 32],
}

pub(crate) fn info<const P: u32>() -> Info {
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
        rate2[i] = to_mont::<P>(mul::<P>(root[i + 2], prod));
        irate2[i] = to_mont::<P>(mul::<P>(iroot[i + 2], iprod));
        prod = mul::<P>(prod, iroot[i + 2]);
        iprod = mul::<P>(iprod, root[i + 2]);
    }
    let (mut prod, mut iprod) = (1u32, 1u32);
    for i in 0..r2.saturating_sub(2) {
        rate3[i] = to_mont::<P>(mul::<P>(root[i + 3], prod));
        irate3[i] = to_mont::<P>(mul::<P>(iroot[i + 3], iprod));
        prod = mul::<P>(prod, iroot[i + 3]);
        iprod = mul::<P>(iprod, root[i + 3]);
    }
    Info { root, iroot, rate2, irate2, rate3, irate3 }
}

/// In-place forward transform of a power-of-two length slice (output in bit-reversed-like order).
/// Runs the AVX2 version when the CPU has it (identical output).
pub fn ntt<const P: u32>(a: &mut [u32]) {
    #[cfg(target_arch = "x86_64")]
    if a.len() >= 32 && crate::simd::available() {
        return crate::simd::ntt_avx2::<P>(a);
    }
    ntt_scalar::<P>(a)
}

/// The portable radix-4 forward transform (see [`ntt`]).
pub fn ntt_scalar<const P: u32>(a: &mut [u32]) {
    let n = a.len();
    assert!(n.is_power_of_two(), "length {n} is not a power of two");
    let h = n.trailing_zeros() as usize;
    assert!(h <= rank2::<P>() as usize, "length 2^{h} exceeds the 2-adic order of {P} - 1");
    let t = info::<P>();
    let one_m = to_mont::<P>(1);
    let imag = to_mont::<P>(t.root[2]);
    let mut len = 0;
    while len < h {
        if h - len == 1 {
            let p = 1 << (h - len - 1);
            let mut rot = one_m;
            for s in 0..1usize << len {
                let offset = s << (h - len);
                let (lo, hi) = a[offset..offset + 2 * p].split_at_mut(p);
                for (l, r) in lo.iter_mut().zip(hi.iter_mut()) {
                    let left = *l;
                    let right = mmul::<P>(*r, rot);
                    *l = add::<P>(left, right);
                    *r = sub::<P>(left, right);
                }
                if s + 1 != 1 << len {
                    rot = mmul::<P>(rot, t.rate2[(!s).trailing_zeros() as usize]);
                }
            }
            len += 1;
        } else {
            let p = 1 << (h - len - 2);
            let mut rot = one_m;
            for s in 0..1usize << len {
                let rot2 = mmul::<P>(rot, rot);
                let rot3 = mmul::<P>(rot2, rot);
                let offset = s << (h - len);
                let block = &mut a[offset..offset + 4 * p];
                let (b01, b23) = block.split_at_mut(2 * p);
                let (b0, b1) = b01.split_at_mut(p);
                let (b2, b3) = b23.split_at_mut(p);
                for i in 0..p {
                    let a0 = b0[i];
                    let a1 = mmul::<P>(b1[i], rot);
                    let a2 = mmul::<P>(b2[i], rot2);
                    let a3 = mmul::<P>(b3[i], rot3);
                    let s02 = add::<P>(a0, a2);
                    let d02 = sub::<P>(a0, a2);
                    let s13 = add::<P>(a1, a3);
                    let d13 = mmul::<P>(sub::<P>(a1, a3), imag);
                    b0[i] = add::<P>(s02, s13);
                    b1[i] = sub::<P>(s02, s13);
                    b2[i] = add::<P>(d02, d13);
                    b3[i] = sub::<P>(d02, d13);
                }
                if s + 1 != 1 << len {
                    rot = mmul::<P>(rot, t.rate3[(!s).trailing_zeros() as usize]);
                }
            }
            len += 2;
        }
    }
}

/// In-place inverse transform (input as produced by [`ntt`]); the result is multiplied by `n`
/// compared to the true inverse, so divide by `n` (see [`intt`]).
pub fn intt_unscaled<const P: u32>(a: &mut [u32]) {
    #[cfg(target_arch = "x86_64")]
    if a.len() >= 32 && crate::simd::available() {
        return crate::simd::intt_unscaled_avx2::<P>(a);
    }
    intt_unscaled_scalar::<P>(a)
}

/// The portable radix-4 inverse transform (see [`intt_unscaled`]).
pub fn intt_unscaled_scalar<const P: u32>(a: &mut [u32]) {
    let n = a.len();
    assert!(n.is_power_of_two(), "length {n} is not a power of two");
    let h = n.trailing_zeros() as usize;
    assert!(h <= rank2::<P>() as usize, "length 2^{h} exceeds the 2-adic order of {P} - 1");
    let t = info::<P>();
    let one_m = to_mont::<P>(1);
    let iimag = to_mont::<P>(t.iroot[2]);
    let mut len = h;
    while len > 0 {
        if len == 1 {
            let p = 1 << (h - len);
            let mut irot = one_m;
            for s in 0..1usize << (len - 1) {
                let offset = s << (h - len + 1);
                let (lo, hi) = a[offset..offset + 2 * p].split_at_mut(p);
                for (l, r) in lo.iter_mut().zip(hi.iter_mut()) {
                    let left = *l;
                    let right = *r;
                    *l = add::<P>(left, right);
                    *r = mmul::<P>(sub::<P>(left, right), irot);
                }
                if s + 1 != 1 << (len - 1) {
                    irot = mmul::<P>(irot, t.irate2[(!s).trailing_zeros() as usize]);
                }
            }
            len -= 1;
        } else {
            let p = 1 << (h - len);
            let mut irot = one_m;
            for s in 0..1usize << (len - 2) {
                let irot2 = mmul::<P>(irot, irot);
                let irot3 = mmul::<P>(irot2, irot);
                let offset = s << (h - len + 2);
                let block = &mut a[offset..offset + 4 * p];
                let (b01, b23) = block.split_at_mut(2 * p);
                let (b0, b1) = b01.split_at_mut(p);
                let (b2, b3) = b23.split_at_mut(p);
                for i in 0..p {
                    let (a0, a1, a2, a3) = (b0[i], b1[i], b2[i], b3[i]);
                    let s01 = add::<P>(a0, a1);
                    let d01 = sub::<P>(a0, a1);
                    let s23 = add::<P>(a2, a3);
                    let d23 = mmul::<P>(sub::<P>(a2, a3), iimag);
                    b0[i] = add::<P>(s01, s23);
                    b1[i] = mmul::<P>(add::<P>(d01, d23), irot);
                    b2[i] = mmul::<P>(sub::<P>(s01, s23), irot2);
                    b3[i] = mmul::<P>(sub::<P>(d01, d23), irot3);
                }
                if s + 1 != 1 << (len - 2) {
                    irot = mmul::<P>(irot, t.irate3[(!s).trailing_zeros() as usize]);
                }
            }
            len -= 2;
        }
    }
}

/// In-place inverse transform including the `1 / n` scaling.
pub fn intt<const P: u32>(a: &mut [u32]) {
    intt_unscaled::<P>(a);
    let n_inv = to_mont::<P>(pow_mod(a.len() as u64, P as u64 - 2, P as u64) as u32);
    for x in a.iter_mut() {
        *x = mmul::<P>(*x, n_inv);
    }
}

/// `a[i] = a[i] * b[i] / n` for the transforms of two length-`n` polynomials: the pointwise
/// product with the `1 / n` of the inverse transform folded in, so that
/// [`intt_unscaled`] of the result is the cyclic convolution.
pub fn pointwise_scaled<const P: u32>(a: &mut [u32], b: &[u32]) {
    assert_eq!(a.len(), b.len(), "spectra of different lengths");
    let n = a.len() as u64;
    let n_inv = pow_mod(n, P as u64 - 2, P as u64) as u32;
    #[cfg(target_arch = "x86_64")]
    if a.len() >= 32 && crate::simd::available() {
        // two Montgomery products: a * b / 2^32, then times (n_inv * 2^64) / 2^32
        return crate::simd::pointwise_avx2::<P>(a, b, to_mont::<P>(to_mont::<P>(n_inv)));
    }
    for (x, &y) in a.iter_mut().zip(b) {
        *x = mul::<P>(mul::<P>(*x, y), n_inv);
    }
}
