//! Experimental AVX-512 version of the transforms of [`crate::simd`] (sixteen Montgomery lanes).
//!
//! Needs the `avx512` cargo feature (the AVX-512 intrinsics are stable from Rust 1.89) and the
//! `avx512f` CPU flag at run time. Same layer structure and the same results as the other
//! engines: the layers with span >= 32 are vectorised along the span, the last five layers
//! (spans 16, 8, 4, 2, 1) run inside two registers per 32 elements with `vpermt2d` shuffles.

use crate::ntt::{info, mmul, neg_inv, rank2, to_mont};
use std::arch::x86_64::*;

/// Whether the AVX-512 transforms can run on this CPU (and were not switched off with
/// [`crate::simd::set_enabled`]).
pub fn available() -> bool {
    crate::simd::enabled() && std::is_x86_feature_detected!("avx512f")
}

struct K {
    p: __m512i,
}

#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn konst<const P: u32>() -> K {
    K { p: _mm512_set1_epi32(P as i32) }
}

#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn bc(x: u32) -> __m512i {
    _mm512_set1_epi32(x as i32)
}

#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn ld(p: *const u32) -> __m512i {
    _mm512_loadu_si512(p as *const _)
}

#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn st(p: *mut u32, v: __m512i) {
    _mm512_storeu_si512(p as *mut _, v)
}

struct W {
    w: __m512i,
    wn: __m512i,
}

#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn mkw<const P: u32>(w: u32) -> W {
    W { w: bc(w), wn: bc(w.wrapping_mul(neg_inv(P))) }
}

/// `a * w / 2^32`, broadcast twiddle with its precomputed quotient multiplier.
#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn vmulc(a: __m512i, w: &W, k: &K) -> __m512i {
    let a1 = _mm512_srli_epi64(a, 32);
    let t0 = _mm512_mul_epu32(a, w.w);
    let t1 = _mm512_mul_epu32(a1, w.w);
    let m0 = _mm512_mul_epu32(a, w.wn);
    let m1 = _mm512_mul_epu32(a1, w.wn);
    let u0 = _mm512_add_epi64(t0, _mm512_mul_epu32(m0, k.p));
    let u1 = _mm512_add_epi64(t1, _mm512_mul_epu32(m1, k.p));
    let r = _mm512_mask_blend_epi32(0xAAAA, _mm512_srli_epi64(u0, 32), u1);
    _mm512_min_epu32(r, _mm512_sub_epi32(r, k.p))
}

/// `a * b / 2^32` for vector twiddles.
#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn vmul<const P: u32>(a: __m512i, b: __m512i, k: &K) -> __m512i {
    let ninv = bc(neg_inv(P));
    let a1 = _mm512_srli_epi64(a, 32);
    let b1 = _mm512_srli_epi64(b, 32);
    let t0 = _mm512_mul_epu32(a, b);
    let t1 = _mm512_mul_epu32(a1, b1);
    let m0 = _mm512_mul_epu32(t0, ninv);
    let m1 = _mm512_mul_epu32(t1, ninv);
    let u0 = _mm512_add_epi64(t0, _mm512_mul_epu32(m0, k.p));
    let u1 = _mm512_add_epi64(t1, _mm512_mul_epu32(m1, k.p));
    let r = _mm512_mask_blend_epi32(0xAAAA, _mm512_srli_epi64(u0, 32), u1);
    _mm512_min_epu32(r, _mm512_sub_epi32(r, k.p))
}

#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn vadd(a: __m512i, b: __m512i, k: &K) -> __m512i {
    let s = _mm512_add_epi32(a, b);
    _mm512_min_epu32(s, _mm512_sub_epi32(s, k.p))
}

#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn vsub(a: __m512i, b: __m512i, k: &K) -> __m512i {
    let d = _mm512_sub_epi32(a, b);
    _mm512_min_epu32(d, _mm512_add_epi32(d, k.p))
}

fn idx(f: impl Fn(usize) -> usize) -> [i32; 16] {
    let mut r = [0i32; 16];
    for (k, x) in r.iter_mut().enumerate() {
        *x = f(k) as i32;
    }
    r
}

/// The shuffle tables of the five in-register layers.
struct Perm {
    // gather L (low half of every pair) and H (high half) of span 8 / 4 / 2 / 1 from (x0, x1)
    l: [__m512i; 4],
    h: [__m512i; 4],
    // twiddle selectors for spans 8 / 4 / 2 (span 1 loads the table directly)
    tw: [__m512i; 3],
    // rebuild x0, x1 from (aa, bb) for spans 8 / 4 / 2 / 1
    o0: [__m512i; 4],
    o1: [__m512i; 4],
}

#[target_feature(enable = "avx512f")]
unsafe fn perms() -> Perm {
    let v = |a: [i32; 16]| _mm512_loadu_si512(a.as_ptr() as *const _);
    // lane j of the 32 elements of (x0, x1) is j for j < 16 and 16 + (j - 16) otherwise, so a
    // "source lane" below 16 reads the first operand and 16.. the second
    let mut l = [v([0; 16]); 4];
    let mut h = [v([0; 16]); 4];
    let mut o0 = [v([0; 16]); 4];
    let mut o1 = [v([0; 16]); 4];
    for (t, s) in [8usize, 4, 2, 1].into_iter().enumerate() {
        // elements of the pair-low / pair-high positions in the 32-element chunk
        let lows: Vec<usize> = (0..32).filter(|j| j & s == 0).collect();
        let highs: Vec<usize> = (0..32).filter(|j| j & s != 0).collect();
        l[t] = v(idx(|k| lows[k]));
        h[t] = v(idx(|k| highs[k]));
        // after the butterfly aa holds the new pair-low elements (in `lows` order) and bb the
        // pair-high ones; element e of the chunk is aa[pos in lows] or bb[pos in highs] + 16
        let src = |e: usize| -> usize {
            if e & s == 0 {
                lows.iter().position(|&x| x == e).unwrap()
            } else {
                16 + highs.iter().position(|&x| x == e).unwrap()
            }
        };
        o0[t] = v(idx(src));
        o1[t] = v(idx(|k| src(16 + k)));
    }
    let tw = [
        v(idx(|k| k / 8)),  // span 8: two blocks of 16 -> [r0 x8 | r1 x8]
        v(idx(|k| k / 4)),  // span 4: four blocks
        v(idx(|k| k / 2)),  // span 2: eight blocks
    ];
    Perm { l, h, tw, o0, o1 }
}

/// `x0, x1 <- the 5 layers` helper: gather, butterfly with the twiddle vector, scatter.
#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn fwd_layer<const P: u32>(x0: &mut __m512i, x1: &mut __m512i, t: usize, tw: __m512i, p: &Perm, k: &K) {
    let l = _mm512_permutex2var_epi32(*x0, p.l[t], *x1);
    let h = _mm512_permutex2var_epi32(*x0, p.h[t], *x1);
    let tt = vmul::<P>(h, tw, k);
    let (aa, bb) = (vadd(l, tt, k), vsub(l, tt, k));
    *x0 = _mm512_permutex2var_epi32(aa, p.o0[t], bb);
    *x1 = _mm512_permutex2var_epi32(aa, p.o1[t], bb);
}

#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn inv_layer<const P: u32>(x0: &mut __m512i, x1: &mut __m512i, t: usize, tw: __m512i, p: &Perm, k: &K) {
    let l = _mm512_permutex2var_epi32(*x0, p.l[t], *x1);
    let h = _mm512_permutex2var_epi32(*x0, p.h[t], *x1);
    let (aa, bb) = (vadd(l, h, k), vmul::<P>(vsub(l, h, k), tw, k));
    *x0 = _mm512_permutex2var_epi32(aa, p.o0[t], bb);
    *x1 = _mm512_permutex2var_epi32(aa, p.o1[t], bb);
}

fn plan(h: usize) -> Vec<(usize, bool)> {
    let mut stages = Vec::new();
    let mut len = 0;
    while len + 5 < h {
        if h - 5 - len == 1 {
            stages.push((len, false));
            len += 1;
        } else {
            stages.push((len, true));
            len += 2;
        }
    }
    stages
}

/// Forward transform, same output as [`crate::ntt::ntt_scalar`].
///
/// # Panics
/// Panics if AVX-512 is unavailable, the length is not a power of two or it exceeds the 2-adic order.
pub fn ntt_avx512<const P: u32>(a: &mut [u32]) {
    assert!(available(), "AVX-512 is not available");
    let n = a.len();
    assert!(n.is_power_of_two(), "length {n} is not a power of two");
    let h = n.trailing_zeros() as usize;
    assert!(h <= rank2::<P>() as usize, "length 2^{h} exceeds the 2-adic order of {P} - 1");
    if h < 6 {
        return crate::ntt::ntt_scalar::<P>(a);
    }
    crate::simd::with_tables_pub::<P, _>(n / 2, |r, _ir| unsafe { ntt_impl::<P>(a, r) })
}

#[target_feature(enable = "avx512f")]
unsafe fn ntt_impl<const P: u32>(a: &mut [u32], r: &[u32]) {
    let n = a.len();
    let h = n.trailing_zeros() as usize;
    let t = info::<P>();
    let k = konst::<P>();
    let one_m = to_mont::<P>(1);
    let imag = mkw::<P>(to_mont::<P>(t.root[2]));
    let ptr = a.as_mut_ptr();
    for (len, is4) in plan(h) {
        let mut rot = one_m;
        if is4 {
            let p = 1usize << (h - len - 2);
            for s in 0..1usize << len {
                let (rot2, rot3) = (mmul::<P>(rot, rot), mmul::<P>(mmul::<P>(rot, rot), rot));
                let (r1, r2, r3) = (mkw::<P>(rot), mkw::<P>(rot2), mkw::<P>(rot3));
                let b = ptr.add(s * 4 * p);
                let mut i = 0;
                while i < p {
                    let a0 = ld(b.add(i));
                    let a1 = vmulc(ld(b.add(p + i)), &r1, &k);
                    let a2 = vmulc(ld(b.add(2 * p + i)), &r2, &k);
                    let a3 = vmulc(ld(b.add(3 * p + i)), &r3, &k);
                    let s02 = vadd(a0, a2, &k);
                    let d02 = vsub(a0, a2, &k);
                    let s13 = vadd(a1, a3, &k);
                    let d13 = vmulc(vsub(a1, a3, &k), &imag, &k);
                    st(b.add(i), vadd(s02, s13, &k));
                    st(b.add(p + i), vsub(s02, s13, &k));
                    st(b.add(2 * p + i), vadd(d02, d13, &k));
                    st(b.add(3 * p + i), vsub(d02, d13, &k));
                    i += 16;
                }
                if s + 1 != 1 << len {
                    rot = mmul::<P>(rot, t.rate3[(!s).trailing_zeros() as usize]);
                }
            }
        } else {
            let p = 1usize << (h - len - 1);
            for s in 0..1usize << len {
                let b = ptr.add(s * 2 * p);
                let w = mkw::<P>(rot);
                let mut i = 0;
                while i < p {
                    let l = ld(b.add(i));
                    let x = vmulc(ld(b.add(p + i)), &w, &k);
                    st(b.add(i), vadd(l, x, &k));
                    st(b.add(p + i), vsub(l, x, &k));
                    i += 16;
                }
                if s + 1 != 1 << len {
                    rot = mmul::<P>(rot, t.rate2[(!s).trailing_zeros() as usize]);
                }
            }
        }
    }
    let p = perms();
    let rp = r.as_ptr();
    for c in 0..n / 32 {
        let q = ptr.add(32 * c);
        let (mut x0, mut x1) = (ld(q), ld(q.add(16)));
        // span 16
        let tt = vmulc(x1, &mkw::<P>(*rp.add(c)), &k);
        (x0, x1) = (vadd(x0, tt, &k), vsub(x0, tt, &k));
        fwd_layer::<P>(&mut x0, &mut x1, 0, _mm512_permutexvar_epi32(p.tw[0], ld(rp.add(2 * c))), &p, &k);
        fwd_layer::<P>(&mut x0, &mut x1, 1, _mm512_permutexvar_epi32(p.tw[1], ld(rp.add(4 * c))), &p, &k);
        fwd_layer::<P>(&mut x0, &mut x1, 2, _mm512_permutexvar_epi32(p.tw[2], ld(rp.add(8 * c))), &p, &k);
        fwd_layer::<P>(&mut x0, &mut x1, 3, ld(rp.add(16 * c)), &p, &k);
        st(q, x0);
        st(q.add(16), x1);
    }
}

/// Inverse transform without the `1 / n` factor.
///
/// # Panics
/// Panics if AVX-512 is unavailable, the length is not a power of two or it exceeds the 2-adic order.
pub fn intt_unscaled_avx512<const P: u32>(a: &mut [u32]) {
    assert!(available(), "AVX-512 is not available");
    let n = a.len();
    assert!(n.is_power_of_two(), "length {n} is not a power of two");
    let h = n.trailing_zeros() as usize;
    assert!(h <= rank2::<P>() as usize, "length 2^{h} exceeds the 2-adic order of {P} - 1");
    if h < 6 {
        return crate::ntt::intt_unscaled_scalar::<P>(a);
    }
    crate::simd::with_tables_pub::<P, _>(n / 2, |_r, ir| unsafe { intt_impl::<P>(a, ir) })
}

#[target_feature(enable = "avx512f")]
unsafe fn intt_impl<const P: u32>(a: &mut [u32], ir: &[u32]) {
    let n = a.len();
    let h = n.trailing_zeros() as usize;
    let t = info::<P>();
    let k = konst::<P>();
    let one_m = to_mont::<P>(1);
    let iimag = mkw::<P>(to_mont::<P>(t.iroot[2]));
    let ptr = a.as_mut_ptr();
    let p = perms();
    let rp = ir.as_ptr();
    for c in 0..n / 32 {
        let q = ptr.add(32 * c);
        let (mut x0, mut x1) = (ld(q), ld(q.add(16)));
        inv_layer::<P>(&mut x0, &mut x1, 3, ld(rp.add(16 * c)), &p, &k);
        inv_layer::<P>(&mut x0, &mut x1, 2, _mm512_permutexvar_epi32(p.tw[2], ld(rp.add(8 * c))), &p, &k);
        inv_layer::<P>(&mut x0, &mut x1, 1, _mm512_permutexvar_epi32(p.tw[1], ld(rp.add(4 * c))), &p, &k);
        inv_layer::<P>(&mut x0, &mut x1, 0, _mm512_permutexvar_epi32(p.tw[0], ld(rp.add(2 * c))), &p, &k);
        let d = vmulc(vsub(x0, x1, &k), &mkw::<P>(*rp.add(c)), &k);
        st(q, vadd(x0, x1, &k));
        st(q.add(16), d);
    }
    for (len, is4) in plan(h).into_iter().rev() {
        let mut irot = one_m;
        if is4 {
            let p = 1usize << (h - len - 2);
            for s in 0..1usize << len {
                let irot2 = mmul::<P>(irot, irot);
                let irot3 = mmul::<P>(irot2, irot);
                let (i1, i2, i3) = (mkw::<P>(irot), mkw::<P>(irot2), mkw::<P>(irot3));
                let b = ptr.add(s * 4 * p);
                let mut i = 0;
                while i < p {
                    let (a0, a1) = (ld(b.add(i)), ld(b.add(p + i)));
                    let (a2, a3) = (ld(b.add(2 * p + i)), ld(b.add(3 * p + i)));
                    let s01 = vadd(a0, a1, &k);
                    let d01 = vsub(a0, a1, &k);
                    let s23 = vadd(a2, a3, &k);
                    let d23 = vmulc(vsub(a2, a3, &k), &iimag, &k);
                    st(b.add(i), vadd(s01, s23, &k));
                    st(b.add(p + i), vmulc(vadd(d01, d23, &k), &i1, &k));
                    st(b.add(2 * p + i), vmulc(vsub(s01, s23, &k), &i2, &k));
                    st(b.add(3 * p + i), vmulc(vsub(d01, d23, &k), &i3, &k));
                    i += 16;
                }
                if s + 1 != 1 << len {
                    irot = mmul::<P>(irot, t.irate3[(!s).trailing_zeros() as usize]);
                }
            }
        } else {
            let p = 1usize << (h - len - 1);
            for s in 0..1usize << len {
                let b = ptr.add(s * 2 * p);
                let w = mkw::<P>(irot);
                let mut i = 0;
                while i < p {
                    let l = ld(b.add(i));
                    let x = ld(b.add(p + i));
                    st(b.add(i), vadd(l, x, &k));
                    st(b.add(p + i), vmulc(vsub(l, x, &k), &w, &k));
                    i += 16;
                }
                if s + 1 != 1 << len {
                    irot = mmul::<P>(irot, t.irate2[(!s).trailing_zeros() as usize]);
                }
            }
        }
    }
}
