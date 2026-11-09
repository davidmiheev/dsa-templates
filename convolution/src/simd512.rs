//! Experimental AVX-512 version of the transforms of [`crate::simd`] (sixteen Montgomery lanes).
//!
//! Needs the `avx512` cargo feature (the AVX-512 intrinsics are stable from Rust 1.89) and the
//! `avx512f` CPU flag at run time. Same layer structure and the same results as the other
//! engines: the layers with span >= 32 are vectorised along the span, the last five layers
//! (spans 16, 8, 4, 2, 1) run inside two registers per 32 elements with `vpermt2d` shuffles.

use crate::ntt::{info, mmul, neg_inv, rank2, to_mont};
use std::arch::x86_64::*;

static ENABLED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(true);

/// Switch the dispatch to the AVX-512 transforms on or off (for benchmarks; the AVX2 and portable
/// code give identical results).
pub fn set_enabled(on: bool) {
    ENABLED.store(on, std::sync::atomic::Ordering::Relaxed);
}

/// Whether the AVX-512 transforms can run on this CPU (and were not switched off with
/// [`set_enabled`] or [`crate::simd::set_enabled`]).
pub fn available() -> bool {
    ENABLED.load(std::sync::atomic::Ordering::Relaxed) && crate::simd::enabled() && std::is_x86_feature_detected!("avx512f")
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

/// The shuffle tables of the five in-register layers (spans 16 handled separately), built at
/// compile time: for span `s` in `[8, 4, 2, 1]`, `l` / `h` gather the pair-low / pair-high
/// elements of the 32-element chunk `(x0, x1)`, `o0` / `o1` rebuild `x0` / `x1` from the new
/// pair-low (`aa`) and pair-high (`bb`) vectors, `tw[t]` selects the twiddle of each lane.
struct PermData {
    l: [[i32; 16]; 4],
    h: [[i32; 16]; 4],
    o0: [[i32; 16]; 4],
    o1: [[i32; 16]; 4],
    tw: [[i32; 16]; 3],
}

const fn perm_data() -> PermData {
    let mut d = PermData { l: [[0; 16]; 4], h: [[0; 16]; 4], o0: [[0; 16]; 4], o1: [[0; 16]; 4], tw: [[0; 16]; 3] };
    let spans = [8usize, 4, 2, 1];
    let mut t = 0;
    while t < 4 {
        let s = spans[t];
        let (mut nl, mut nh) = (0, 0);
        let mut j = 0;
        while j < 32 {
            if j & s == 0 {
                d.l[t][nl] = j as i32;
                nl += 1;
            } else {
                d.h[t][nh] = j as i32;
                nh += 1;
            }
            j += 1;
        }
        // element e of the chunk comes from aa[rank among the pair-low] or bb[rank among the pair-high] + 16
        let (mut rl, mut rh) = (0, 0);
        let mut e = 0;
        while e < 32 {
            let src = if e & s == 0 {
                rl += 1;
                rl - 1
            } else {
                rh += 1;
                16 + rh - 1
            };
            if e < 16 {
                d.o0[t][e] = src;
            } else {
                d.o1[t][e - 16] = src;
            }
            e += 1;
        }
        t += 1;
    }
    let mut k = 0;
    while k < 16 {
        d.tw[0][k] = (k / 8) as i32;
        d.tw[1][k] = (k / 4) as i32;
        d.tw[2][k] = (k / 2) as i32;
        k += 1;
    }
    d
}

static PERM: PermData = perm_data();

struct Perm {
    l: [__m512i; 4],
    h: [__m512i; 4],
    tw: [__m512i; 3],
    o0: [__m512i; 4],
    o1: [__m512i; 4],
}

#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn perms() -> Perm {
    let v = |a: &[i32; 16]| _mm512_loadu_si512(a.as_ptr() as *const _);
    Perm {
        l: [v(&PERM.l[0]), v(&PERM.l[1]), v(&PERM.l[2]), v(&PERM.l[3])],
        h: [v(&PERM.h[0]), v(&PERM.h[1]), v(&PERM.h[2]), v(&PERM.h[3])],
        tw: [v(&PERM.tw[0]), v(&PERM.tw[1]), v(&PERM.tw[2])],
        o0: [v(&PERM.o0[0]), v(&PERM.o0[1]), v(&PERM.o0[2]), v(&PERM.o0[3])],
        o1: [v(&PERM.o1[0]), v(&PERM.o1[1]), v(&PERM.o1[2]), v(&PERM.o1[3])],
    }
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
