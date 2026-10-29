//! AVX2 versions of the radix-4 transforms of [`crate::ntt`] (eight Montgomery lanes at a time).
//!
//! The arithmetic and the layer structure are the same as in the scalar code (so the outputs are
//! identical word for word): the layers whose butterfly span is at least 16 are vectorised along
//! the span with broadcast twiddles, and the last four layers (spans 8, 4, 2, 1) run inside two
//! registers per 16 elements with lane shuffles. The twiddles of the last four layers come from a
//! per-thread table `R[s]` of the block rotations (shared by all sizes).

use crate::ntt::{info, mmul, neg_inv, rank2, to_mont};
use std::arch::x86_64::*;
use std::cell::RefCell;

static ENABLED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(true);

/// Switch the dispatch to the AVX2 transforms on or off (for benchmarks: the portable code
/// gives identical results).
pub fn set_enabled(on: bool) {
    ENABLED.store(on, std::sync::atomic::Ordering::Relaxed);
}

pub(crate) fn enabled() -> bool {
    ENABLED.load(std::sync::atomic::Ordering::Relaxed)
}

/// Whether the AVX2 transforms can run on this CPU (and were not switched off).
pub fn available() -> bool {
    enabled() && std::is_x86_feature_detected!("avx2")
}

/// Rotations `R[s]` of the radix-2 layers (and their inverses), Montgomery form.
struct Tables {
    r: Vec<u32>,
    ir: Vec<u32>,
    info: crate::ntt::Info,
    c: Consts,
}

thread_local! {
    static CACHE: RefCell<Vec<(u32, Tables)>> = const { RefCell::new(Vec::new()) };
}

fn build_tables<const P: u32>(len: usize) -> Tables {
    let t = info::<P>();
    let one = to_mont::<P>(1);
    let (mut r, mut ir) = (Vec::with_capacity(len), Vec::with_capacity(len));
    let (mut x, mut y) = (one, one);
    for s in 0..len {
        r.push(x);
        ir.push(y);
        let z = (!s).trailing_zeros() as usize;
        x = mmul::<P>(x, t.rate2[z.min(31)]);
        y = mmul::<P>(y, t.irate2[z.min(31)]);
    }
    Tables { r, ir, info: info::<P>(), c: consts::<P>() }
}

/// The rotation tables `(R, IR)` for other engines.
pub(crate) fn with_tables_pub<const P: u32, T>(len: usize, f: impl FnOnce(&[u32], &[u32]) -> T) -> T {
    with_tables::<P, _>(len, |t| f(&t.r, &t.ir))
}

/// Runs `f` with tables holding at least `len` rotations for the prime `P`.
fn with_tables<const P: u32, T>(len: usize, f: impl FnOnce(&Tables) -> T) -> T {
    CACHE.with(|c| {
        let mut c = c.borrow_mut();
        let idx = match c.iter().position(|(p, _)| *p == P) {
            Some(i) => i,
            None => {
                c.push((P, build_tables::<P>(0)));
                c.len() - 1
            }
        };
        if c[idx].1.r.len() < len {
            c[idx].1 = build_tables::<P>(len.next_power_of_two());
        }
        f(&c[idx].1)
    })
}

struct K {
    p: __m256i,
    ninv: __m256i,
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn konst<const P: u32>() -> K {
    K { p: _mm256_set1_epi32(P as i32), ninv: _mm256_set1_epi32(neg_inv(P) as i32) }
}

/// Montgomery product `a * b / 2^32` of lanes in `[0, P)`.
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn vmul(a: __m256i, b: __m256i, k: &K) -> __m256i {
    let a1 = _mm256_srli_epi64(a, 32);
    let b1 = _mm256_srli_epi64(b, 32);
    let t0 = _mm256_mul_epu32(a, b);
    let t1 = _mm256_mul_epu32(a1, b1);
    let m0 = _mm256_mul_epu32(t0, k.ninv);
    let m1 = _mm256_mul_epu32(t1, k.ninv);
    let u0 = _mm256_add_epi64(t0, _mm256_mul_epu32(m0, k.p));
    let u1 = _mm256_add_epi64(t1, _mm256_mul_epu32(m1, k.p));
    let r = _mm256_blend_epi32(_mm256_srli_epi64(u0, 32), u1, 0b1010_1010);
    _mm256_min_epu32(r, _mm256_sub_epi32(r, k.p))
}

/// A broadcast twiddle with its precomputed `w * (-P^-1) mod 2^32`.
struct W {
    w: __m256i,
    wn: __m256i,
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn mkw<const P: u32>(w: u32) -> W {
    W { w: bc(w), wn: bc(w.wrapping_mul(neg_inv(P))) }
}

/// `a * w / 2^32` for a broadcast twiddle `w` (the quotient `m` comes straight from `a`).
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn vmulc(a: __m256i, w: &W, k: &K) -> __m256i {
    let a1 = _mm256_srli_epi64(a, 32);
    let t0 = _mm256_mul_epu32(a, w.w);
    let t1 = _mm256_mul_epu32(a1, w.w);
    let m0 = _mm256_mul_epu32(a, w.wn);
    let m1 = _mm256_mul_epu32(a1, w.wn);
    let u0 = _mm256_add_epi64(t0, _mm256_mul_epu32(m0, k.p));
    let u1 = _mm256_add_epi64(t1, _mm256_mul_epu32(m1, k.p));
    let r = _mm256_blend_epi32(_mm256_srli_epi64(u0, 32), u1, 0b1010_1010);
    _mm256_min_epu32(r, _mm256_sub_epi32(r, k.p))
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn vadd(a: __m256i, b: __m256i, k: &K) -> __m256i {
    let s = _mm256_add_epi32(a, b);
    _mm256_min_epu32(s, _mm256_sub_epi32(s, k.p))
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn vsub(a: __m256i, b: __m256i, k: &K) -> __m256i {
    let d = _mm256_sub_epi32(a, b);
    _mm256_min_epu32(d, _mm256_add_epi32(d, k.p))
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn ld(p: *const u32) -> __m256i {
    _mm256_loadu_si256(p as *const __m256i)
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn st(p: *mut u32, v: __m256i) {
    _mm256_storeu_si256(p as *mut __m256i, v)
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn bc(x: u32) -> __m256i {
    _mm256_set1_epi32(x as i32)
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn shuf(a: __m256i, b: __m256i, imm: i32) -> __m256i {
    // shuffle_ps with a compile-time immediate is required: dispatch on the three used values
    let (fa, fb) = (_mm256_castsi256_ps(a), _mm256_castsi256_ps(b));
    let r = match imm {
        0x44 => _mm256_shuffle_ps::<0x44>(fa, fb),
        0xEE => _mm256_shuffle_ps::<0xEE>(fa, fb),
        0x88 => _mm256_shuffle_ps::<0x88>(fa, fb),
        _ => _mm256_shuffle_ps::<0xDD>(fa, fb),
    };
    _mm256_castps_si256(r)
}

/// Elements per cache block: the layers whose block fits are run block by block.
const BLOCK: usize = 1 << 14;

/// The radix of each generic layer group: `(len, is_radix_4)`, for the layers with span >= 16.
fn plan(h: usize) -> Vec<(usize, bool)> {
    let mut stages = Vec::new();
    let mut len = 0;
    while len + 4 < h {
        if h - 4 - len == 1 {
            stages.push((len, false));
            len += 1;
        } else {
            stages.push((len, true));
            len += 2;
        }
    }
    stages
}

/// `x^e` for `x` in Montgomery form.
fn mpow<const P: u32>(x: u32, mut e: u64) -> u32 {
    let (mut r, mut b) = (to_mont::<P>(1), x);
    while e > 0 {
        if e & 1 == 1 {
            r = mmul::<P>(r, b);
        }
        b = mmul::<P>(b, b);
        e >>= 1;
    }
    r
}

/// `X[j]` = the rotation of block `2^j` of the sequence `rot[s+1] = rot[s] * rate[ctz(!s)]`;
/// the sequence is multiplicative in the bits of `s`: `rot[s] = prod_{j in bits(s)} X[j]`.
fn bit_rotations<const P: u32>(rate: &[u32; 32]) -> [u32; 32] {
    let mut x = [0u32; 32];
    for j in 0..26 {
        // rot[2^j] = rate[j] * prod_{i < j} rate[i]^(2^(j-1-i))
        let mut v = rate[j];
        for (i, &r) in rate.iter().enumerate().take(j) {
            v = mmul::<P>(v, mpow::<P>(r, 1u64 << (j - 1 - i)));
        }
        x[j] = v;
    }
    x
}

/// The rotation of block `s0 * 2^shift` of a sequence with the bit rotations `x`.
fn rot_at<const P: u32>(s0: usize, shift: usize, x: &[u32; 32]) -> u32 {
    let mut r = to_mont::<P>(1);
    let mut bits = s0;
    while bits != 0 {
        let j = bits.trailing_zeros() as usize;
        r = mmul::<P>(r, x[j + shift]);
        bits &= bits - 1;
    }
    r
}

struct Consts {
    x2: [u32; 32],
    x4: [u32; 32],
    ix2: [u32; 32],
    ix4: [u32; 32],
}

fn consts<const P: u32>() -> Consts {
    let t = info::<P>();
    Consts {
        x2: bit_rotations::<P>(&t.rate2),
        x4: bit_rotations::<P>(&t.rate3),
        ix2: bit_rotations::<P>(&t.irate2),
        ix4: bit_rotations::<P>(&t.irate3),
    }
}

/// Forward transform, same output as [`crate::ntt::ntt`]. Requires AVX2.
///
/// # Panics
/// Panics if AVX2 is unavailable, the length is not a power of two, or it exceeds the 2-adic order.
pub fn ntt_avx2<const P: u32>(a: &mut [u32]) {
    assert!(available(), "AVX2 is not available");
    let n = a.len();
    assert!(n.is_power_of_two(), "length {n} is not a power of two");
    let h = n.trailing_zeros() as usize;
    assert!(h <= rank2::<P>() as usize, "length 2^{h} exceeds the 2-adic order of {P} - 1");
    if h < 5 {
        return crate::ntt::ntt::<P>(a);
    }
    with_tables::<P, _>(n / 2, |tw| unsafe { ntt_impl::<P>(a, tw) })
}

/// One radix-4 layer group on `nb` consecutive blocks of `4p` elements.
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn fwd4<const P: u32>(base: *mut u32, p: usize, nb: usize, mut rot: u32, t: &crate::ntt::Info, k: &K, imagv: &W) {
    for s in 0..nb {
        let rot2 = mmul::<P>(rot, rot);
        let rot3 = mmul::<P>(rot2, rot);
        let (r1, r2, r3) = (mkw::<P>(rot), mkw::<P>(rot2), mkw::<P>(rot3));
        let b = base.add(s * 4 * p);
        let mut i = 0;
        while i < p {
            let a0 = ld(b.add(i));
            let a1 = vmulc(ld(b.add(p + i)), &r1, k);
            let a2 = vmulc(ld(b.add(2 * p + i)), &r2, k);
            let a3 = vmulc(ld(b.add(3 * p + i)), &r3, k);
            let s02 = vadd(a0, a2, k);
            let d02 = vsub(a0, a2, k);
            let s13 = vadd(a1, a3, k);
            let d13 = vmulc(vsub(a1, a3, k), imagv, k);
            st(b.add(i), vadd(s02, s13, k));
            st(b.add(p + i), vsub(s02, s13, k));
            st(b.add(2 * p + i), vadd(d02, d13, k));
            st(b.add(3 * p + i), vsub(d02, d13, k));
            i += 8;
        }
        if s + 1 != nb {
            rot = mmul::<P>(rot, t.rate3[(!s).trailing_zeros() as usize]);
        }
    }
}

/// One radix-2 layer on `nb` consecutive blocks of `2p` elements.
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn fwd2<const P: u32>(base: *mut u32, p: usize, nb: usize, mut rot: u32, t: &crate::ntt::Info, k: &K) {
    for s in 0..nb {
        let b = base.add(s * 2 * p);
        let r = mkw::<P>(rot);
        let mut i = 0;
        while i < p {
            let l = ld(b.add(i));
            let x = vmulc(ld(b.add(p + i)), &r, k);
            st(b.add(i), vadd(l, x, k));
            st(b.add(p + i), vsub(l, x, k));
            i += 8;
        }
        if s + 1 != nb {
            rot = mmul::<P>(rot, t.rate2[(!s).trailing_zeros() as usize]);
        }
    }
}

/// The last four layers (spans 8, 4, 2, 1) on chunks `c0 .. c0 + nc` of 16 elements.
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn fwd_tail<const P: u32>(ptr: *mut u32, r: *const u32, c0: usize, nc: usize, k: &K) {
    let idx1 = _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7);
    let idx2 = _mm256_setr_epi32(0, 0, 2, 2, 1, 1, 3, 3);
    for c in c0..c0 + nc {
        let q = ptr.add(16 * c);
        let (mut x0, mut x1) = (ld(q), ld(q.add(8)));
        // span 8
        let tt = vmulc(x1, &mkw::<P>(*r.add(c)), k);
        (x0, x1) = (vadd(x0, tt, k), vsub(x0, tt, k));
        // span 4
        let hi = _mm256_permute2x128_si256::<0x31>(x0, x1);
        let lo = _mm256_permute2x128_si256::<0x20>(x0, x1);
        let (r0, r1) = (*r.add(2 * c) as i32, *r.add(2 * c + 1) as i32);
        let tw4 = _mm256_setr_epi32(r0, r0, r0, r0, r1, r1, r1, r1);
        let tt = vmul(hi, tw4, k);
        let (aa, bb) = (vadd(lo, tt, k), vsub(lo, tt, k));
        x0 = _mm256_permute2x128_si256::<0x20>(aa, bb);
        x1 = _mm256_permute2x128_si256::<0x31>(aa, bb);
        // span 2
        let l = shuf(x0, x1, 0x44);
        let hh = shuf(x0, x1, 0xEE);
        let tw2 = _mm256_permutevar8x32_epi32(ld(r.add(4 * c)), idx2);
        let tt = vmul(hh, tw2, k);
        let (aa, bb) = (vadd(l, tt, k), vsub(l, tt, k));
        x0 = _mm256_unpacklo_epi64(aa, bb);
        x1 = _mm256_unpackhi_epi64(aa, bb);
        // span 1
        let l = shuf(x0, x1, 0x88);
        let hh = shuf(x0, x1, 0xDD);
        let tw1 = _mm256_permutevar8x32_epi32(ld(r.add(8 * c)), idx1);
        let tt = vmul(hh, tw1, k);
        let (aa, bb) = (vadd(l, tt, k), vsub(l, tt, k));
        st(q, _mm256_unpacklo_epi32(aa, bb));
        st(q.add(8), _mm256_unpackhi_epi32(aa, bb));
    }
}

#[target_feature(enable = "avx2")]
unsafe fn ntt_impl<const P: u32>(a: &mut [u32], tw: &Tables) {
    let n = a.len();
    let h = n.trailing_zeros() as usize;
    let t = &tw.info;
    let k = konst::<P>();
    let one_m = to_mont::<P>(1);
    let imagv = mkw::<P>(to_mont::<P>(t.root[2]));
    let ptr = a.as_mut_ptr();
    let stages = plan(h);
    // layers whose block exceeds the cache block run over the whole array, the rest block by block
    let split = stages.iter().position(|&(len, _)| (1usize << (h - len)) <= BLOCK).unwrap_or(stages.len());
    for &(len, is4) in &stages[..split] {
        if is4 {
            fwd4::<P>(ptr, 1 << (h - len - 2), 1 << len, one_m, t, &k, &imagv);
        } else {
            fwd2::<P>(ptr, 1 << (h - len - 1), 1 << len, one_m, t, &k);
        }
    }
    if split == stages.len() {
        fwd_tail::<P>(ptr, tw.r.as_ptr(), 0, n / 16, &k);
        return;
    }
    let c = &tw.c;
    let lf = stages[split].0; // sub-blocks are the blocks of that layer
    let sub = 1usize << (h - lf);
    for s0 in 0..1usize << lf {
        let base = ptr.add(s0 * sub);
        for &(len, is4) in &stages[split..] {
            let (nb, shift) = (1usize << (len - lf), len - lf);
            if is4 {
                let rot0 = rot_at::<P>(s0, shift, &c.x4);
                fwd4::<P>(base, 1 << (h - len - 2), nb, rot0, t, &k, &imagv);
            } else {
                let rot0 = rot_at::<P>(s0, shift, &c.x2);
                fwd2::<P>(base, 1 << (h - len - 1), nb, rot0, t, &k);
            }
        }
        fwd_tail::<P>(ptr, tw.r.as_ptr(), s0 * (sub / 16), sub / 16, &k);
    }
}

/// Inverse transform without the `1 / n` factor, same output as [`crate::ntt::intt_unscaled`].
///
/// # Panics
/// Panics if AVX2 is unavailable, the length is not a power of two, or it exceeds the 2-adic order.
pub fn intt_unscaled_avx2<const P: u32>(a: &mut [u32]) {
    assert!(available(), "AVX2 is not available");
    let n = a.len();
    assert!(n.is_power_of_two(), "length {n} is not a power of two");
    let h = n.trailing_zeros() as usize;
    assert!(h <= rank2::<P>() as usize, "length 2^{h} exceeds the 2-adic order of {P} - 1");
    if h < 5 {
        return crate::ntt::intt_unscaled::<P>(a);
    }
    with_tables::<P, _>(n / 2, |tw| unsafe { intt_impl::<P>(a, tw) })
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn inv4<const P: u32>(base: *mut u32, p: usize, nb: usize, mut irot: u32, t: &crate::ntt::Info, k: &K, iimag: &W) {
    for s in 0..nb {
        let irot2 = mmul::<P>(irot, irot);
        let irot3 = mmul::<P>(irot2, irot);
        let (i1, i2, i3) = (mkw::<P>(irot), mkw::<P>(irot2), mkw::<P>(irot3));
        let b = base.add(s * 4 * p);
        let mut i = 0;
        while i < p {
            let (a0, a1) = (ld(b.add(i)), ld(b.add(p + i)));
            let (a2, a3) = (ld(b.add(2 * p + i)), ld(b.add(3 * p + i)));
            let s01 = vadd(a0, a1, k);
            let d01 = vsub(a0, a1, k);
            let s23 = vadd(a2, a3, k);
            let d23 = vmulc(vsub(a2, a3, k), iimag, k);
            st(b.add(i), vadd(s01, s23, k));
            st(b.add(p + i), vmulc(vadd(d01, d23, k), &i1, k));
            st(b.add(2 * p + i), vmulc(vsub(s01, s23, k), &i2, k));
            st(b.add(3 * p + i), vmulc(vsub(d01, d23, k), &i3, k));
            i += 8;
        }
        if s + 1 != nb {
            irot = mmul::<P>(irot, t.irate3[(!s).trailing_zeros() as usize]);
        }
    }
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn inv2<const P: u32>(base: *mut u32, p: usize, nb: usize, mut irot: u32, t: &crate::ntt::Info, k: &K) {
    for s in 0..nb {
        let b = base.add(s * 2 * p);
        let ir = mkw::<P>(irot);
        let mut i = 0;
        while i < p {
            let l = ld(b.add(i));
            let x = ld(b.add(p + i));
            st(b.add(i), vadd(l, x, k));
            st(b.add(p + i), vmulc(vsub(l, x, k), &ir, k));
            i += 8;
        }
        if s + 1 != nb {
            irot = mmul::<P>(irot, t.irate2[(!s).trailing_zeros() as usize]);
        }
    }
}

/// The inverses of the last four layers (spans 1, 2, 4, 8) on chunks `c0 .. c0 + nc`.
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn inv_tail<const P: u32>(ptr: *mut u32, r: *const u32, c0: usize, nc: usize, k: &K) {
    let idx1 = _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7);
    let idx2 = _mm256_setr_epi32(0, 0, 2, 2, 1, 1, 3, 3);
    for c in c0..c0 + nc {
        let q = ptr.add(16 * c);
        let (mut x0, mut x1) = (ld(q), ld(q.add(8)));
        // span 1
        let l = shuf(x0, x1, 0x88);
        let hh = shuf(x0, x1, 0xDD);
        let tw1 = _mm256_permutevar8x32_epi32(ld(r.add(8 * c)), idx1);
        let (aa, bb) = (vadd(l, hh, k), vmul(vsub(l, hh, k), tw1, k));
        x0 = _mm256_unpacklo_epi32(aa, bb);
        x1 = _mm256_unpackhi_epi32(aa, bb);
        // span 2
        let l = shuf(x0, x1, 0x44);
        let hh = shuf(x0, x1, 0xEE);
        let tw2 = _mm256_permutevar8x32_epi32(ld(r.add(4 * c)), idx2);
        let (aa, bb) = (vadd(l, hh, k), vmul(vsub(l, hh, k), tw2, k));
        x0 = _mm256_unpacklo_epi64(aa, bb);
        x1 = _mm256_unpackhi_epi64(aa, bb);
        // span 4
        let hi = _mm256_permute2x128_si256::<0x31>(x0, x1);
        let lo = _mm256_permute2x128_si256::<0x20>(x0, x1);
        let (r0, r1) = (*r.add(2 * c) as i32, *r.add(2 * c + 1) as i32);
        let tw4 = _mm256_setr_epi32(r0, r0, r0, r0, r1, r1, r1, r1);
        let (aa, bb) = (vadd(lo, hi, k), vmul(vsub(lo, hi, k), tw4, k));
        x0 = _mm256_permute2x128_si256::<0x20>(aa, bb);
        x1 = _mm256_permute2x128_si256::<0x31>(aa, bb);
        // span 8
        let d = vmulc(vsub(x0, x1, k), &mkw::<P>(*r.add(c)), k);
        st(q, vadd(x0, x1, k));
        st(q.add(8), d);
    }
}

#[target_feature(enable = "avx2")]
unsafe fn intt_impl<const P: u32>(a: &mut [u32], tw: &Tables) {
    let n = a.len();
    let h = n.trailing_zeros() as usize;
    let t = &tw.info;
    let k = konst::<P>();
    let one_m = to_mont::<P>(1);
    let iimag = mkw::<P>(to_mont::<P>(t.iroot[2]));
    let ptr = a.as_mut_ptr();
    let stages = plan(h);
    let split = stages.iter().position(|&(len, _)| (1usize << (h - len)) <= BLOCK).unwrap_or(stages.len());
    if split == stages.len() {
        inv_tail::<P>(ptr, tw.ir.as_ptr(), 0, n / 16, &k);
    } else {
        let c = &tw.c;
        let lf = stages[split].0;
        let sub = 1usize << (h - lf);
        for s0 in 0..1usize << lf {
            let base = ptr.add(s0 * sub);
            inv_tail::<P>(ptr, tw.ir.as_ptr(), s0 * (sub / 16), sub / 16, &k);
            for &(len, is4) in stages[split..].iter().rev() {
                let (nb, shift) = (1usize << (len - lf), len - lf);
                if is4 {
                    let rot0 = rot_at::<P>(s0, shift, &c.ix4);
                    inv4::<P>(base, 1 << (h - len - 2), nb, rot0, t, &k, &iimag);
                } else {
                    let rot0 = rot_at::<P>(s0, shift, &c.ix2);
                    inv2::<P>(base, 1 << (h - len - 1), nb, rot0, t, &k);
                }
            }
        }
    }
    for &(len, is4) in stages[..split].iter().rev() {
        if is4 {
            inv4::<P>(ptr, 1 << (h - len - 2), 1 << len, one_m, t, &k, &iimag);
        } else {
            inv2::<P>(ptr, 1 << (h - len - 1), 1 << len, one_m, t, &k);
        }
    }
}

/// `a[i] = a[i] * b[i] * scale` for lanes in `[0, P)`, where `scale_m` is the Montgomery form of
/// `scale * 2^32`, i.e. `to_mont(to_mont(scale))`-style input: pass `to_mont::<P>(scale) ` mapped
/// through one more Montgomery form (see [`crate::pointwise_scaled`]).
pub fn pointwise_avx2<const P: u32>(a: &mut [u32], b: &[u32], scale_mm: u32) {
    assert!(available(), "AVX2 is not available");
    assert_eq!(a.len(), b.len());
    unsafe { pointwise_impl::<P>(a, b, scale_mm) }
}

#[target_feature(enable = "avx2")]
unsafe fn pointwise_impl<const P: u32>(a: &mut [u32], b: &[u32], scale_mm: u32) {
    let k = konst::<P>();
    let sc = bc(scale_mm);
    let n = a.len();
    let (pa, pb) = (a.as_mut_ptr(), b.as_ptr());
    let mut i = 0;
    while i + 8 <= n {
        // (a * b / 2^32) * scale_mm / 2^32 with scale_mm = scale * 2^64
        st(pa.add(i), vmul(vmul(ld(pa.add(i)), ld(pb.add(i)), &k), sc, &k));
        i += 8;
    }
    while i < n {
        *pa.add(i) = mmul::<P>(mmul::<P>(*pa.add(i), *pb.add(i)), scale_mm);
        i += 1;
    }
}
