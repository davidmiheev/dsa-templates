//! Arithmetic and fast polynomial multiplication in the field `GF(2^64) = F_2[x] / (x^64 + x^4 + x^3 + x + 1)`.
//! Elements are `u64` bit patterns (bit `i` = coefficient of `x^i`); addition is XOR.
//!
//! [`mul`] uses the `PCLMULQDQ` carry-less multiplication instruction when the CPU has it (detected
//! at run time) and a portable shift-and-xor loop otherwise.
//!
//! [`convolution`] multiplies two polynomials over the field with an **additive FFT** over a
//! *Cantor basis* `b_1 = 1`, `b_{i+1}^2 + b_{i+1} = b_i`: the span `V_k` of the first `k` basis
//! elements is mapped by `L(x) = x^2 + x` onto `V_{k-1}` with kernel `{0, 1}`, so with
//! `f(x) = f_0(x^2 + x) + x f_1(x^2 + x)` the values of `f` on `V_k` follow from the values of
//! `f_0`, `f_1` on `V_{k-1}` in one multiplication per pair of points. The decomposition into
//! `f_0, f_1` costs only XORs (it is a halving recursion on the quarters of the coefficient
//! array). Total `O(n log n)` field multiplications and `O(n log^2 n)` XORs.
//!
//! # Example
//!
//! ```
//! use gf2_64::{convolution, mul};
//!
//! assert_eq!(mul(1 << 63, 2), 0b11011); // x^64 = x^4 + x^3 + x + 1
//! let c = convolution(&[1, 1], &[1, 1]); // (1 + x)^2 = 1 + x^2 over a field of characteristic 2
//! assert_eq!(c, vec![1, 0, 1]);
//! ```

use std::sync::OnceLock;

/// Reduce the 128-bit carry-less product `hi * x^64 + lo`.
#[inline(always)]
fn reduce(lo: u64, hi: u64) -> u64 {
    // hi * (x^4 + x^3 + x + 1), with the bits shifted out of 64 folded in once more
    let overflow = (hi >> 63) ^ (hi >> 61) ^ (hi >> 60);
    let folded = overflow ^ (overflow << 1) ^ (overflow << 3) ^ (overflow << 4);
    lo ^ hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4) ^ folded
}

/// Portable carry-less product.
#[inline(always)]
fn clmul_soft(a: u64, b: u64) -> (u64, u64) {
    let (mut lo, mut hi) = (0u64, 0u64);
    for i in 0..64 {
        if (b >> i) & 1 == 1 {
            lo ^= a << i;
            if i > 0 {
                hi ^= a >> (64 - i);
            }
        }
    }
    (lo, hi)
}

/// Product in `GF(2^64)` without the instruction-set shortcut.
pub fn mul_portable(a: u64, b: u64) -> u64 {
    let (lo, hi) = clmul_soft(a, b);
    reduce(lo, hi)
}

#[cfg(target_arch = "x86_64")]
mod hw {
    use std::arch::x86_64::*;

    #[inline]
    #[target_feature(enable = "pclmulqdq")]
    pub unsafe fn clmul(a: u64, b: u64) -> (u64, u64) {
        let r = _mm_clmulepi64_si128(_mm_set_epi64x(0, a as i64), _mm_set_epi64x(0, b as i64), 0);
        (_mm_cvtsi128_si64(r) as u64, _mm_extract_epi64(r, 1) as u64)
    }
}

fn has_clmul() -> bool {
    static FLAG: OnceLock<bool> = OnceLock::new();
    *FLAG.get_or_init(|| {
        #[cfg(target_arch = "x86_64")]
        {
            std::is_x86_feature_detected!("pclmulqdq") && std::is_x86_feature_detected!("sse4.1")
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            false
        }
    })
}

/// Product in `GF(2^64)`.
pub fn mul(a: u64, b: u64) -> u64 {
    #[cfg(target_arch = "x86_64")]
    if has_clmul() {
        // SAFETY: the instruction set was detected at run time.
        let (lo, hi) = unsafe { hw::clmul(a, b) };
        return reduce(lo, hi);
    }
    mul_portable(a, b)
}

/// Square of an element.
pub fn square(a: u64) -> u64 {
    mul(a, a)
}

/// `x` with `x^2 + x = c`, found by Gaussian elimination over `F_2` (the map is linear).
///
/// # Panics
/// Panics if no solution exists (the trace of `c` is 1).
fn solve_quadratic(c: u64) -> u64 {
    // images of the unit vectors
    let mut rows: Vec<(u64, u64)> = Vec::new(); // (image, combination of unit vectors)
    for i in 0..64 {
        let e = 1u64 << i;
        let (mut v, mut combo) = (square(e) ^ e, e);
        for &(r, rc) in &rows {
            if v & (r & r.wrapping_neg()) != 0 {
                v ^= r;
                combo ^= rc;
            }
        }
        if v != 0 {
            rows.push((v, combo));
        }
    }
    let (mut v, mut combo) = (c, 0u64);
    for &(r, rc) in &rows {
        if v & (r & r.wrapping_neg()) != 0 {
            v ^= r;
            combo ^= rc;
        }
    }
    assert!(v == 0, "x^2 + x = {c:#x} has no solution in GF(2^64)");
    combo
}

/// The first `len` elements of the Cantor basis: `b[0] = 1`, `b[i + 1]^2 + b[i + 1] = b[i]`.
pub fn cantor_basis(len: usize) -> Vec<u64> {
    assert!(len <= 64, "GF(2^64) has dimension 64");
    let mut b = vec![1u64];
    while b.len() < len {
        let next = solve_quadratic(*b.last().unwrap());
        b.push(next);
    }
    b.truncate(len);
    b
}

/// Tables for transforms of length up to `2^m`.
struct Tables {
    /// `u[k][p] = point(2 * rev_{k-1}(p))` for `p < 2^{k-1}`: the evaluation point used when
    /// combining at level `k` (position `p` of the child arrays).
    u: Vec<Vec<u64>>,
}

impl Tables {
    fn new(m: usize) -> Self {
        let basis = cantor_basis(m.max(1));
        // points[j] = sum of basis elements selected by the bits of j
        let mut points = vec![0u64; 1 << m];
        for (i, &b) in basis.iter().enumerate().take(m) {
            for j in 0..1usize << i {
                points[(1 << i) + j] = points[j] ^ b;
            }
        }
        let mut u = vec![Vec::new(); m + 1];
        for (k, uk) in u.iter_mut().enumerate().skip(1) {
            let half = 1usize << (k - 1);
            *uk = (0..half)
                .map(|p| {
                    let rev = if k == 1 { 0 } else { p.reverse_bits() >> (usize::BITS as usize - (k - 1)) };
                    points[2 * rev]
                })
                .collect();
        }
        Tables { u }
    }
}

/// Split `f` (length `2^s`, `s >= 1`) into `(f_0 | f_1)` with `f(x) = f_0(x^2 + x) + x f_1(x^2 + x)`.
fn split(f: &mut [u64]) {
    let len = f.len();
    if len <= 2 {
        return;
    }
    let (m, q) = (len / 2, len / 4);
    {
        let (_, rest) = f.split_at_mut(q);
        let (q1, rest) = rest.split_at_mut(q);
        let (q2, q3) = rest.split_at_mut(q);
        for i in 0..q {
            q2[i] ^= q3[i];
            q1[i] ^= q2[i];
        }
    }
    {
        let (a, b) = f.split_at_mut(m);
        split(a);
        split(b);
    }
    // layout [a0 a1 b0 b1] -> [a0 b0 a1 b1]
    let (_, rest) = f.split_at_mut(q);
    let (x, rest) = rest.split_at_mut(q);
    let (y, _) = rest.split_at_mut(q);
    x.swap_with_slice(y);
}

/// Inverse of [`split`].
fn unsplit(f: &mut [u64]) {
    let len = f.len();
    if len <= 2 {
        return;
    }
    let (m, q) = (len / 2, len / 4);
    {
        let (_, rest) = f.split_at_mut(q);
        let (x, rest) = rest.split_at_mut(q);
        let (y, _) = rest.split_at_mut(q);
        x.swap_with_slice(y);
    }
    {
        let (a, b) = f.split_at_mut(m);
        unsplit(a);
        unsplit(b);
    }
    let (_, rest) = f.split_at_mut(q);
    let (q1, rest) = rest.split_at_mut(q);
    let (q2, q3) = rest.split_at_mut(q);
    for i in 0..q {
        q1[i] ^= q2[i];
        q2[i] ^= q3[i];
    }
}

trait Mul {
    fn mul(a: u64, b: u64) -> u64;
}

struct Soft;
impl Mul for Soft {
    #[inline(always)]
    fn mul(a: u64, b: u64) -> u64 {
        mul_portable(a, b)
    }
}

#[cfg(target_arch = "x86_64")]
struct Hard;
#[cfg(target_arch = "x86_64")]
impl Mul for Hard {
    #[inline(always)]
    fn mul(a: u64, b: u64) -> u64 {
        // only instantiated inside functions compiled with the pclmulqdq feature
        let (lo, hi) = unsafe { hw::clmul(a, b) };
        reduce(lo, hi)
    }
}

/// Values of `f` (length `2^k`, monomial coefficients) at the Cantor-basis points, in bit-reversed order.
#[inline(always)]
fn fft<M: Mul>(f: &mut [u64], k: usize, t: &Tables) {
    if k == 0 {
        return;
    }
    split(f);
    let h = f.len() / 2;
    let (f0, f1) = f.split_at_mut(h);
    fft::<M>(f0, k - 1, t);
    fft::<M>(f1, k - 1, t);
    for ((a, b), &u) in f0.iter_mut().zip(f1.iter_mut()).zip(&t.u[k]) {
        let r0 = *a ^ M::mul(*b, u);
        *b ^= r0;
        *a = r0;
    }
}

/// Inverse of [`fft`].
#[inline(always)]
fn ifft<M: Mul>(v: &mut [u64], k: usize, t: &Tables) {
    if k == 0 {
        return;
    }
    let h = v.len() / 2;
    {
        let (v0, v1) = v.split_at_mut(h);
        for ((a, b), &u) in v0.iter_mut().zip(v1.iter_mut()).zip(&t.u[k]) {
            let r0 = *a;
            let f1 = r0 ^ *b;
            *b = f1;
            *a = r0 ^ M::mul(f1, u);
        }
        ifft::<M>(v0, k - 1, t);
        ifft::<M>(v1, k - 1, t);
    }
    unsplit(v);
}

#[inline(always)]
fn convolve_with<M: Mul>(fa: &mut [u64], fb: &mut [u64], k: usize, t: &Tables) {
    fft::<M>(fa, k, t);
    fft::<M>(fb, k, t);
    for (x, &y) in fa.iter_mut().zip(fb.iter()) {
        *x = M::mul(*x, y);
    }
    ifft::<M>(fa, k, t);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "pclmulqdq,sse4.1")]
unsafe fn convolve_hw(fa: &mut [u64], fb: &mut [u64], k: usize, t: &Tables) {
    convolve_with::<Hard>(fa, fb, k, t);
}

/// The product of two polynomials over `GF(2^64)` (coefficient `i` of the result is
/// `sum_{j + l = i} a[j] * b[l]` in the field). Result length `a.len() + b.len() - 1`.
///
/// # Panics
/// Panics if the result is longer than `2^40` (the Cantor basis tables are bounded by memory long before).
pub fn convolution(a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let len = a.len() + b.len() - 1;
    if a.len().min(b.len()) <= 8 {
        let mut c = vec![0u64; len];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                c[i + j] ^= mul(x, y);
            }
        }
        return c;
    }
    let k = len.next_power_of_two().trailing_zeros() as usize;
    let t = Tables::new(k);
    let (mut fa, mut fb) = (vec![0u64; 1 << k], vec![0u64; 1 << k]);
    fa[..a.len()].copy_from_slice(a);
    fb[..b.len()].copy_from_slice(b);
    #[cfg(target_arch = "x86_64")]
    if has_clmul() {
        // SAFETY: the instruction set was detected at run time.
        unsafe { convolve_hw(&mut fa, &mut fb, k, &t) };
        fa.truncate(len);
        return fa;
    }
    convolve_with::<Soft>(&mut fa, &mut fb, k, &t);
    fa.truncate(len);
    fa
}
