//! Terms of a linear recurrence at huge indices.

use modint::StaticModInt;
use convolution::{intt, ntt};
use polynomial::fps;

type Z<const M: u32> = StaticModInt<M>;

/// `1 - c_1 x - ... - c_d x^d`.
fn denominator<const M: u32>(c: &[Z<M>]) -> Vec<Z<M>> {
    let mut q = vec![Z::<M>::one()];
    q.extend(c.iter().map(|&x| -x));
    q
}

/// The numerator `A(x) Q(x) mod x^d` of the generating function `A / Q` of `a`.
fn numerator<const M: u32>(a: &[Z<M>], q: &[Z<M>]) -> Vec<Z<M>> {
    let d = a.len();
    fps::mul_trunc(a, q, d)
}

/// `s_k` for the sequence with initial terms `a_0, ..., a_{d-1}` and
/// `s_i = c_1 s_{i-1} + ... + c_d s_{i-d}` (`i >= d`), by Bostan–Mori: `O(d log d log k)`.
///
/// # Panics
/// Panics if `a.len() != c.len()`.
pub fn kth_term<const M: u32>(a: &[Z<M>], c: &[Z<M>], mut k: u64) -> Z<M> {
    let d = a.len();
    assert_eq!(d, c.len(), "need as many initial terms as coefficients");
    if d == 0 {
        return Z::zero();
    }
    let mut q = denominator(c);
    let mut p = numerator(a, &q);
    if d >= SPECTRAL_THRESHOLD {
        return kth_term_spectral(p, q, k);
    }
    while k > 0 {
        let q_neg: Vec<Z<M>> = q.iter().enumerate().map(|(i, &x)| if i % 2 == 1 { -x } else { x }).collect();
        let u = fps::mul(&p, &q_neg);
        let v = fps::mul(&q, &q_neg);
        let parity = (k & 1) as usize;
        p = u.iter().skip(parity).step_by(2).copied().collect();
        q = v.iter().step_by(2).copied().collect();
        p.truncate(d);
        q.truncate(d + 1);
        k >>= 1;
    }
    p.first().copied().unwrap_or_else(Z::zero) / q[0]
}

/// Below this order the plain polynomial products are as fast as the spectral loop.
const SPECTRAL_THRESHOLD: usize = 64;

/// Bostan–Mori on spectra: per step two forward transforms of size `n >= 2d + 2` and two inverse
/// transforms of size `n / 2`. In the transform's bit-reversed order the entries `2m` and
/// `2m + 1` are the values at `x` and `-x`, so `Q(-x)` is a neighbour swap and the even / odd
/// part of a product is read off directly at the squared points `x^2` (the position `m` of the
/// half-size transform).
fn kth_term_spectral<const M: u32>(p: Vec<Z<M>>, q: Vec<Z<M>>, mut k: u64) -> Z<M> {
    let d = p.len();
    let n = (2 * d + 1).next_power_of_two();
    let h = n / 2;
    let m64 = M as u64;
    // 1 / x at the points of the pairs: the spectrum of the monomial x^{n-1}
    let mut inv_x = vec![0u32; n];
    inv_x[n - 1] = 1;
    ntt::<M>(&mut inv_x);
    let half = (M as u64).div_ceil(2);
    let (mut fp, mut fq) = (vec![0u32; n], vec![0u32; n]);
    for (dst, src) in fp.iter_mut().zip(&p) {
        *dst = src.val();
    }
    for (dst, src) in fq.iter_mut().zip(&q) {
        *dst = src.val();
    }
    let (mut pe, mut qe) = (vec![0u32; h], vec![0u32; h]);
    while k > 0 {
        ntt::<M>(&mut fp);
        ntt::<M>(&mut fq);
        let parity = k & 1;
        for m in 0..h {
            let (a, b) = (fq[2 * m] as u64, fq[2 * m + 1] as u64);
            let (u0, u1) = (fp[2 * m] as u64 * b % m64, fp[2 * m + 1] as u64 * a % m64);
            pe[m] = if parity == 0 {
                ((u0 + u1) % m64 * half % m64) as u32
            } else {
                ((u0 + m64 - u1) % m64 * half % m64 * inv_x[2 * m] as u64 % m64) as u32
            };
            qe[m] = (a * b % m64) as u32;
        }
        intt::<M>(&mut pe);
        intt::<M>(&mut qe);
        fp[..d].copy_from_slice(&pe[..d]);
        fp[d..].fill(0);
        fq[..=d].copy_from_slice(&qe[..=d]);
        fq[d + 1..].fill(0);
        k >>= 1;
    }
    Z::<M>::raw(fp[0]) / Z::<M>::raw(fq[0])
}

/// `x^k mod g` for a monic `g` of degree `d >= 1`, by square-and-multiply with a precomputed
/// reciprocal for the divisions.
fn pow_x_mod<const M: u32>(g: &[Z<M>], mut k: u64) -> Vec<Z<M>> {
    let d = g.len() - 1;
    if d == 1 {
        // x = -g_0 modulo (x + g_0): a scalar power
        let (mut r, mut b) = (Z::<M>::one(), -g[0]);
        while k > 0 {
            if k & 1 == 1 {
                r *= b;
            }
            b *= b;
            k >>= 1;
        }
        return vec![r];
    }
    // reciprocal of the reversed g, for quotients of degree up to d - 2
    let mut g_rev = g.to_vec();
    g_rev.reverse();
    let inv = fps::inv(&g_rev, d);
    let reduce = |f: Vec<Z<M>>| -> Vec<Z<M>> {
        if f.len() <= d {
            return f;
        }
        let qn = f.len() - d; // number of quotient coefficients
        let mut f_rev = f.clone();
        f_rev.reverse();
        let mut qr = fps::mul_trunc(&f_rev[..qn], &inv[..qn], qn);
        qr.resize(qn, Z::zero());
        qr.reverse();
        let prod = fps::mul(&qr, g);
        (0..d).map(|i| f[i] - prod[i]).collect()
    };
    // left to right: square, then multiply by x with a shift and one cheap reduction
    let mut result: Vec<Z<M>> = vec![Z::one()];
    for bit in (0..64 - k.leading_zeros()).rev() {
        if result.len() > 1 {
            result = reduce(fps::mul(&result, &result));
        }
        if (k >> bit) & 1 == 1 {
            result.insert(0, Z::zero());
            if result.len() > d {
                let top = result[d];
                for i in 0..d {
                    result[i] -= top * g[i];
                }
                result.truncate(d);
            }
        }
    }
    result.resize(d, Z::zero());
    result
}

/// `s_k, s_{k+1}, ..., s_{k+m-1}` for the recurrence of [`kth_term`], in
/// `O(d log d log k + (d + m) log (d + m))`.
///
/// # Panics
/// Panics if `a.len() != c.len()` or `d = 0`.
pub fn consecutive_terms<const M: u32>(a: &[Z<M>], c: &[Z<M>], k: u64, m: usize) -> Vec<Z<M>> {
    let d = a.len();
    assert!(d >= 1 && d == c.len(), "need as many initial terms as coefficients (at least one)");
    let q = denominator(c);
    let num = numerator(a, &q);
    // the sequence for 2d - 1 terms, as the series num / q
    let series = fps::mul_trunc(&num, &fps::inv(&q, 2 * d), 2 * d);
    // characteristic polynomial g(x) = x^d - c_1 x^{d-1} - ... - c_d, monic
    let mut g: Vec<Z<M>> = vec![Z::zero(); d + 1];
    g[d] = Z::one();
    for (j, &cj) in c.iter().enumerate() {
        g[d - 1 - j] = -cj;
    }
    let r = pow_x_mod(&g, k); // x^k = sum r_j x^j (mod g)
    // s_{k+t} = sum_j r_j s_{j+t}
    let mut r_rev = r.clone();
    r_rev.reverse();
    let corr = fps::mul(&r_rev, &series[..2 * d - 1]);
    let head: Vec<Z<M>> = (0..d).map(|t| corr[d - 1 + t]).collect();
    if m <= d {
        return head[..m].to_vec();
    }
    // continue with the recurrence: head is a new initial segment
    let num2 = numerator(&head, &q);
    fps::mul_trunc(&num2, &fps::inv(&q, m), m)
}
