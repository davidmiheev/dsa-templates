//! Ranked (polynomial in the set size) representation of set power series: entry `[S][r]` of the
//! zeta transform holds `sum_{T subset S, |T| = r} a[T]`, stored with stride `n + 1` so that the
//! transforms are vector additions and the per-mask work touches contiguous memory. Products,
//! `exp` and `log` are then polynomial operations in `r`, truncated at degree `n`.

use super::{add, pow_mod, sub};

#[inline(always)]
fn reduce<const M: u32>(acc: u64, lim: u64) -> u64 {
    let t = acc.wrapping_sub(lim);
    acc.min(t)
}

fn zeta_ranked<const M: u32>(a: &mut [u32], n: usize) {
    let w = n + 1;
    for bit in 0..n {
        let step = w << bit;
        for chunk in a.chunks_exact_mut(2 * step) {
            let (lo, hi) = chunk.split_at_mut(step);
            for (h, &l) in hi.iter_mut().zip(lo.iter()) {
                *h = add::<M>(*h, l);
            }
        }
    }
}

fn mobius_ranked<const M: u32>(a: &mut [u32], n: usize) {
    let w = n + 1;
    for bit in 0..n {
        let step = w << bit;
        for chunk in a.chunks_exact_mut(2 * step) {
            let (lo, hi) = chunk.split_at_mut(step);
            for (h, &l) in hi.iter_mut().zip(lo.iter()) {
                *h = sub::<M>(*h, l);
            }
        }
    }
}

/// Ranked zeta transform of `a` (length `2^n`).
fn to_ranked<const M: u32>(a: &[u32], n: usize) -> Vec<u32> {
    let w = n + 1;
    let mut r = vec![0u32; w << n];
    for (s, &x) in a.iter().enumerate() {
        r[s * w + s.count_ones() as usize] = x;
    }
    zeta_ranked::<M>(&mut r, n);
    r
}

/// Inverse of [`to_ranked`] on the image of a product / `exp` / `log`: reads the rank `|S|`.
fn from_ranked<const M: u32>(mut r: Vec<u32>, n: usize) -> Vec<u32> {
    let w = n + 1;
    mobius_ranked::<M>(&mut r, n);
    (0..1usize << n).map(|s| r[s * w + s.count_ones() as usize]).collect()
}

/// `out[k] = sum_{i + j = k} x[i] y[j]` for `k <= n`, where `x`, `y` are nonzero only up to
/// index `p`.
#[inline(always)]
fn mul_trunc<const M: u32>(x: &[u32], y: &[u32], p: usize, n: usize, out: &mut [u32]) {
    let lim = (M as u64) * (M as u64);
    for (k, o) in out.iter_mut().enumerate().take(n + 1) {
        if k > 2 * p {
            *o = 0;
            continue;
        }
        let mut acc = 0u64;
        for i in k.saturating_sub(p)..=k.min(p) {
            acc = reduce::<M>(acc + x[i] as u64 * y[k - i] as u64, lim);
        }
        *o = (acc % M as u64) as u32;
    }
}

/// Multiplication by a fixed set power series.
struct Multiplier<const M: u32> {
    n: usize,
    hat: Vec<u32>,
}

impl<const M: u32> Multiplier<M> {
    fn new(b: &[u32]) -> Self {
        assert!(b.len().is_power_of_two(), "length is not a power of two");
        let n = b.len().trailing_zeros() as usize;
        Self { n, hat: to_ranked::<M>(b, n) }
    }

    /// `b * a` in the subset convolution ring.
    fn mul(&self, a: &[u32]) -> Vec<u32> {
        let (n, w) = (self.n, self.n + 1);
        assert_eq!(a.len(), 1 << n, "length mismatch");
        let mut r = to_ranked::<M>(a, n);
        let mut tmp = vec![0u32; w];
        for (s, blk) in r.chunks_exact_mut(w).enumerate() {
            let p = s.count_ones() as usize;
            mul_trunc::<M>(&self.hat[s * w..(s + 1) * w], blk, p, n, &mut tmp);
            blk.copy_from_slice(&tmp);
        }
        from_ranked::<M>(r, n)
    }

    /// The transpose of `mul`: `v[T] = sum_{U disjoint from T} b[U] a[T | U]`.
    fn mul_transposed(&self, a: &[u32]) -> Vec<u32> {
        let mut rev: Vec<u32> = a.to_vec();
        rev.reverse();
        let mut v = self.mul(&rev);
        v.reverse();
        v
    }
}

/// `c[S] = sum_{T subset S} a[T] b[S \ T]`. Both slices need the same power-of-two length.
pub fn subset_convolution<const M: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len(), "length mismatch");
    Multiplier::<M>::new(b).mul(a)
}

fn inverses<const M: u32>(n: usize) -> Vec<u32> {
    assert!((M as usize) > n, "modulus too small");
    (0..=n)
        .map(|k| if k == 0 { 0 } else { pow_mod(k as u64, M as u64 - 2, M as u64) as u32 })
        .collect()
}

/// `exp(s) = sum_m s^m / m!` in the subset convolution ring. Requires `s[0] = 0`.
pub fn subset_exp<const M: u32>(s: &[u32]) -> Vec<u32> {
    assert!(s.len().is_power_of_two(), "length is not a power of two");
    assert_eq!(s[0], 0, "constant term must be zero");
    let n = s.len().trailing_zeros() as usize;
    let w = n + 1;
    let inv = inverses::<M>(n);
    let lim = (M as u64) * (M as u64);
    let mut r = to_ranked::<M>(s, n);
    let mut jf = vec![0u64; w];
    let mut g = vec![0u32; w];
    for (mask, blk) in r.chunks_exact_mut(w).enumerate() {
        let p = mask.count_ones() as usize;
        for j in 1..=p {
            jf[j] = blk[j] as u64 * j as u64 % M as u64;
        }
        g[0] = 1;
        for k in 1..=n {
            let mut acc = 0u64;
            for j in 1..=k.min(p) {
                acc = reduce::<M>(acc + jf[j] * g[k - j] as u64, lim);
            }
            g[k] = ((acc % M as u64) * inv[k] as u64 % M as u64) as u32;
        }
        blk.copy_from_slice(&g);
    }
    from_ranked::<M>(r, n)
}

/// `log(s)` in the subset convolution ring, the inverse of [`subset_exp`]. Requires `s[0] = 1`.
pub fn subset_log<const M: u32>(s: &[u32]) -> Vec<u32> {
    assert!(s.len().is_power_of_two(), "length is not a power of two");
    assert_eq!(s[0], 1, "constant term must be one");
    let n = s.len().trailing_zeros() as usize;
    let w = n + 1;
    let inv = inverses::<M>(n);
    let lim = (M as u64) * (M as u64);
    let mut r = to_ranked::<M>(s, n);
    let mut jg = vec![0u64; w];
    let mut g = vec![0u32; w];
    for (mask, blk) in r.chunks_exact_mut(w).enumerate() {
        let p = mask.count_ones() as usize;
        g[0] = 0;
        for k in 1..=n {
            let mut acc = 0u64;
            for j in k.saturating_sub(p).max(1)..k {
                acc = reduce::<M>(acc + jg[j] * blk[k - j] as u64, lim);
            }
            let t = (acc % M as u64) * inv[k] as u64 % M as u64;
            let fk = if k <= p { blk[k] as u64 } else { 0 };
            g[k] = ((fk + M as u64 - t) % M as u64) as u32;
            jg[k] = g[k] as u64 * k as u64 % M as u64;
        }
        blk.copy_from_slice(&g);
    }
    from_ranked::<M>(r, n)
}

/// `d[k] = f^{(k)}(x0)` for `k = 0..=n` (Taylor coefficients times `k!`).
fn derivatives<const M: u32>(f: &[u32], x0: u32, n: usize) -> Vec<u32> {
    let mut cur: Vec<u64> = f.iter().map(|&x| x as u64).collect();
    let mut out = Vec::with_capacity(n + 1);
    let mut fact = 1u64;
    for k in 0..=n {
        if k > 0 {
            fact = fact * k as u64 % M as u64;
        }
        if cur.is_empty() {
            out.push(0);
            continue;
        }
        // synthetic division by (x - x0): afterwards cur[0] is the remainder f(x0)
        let mut t = 0u64;
        for c in cur.iter_mut().rev() {
            t = (*c + t * x0 as u64) % M as u64;
            *c = t;
        }
        out.push((cur[0] * fact % M as u64) as u32);
        cur.remove(0);
    }
    out
}

/// `f(s)` where `f = sum f[i] x^i` and `s` is a set power series; `O(n^2 2^n + |f| n)`.
pub fn subset_compose<const M: u32>(f: &[u32], s: &[u32]) -> Vec<u32> {
    assert!(s.len().is_power_of_two(), "length is not a power of two");
    let n = s.len().trailing_zeros() as usize;
    // h[k] = f^{(k)}(s restricted to the first j variables), k <= n - j
    let mut h: Vec<Vec<u32>> = derivatives::<M>(f, s[0], n).into_iter().map(|x| vec![x]).collect();
    for j in 1..=n {
        let half = 1usize << (j - 1);
        let mult = Multiplier::<M>::new(&s[half..2 * half]);
        let mut next = Vec::with_capacity(n - j + 1);
        for k in 0..=n - j {
            let mut v = std::mem::take(&mut h[k]);
            let hi = mult.mul(&h[k + 1]);
            v.extend_from_slice(&hi);
            next.push(v);
        }
        h = next;
    }
    h.swap_remove(0)
}

/// `ans[m] = sum_S w[S] * (s^m)[S]` for `m < m_count`; the transpose of [`subset_compose`].
pub fn subset_power_projection<const M: u32>(s: &[u32], w: &[u32], m_count: usize) -> Vec<u32> {
    assert!(s.len().is_power_of_two(), "length is not a power of two");
    assert_eq!(s.len(), w.len(), "length mismatch");
    let n = s.len().trailing_zeros() as usize;
    let mut dual: Vec<Vec<u32>> = vec![w.to_vec()];
    for j in (1..=n).rev() {
        let half = 1usize << (j - 1);
        let mult = Multiplier::<M>::new(&s[half..2 * half]);
        let mut next = vec![vec![0u32; half]; n - j + 2];
        for k in 0..=n - j {
            let (lo, hi) = dual[k].split_at(half);
            for (x, &y) in next[k].iter_mut().zip(lo) {
                *x = add::<M>(*x, y);
            }
            let t = mult.mul_transposed(hi);
            for (x, &y) in next[k + 1].iter_mut().zip(&t) {
                *x = add::<M>(*x, y);
            }
        }
        dual = next;
    }
    // dual[k][0] pairs with f^{(k)}(s0) = sum_i f[i] i!/(i-k)! s0^{i-k}
    let s0 = s[0] as u64;
    let m = M as u64;
    let mut pows = vec![1u64; m_count.max(1)];
    for e in 1..m_count {
        pows[e] = pows[e - 1] * s0 % m;
    }
    (0..m_count)
        .map(|i| {
            let mut ff = 1u64; // i (i-1) ... (i-k+1)
            let mut sum = 0u64;
            for k in 0..=n.min(i) {
                if k > 0 {
                    ff = ff * (i - k + 1) as u64 % m;
                }
                sum = (sum + dual[k][0] as u64 * ff % m * pows[i - k]) % m;
            }
            sum as u32
        })
        .collect()
}
