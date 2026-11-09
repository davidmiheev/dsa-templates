//! Core operations on formal power series / polynomials over `Z/MZ` (see the crate docs).

use convolution::convolution_mod;
use modint::StaticModInt;
use number_theory::{mod_sqrt, NONE};

type Mint<const M: u32> = StaticModInt<M>;

/// The residues as plain `u32`s, without copying.
fn as_raw<const M: u32>(a: &[Mint<M>]) -> &[u32] {
    // SAFETY: `StaticModInt<M>` is `repr(transparent)` over `u32`.
    unsafe { std::slice::from_raw_parts(a.as_ptr().cast::<u32>(), a.len()) }
}

/// Residues (each `< M`) as `Mint`s, reusing the allocation.
fn from_raw<const M: u32>(a: Vec<u32>) -> Vec<Mint<M>> {
    debug_assert!(a.iter().all(|&x| x < M));
    let mut a = std::mem::ManuallyDrop::new(a);
    // SAFETY: same size and alignment (`repr(transparent)`), the parts come from a live `Vec`.
    unsafe { Vec::from_raw_parts(a.as_mut_ptr().cast::<Mint<M>>(), a.len(), a.capacity()) }
}

/// Product of two polynomials (length `a.len() + b.len() - 1`, empty if either is empty).
pub fn mul<const M: u32>(a: &[Mint<M>], b: &[Mint<M>]) -> Vec<Mint<M>> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    from_raw(convolution_mod::<M>(as_raw(a), as_raw(b)))
}

/// The first `n` coefficients of `a * b`.
pub fn mul_trunc<const M: u32>(a: &[Mint<M>], b: &[Mint<M>], n: usize) -> Vec<Mint<M>> {
    let mut c = mul(&a[..a.len().min(n)], &b[..b.len().min(n)]);
    c.resize(n, Mint::zero());
    c
}

/// `a` truncated or zero-padded to exactly `n` coefficients.
pub fn resized<const M: u32>(a: &[Mint<M>], n: usize) -> Vec<Mint<M>> {
    let mut v = a[..a.len().min(n)].to_vec();
    v.resize(n, Mint::zero());
    v
}

/// Remove trailing zero coefficients.
pub fn trim<const M: u32>(a: &mut Vec<Mint<M>>) {
    while a.last().is_some_and(|x| x.val() == 0) {
        a.pop();
    }
}

/// Formal derivative (one coefficient shorter, empty stays empty).
pub fn derivative<const M: u32>(a: &[Mint<M>]) -> Vec<Mint<M>> {
    (1..a.len()).map(|i| a[i] * Mint::from(i)).collect()
}

/// `1, 1/2, ..., 1/(n-1)` as a table indexed by the denominator (`t[0] = 0`).
pub fn inverses<const M: u32>(n: usize) -> Vec<Mint<M>> {
    let mut t = vec![Mint::<M>::zero(); n.max(2)];
    t[1] = Mint::one();
    for i in 2..n {
        // M = q i + r  =>  1/i = -q / r
        let (q, r) = (M as usize / i, M as usize % i);
        t[i] = -Mint::from(q) * t[r];
    }
    t.truncate(n);
    t
}

/// Antiderivative with zero constant term (one coefficient longer).
pub fn integral<const M: u32>(a: &[Mint<M>]) -> Vec<Mint<M>> {
    let inv = inverses::<M>(a.len() + 1);
    let mut r = vec![Mint::zero(); a.len() + 1];
    for (i, &x) in a.iter().enumerate() {
        r[i + 1] = x * inv[i + 1];
    }
    r
}

/// `1 / f` modulo `x^n` (Newton iteration); `f[0]` must be non-zero.
pub fn inv<const M: u32>(f: &[Mint<M>], n: usize) -> Vec<Mint<M>> {
    assert!(!f.is_empty() && f[0].val() != 0, "the constant term must be invertible");
    let mut g = vec![f[0].inv()];
    let mut m = 1;
    while m < n {
        m *= 2;
        // g <- g (2 - f g) mod x^m
        let fg = mul_trunc(f, &g, m);
        let mut t = vec![Mint::<M>::zero(); m];
        for (i, x) in fg.iter().enumerate() {
            t[i] = -*x;
        }
        t[0] += Mint::from(2u32);
        g = mul_trunc(&g, &t, m);
    }
    g.truncate(n);
    g.resize(n, Mint::zero());
    g
}

/// `log f` modulo `x^n`; `f[0]` must be 1.
pub fn log<const M: u32>(f: &[Mint<M>], n: usize) -> Vec<Mint<M>> {
    assert!(!f.is_empty() && f[0].val() == 1, "log needs constant term 1");
    if n <= 1 {
        return vec![Mint::zero(); n];
    }
    let d = derivative(&f[..f.len().min(n)]);
    let q = mul_trunc(&d, &inv(f, n - 1), n - 1);
    let mut r = integral(&q);
    r.resize(n, Mint::zero());
    r
}

/// `exp f` modulo `x^n`; `f[0]` must be 0.
pub fn exp<const M: u32>(f: &[Mint<M>], n: usize) -> Vec<Mint<M>> {
    assert!(f.is_empty() || f[0].val() == 0, "exp needs constant term 0");
    let mut g = vec![Mint::<M>::one()];
    let mut m = 1;
    while m < n {
        m *= 2;
        // g <- g (1 - log g + f) mod x^m
        let lg = log(&resized(&g, m), m);
        let mut t = vec![Mint::<M>::zero(); m];
        for i in 0..m {
            let fi = if i < f.len() { f[i] } else { Mint::zero() };
            t[i] = fi - lg[i];
        }
        t[0] += Mint::one();
        g = mul_trunc(&g, &t, m);
    }
    g.truncate(n);
    g.resize(n, Mint::zero());
    g
}

/// `f^e` modulo `x^n` for any exponent `e >= 0` (including when `f` has a zero constant term).
pub fn pow<const M: u32>(f: &[Mint<M>], e: u64, n: usize) -> Vec<Mint<M>> {
    let mut result = vec![Mint::<M>::zero(); n];
    if n == 0 {
        return result;
    }
    if e == 0 {
        result[0] = Mint::one();
        return result;
    }
    let Some(k) = f.iter().position(|x| x.val() != 0) else { return result };
    // f = x^k * c * h with h(0) = 1
    if (k as u128) * (e as u128) >= n as u128 {
        return result;
    }
    let shift = k * e as usize;
    let c = f[k];
    let c_inv = c.inv();
    let h: Vec<Mint<M>> = f[k..].iter().map(|&x| x * c_inv).collect();
    let len = n - shift;
    let mut lg = log(&resized(&h, len), len);
    let em = Mint::<M>::from(e % M as u64);
    for x in lg.iter_mut() {
        *x *= em;
    }
    let g = exp(&lg, len);
    let scale = c.pow(e % (M as u64 - 1));
    for i in 0..len {
        result[shift + i] = g[i] * scale;
    }
    result
}

/// A square root of `f` modulo `x^n`, or an **empty vector** if none exists (`n >= 1`).
pub fn sqrt<const M: u32>(f: &[Mint<M>], n: usize) -> Vec<Mint<M>> {
    assert!(n >= 1, "need at least one term");
    let Some(k) = f.iter().position(|x| x.val() != 0) else { return vec![Mint::zero(); n] };
    if k >= n {
        return vec![Mint::zero(); n]; // f = 0 modulo x^n
    }
    if k % 2 == 1 {
        return Vec::new();
    }
    let r = mod_sqrt(f[k].val() as u64, M as u64);
    if r == NONE {
        return Vec::new();
    }
    let h: Vec<Mint<M>> = f[k..].to_vec();
    let len = n - k / 2;
    let mut g = vec![Mint::<M>::new(r)];
    let half = Mint::<M>::new(2).inv();
    let mut m = 1;
    while m < len {
        m *= 2;
        // g <- (g + h / g) / 2 mod x^m
        let ig = inv(&resized(&g, m), m);
        let q = mul_trunc(&h, &ig, m);
        let gg = resized(&g, m);
        g = (0..m).map(|i| (gg[i] + q[i]) * half).collect();
    }
    g.truncate(len);
    let mut result = vec![Mint::zero(); n];
    result[k / 2..k / 2 + len].copy_from_slice(&g);
    result
}

/// Quotient and remainder of `f / g` (`g` non-zero, trailing zeros are ignored). The quotient has
/// `max(deg f - deg g + 1, 0)` coefficients and the remainder fewer than `deg g` (both untrimmed
/// beyond that).
pub fn div_rem<const M: u32>(f: &[Mint<M>], g: &[Mint<M>]) -> (Vec<Mint<M>>, Vec<Mint<M>>) {
    let mut gt = g.to_vec();
    trim(&mut gt);
    assert!(!gt.is_empty(), "division by the zero polynomial");
    let mut ft = f.to_vec();
    trim(&mut ft);
    let (n, m) = (ft.len(), gt.len());
    if n < m {
        return (Vec::new(), ft);
    }
    let qn = n - m + 1;
    // reversed: q^R = f^R / g^R mod x^qn
    let fr: Vec<Mint<M>> = ft.iter().rev().take(qn).copied().collect();
    let gr: Vec<Mint<M>> = gt.iter().rev().copied().collect();
    let mut q = mul_trunc(&fr, &inv(&gr, qn), qn);
    q.reverse();
    let qg = mul(&q, &gt);
    let mut r: Vec<Mint<M>> = (0..m - 1).map(|i| ft[i] - qg[i]).collect();
    trim(&mut r);
    (q, r)
}

/// Product of two bivariate polynomials stored row-major (`rows x cols`, entry `(i, j)` is the
/// coefficient of `x^i y^j`), truncated to `out_rows` rows and `cols` columns, by Kronecker
/// substitution `x = y^(2 cols - 1)`.
pub fn mul_bivariate_trunc<const M: u32>(
    a: &[Mint<M>],
    a_rows: usize,
    b: &[Mint<M>],
    b_rows: usize,
    cols: usize,
    out_rows: usize,
) -> Vec<Mint<M>> {
    let w = 2 * cols - 1;
    let a_rows = a_rows.min(out_rows);
    let b_rows = b_rows.min(out_rows);
    let spread = |v: &[Mint<M>], rows: usize| -> Vec<Mint<M>> {
        let mut out = vec![Mint::<M>::zero(); (rows - 1) * w + cols];
        for i in 0..rows {
            out[i * w..i * w + cols].copy_from_slice(&v[i * cols..(i + 1) * cols]);
        }
        out
    };
    let c = mul(&spread(a, a_rows), &spread(b, b_rows));
    let mut out = vec![Mint::<M>::zero(); out_rows * cols];
    for i in 0..out_rows {
        for j in 0..cols {
            if i * w + j < c.len() {
                out[i * cols + j] = c[i * w + j];
            }
        }
    }
    out
}

/// `1 / f` for a bivariate series `f` (`rows x cols`, row-major, `f(0,0) != 0`) truncated to the
/// same shape: Newton iteration in `x` with Kronecker-substituted products.
pub fn inv_bivariate<const M: u32>(f: &[Mint<M>], rows: usize, cols: usize) -> Vec<Mint<M>> {
    assert_eq!(f.len(), rows * cols, "shape mismatch");
    assert!(f[0].val() != 0, "the constant term must be invertible");
    let mut g = inv(&f[..cols], cols); // 1 row
    let mut have = 1;
    while have < rows {
        let r = (2 * have).min(rows);
        // g <- g (2 - f g) mod (x^r, y^cols)
        let fg = mul_bivariate_trunc(&f[..r.min(rows) * cols], r, &g, have, cols, r);
        let mut t: Vec<Mint<M>> = fg.iter().map(|&x| -x).collect();
        t[0] += Mint::from(2u32);
        g = mul_bivariate_trunc(&g, have, &t, r, cols, r);
        have = r;
    }
    g.truncate(rows * cols);
    g
}
