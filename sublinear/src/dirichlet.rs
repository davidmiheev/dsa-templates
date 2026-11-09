//! Dirichlet products and inverses of functions given by their prefix sums at the floor quotients
//! of `n`, in `O(n^{3/4})`.

use super::quotient::{isqrt, Quotients};

fn pow_mod(mut b: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1 % m;
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

/// The values `f(i) = F(i) - F(i - 1)` from the prefix sums `F(1), F(2), ...` (`F(0)` is taken as 0
/// whatever the table holds in slot 0).
fn differences<const M: u64>(prefix: &[u64]) -> Vec<u32> {
    (0..prefix.len())
        .map(|i| match i {
            0 => 0,
            1 => prefix[1] as u32,
            _ => ((prefix[i] + M - prefix[i - 1]) % M) as u32,
        })
        .collect()
}

/// `[F, G]` pairs stored together: one cache line serves both lookups of a term.
struct Pairs {
    r: u64,
    small: Vec<[u32; 2]>,
    large: Vec<[u32; 2]>,
    /// `inv[i] = 1 / i`
    inv: Vec<f64>,
}

impl Pairs {
    fn new(n: u64) -> Self {
        let r = isqrt(n);
        let inv = (0..=r).map(|i| if i == 0 { 0.0 } else { 1.0 / i as f64 }).collect();
        Pairs { r, small: vec![[0, 0]; r as usize + 1], large: vec![[0, 0]; r as usize + 1], inv }
    }
}

/// The sum over `i <= s = isqrt(x)` of `df[i] * G(x / i) + dg[i] * F(x / i)`, minus
/// `F(s) G(s)`, i.e. `H(x)` of the product; `k` is the index of the large entry (`x = n / k`) or 0
/// for `x <= r`. With `skip_first`, the term `i = 1` is left out of the first sum (used by the
/// inverse, which solves for it).
#[inline(always)]
fn hyperbola<const M: u64>(t: &Pairs, df: &[u32], dg: &[u32], x: u64, k: usize, skip_first: bool) -> u64 {
    let s = isqrt(x) as usize;
    let lim = 2 * M * M;
    let xf = x as f64 + 0.5;
    let mut acc = 0u64;
    let direct = (t.r as usize).checked_div(k).map_or(0, |d| d.min(s));
    let first = if skip_first { 2 } else { 1 };
    // i <= direct: x / i = n / (k i) is the large entry k * i
    for i in 1..=direct {
        let [f, g] = unsafe { *t.large.get_unchecked(k * i) };
        let (a, b) = unsafe { (*df.get_unchecked(i), *dg.get_unchecked(i)) };
        let term = if i < first { b as u64 * f as u64 } else { a as u64 * g as u64 + b as u64 * f as u64 };
        acc += term;
        acc = acc.min(acc.wrapping_sub(lim));
    }
    for i in direct + 1..=s {
        let q = (xf * unsafe { *t.inv.get_unchecked(i) }) as usize;
        let [f, g] = unsafe { *t.small.get_unchecked(q) };
        let (a, b) = unsafe { (*df.get_unchecked(i), *dg.get_unchecked(i)) };
        let term = if i < first { b as u64 * f as u64 } else { a as u64 * g as u64 + b as u64 * f as u64 };
        acc += term;
        acc = acc.min(acc.wrapping_sub(lim));
    }
    acc % M
}

/// The prefix sums of the Dirichlet product `f * g` at all of `Q_n`, given those of `f` and `g`
/// (all modulo the prime `M < 2^31`).
pub fn dirichlet_product<const M: u64>(f: &Quotients<u64>, g: &Quotients<u64>) -> Quotients<u64> {
    assert_eq!(f.n(), g.n(), "tables for different n");
    assert!(M < 1 << 31, "the modulus must be below 2^31");
    let n = f.n();
    let mut t = Pairs::new(n);
    let r = t.r as usize;
    for i in 0..=r {
        t.small[i] = [f.small[i] as u32, g.small[i] as u32];
        t.large[i] = [f.large[i] as u32, g.large[i] as u32];
    }
    let df = differences::<M>(&f.small);
    let dg = differences::<M>(&g.small);
    let mut out = Quotients::new(n, 0u64);
    for v in 1..=r {
        let s = isqrt(v as u64) as usize;
        let sub = f.small[s] * g.small[s] % M;
        out.small[v] = (hyperbola::<M>(&t, &df, &dg, v as u64, 0, false) + M - sub) % M;
    }
    for k in 1..=r {
        let x = n / k as u64;
        let s = isqrt(x) as usize;
        let sub = f.small[s] * g.small[s] % M;
        out.large[k] = (hyperbola::<M>(&t, &df, &dg, x, k, false) + M - sub) % M;
    }
    out
}

/// The prefix sums of the Dirichlet inverse `g` of `f` (`f(1) != 0`) at all of `Q_n`, given
/// those of `f`, modulo the prime `M < 2^31`.
pub fn dirichlet_inverse<const M: u64>(f: &Quotients<u64>) -> Quotients<u64> {
    assert!(M < 1 << 31, "the modulus must be below 2^31");
    let n = f.n();
    assert!(!f.small[1].is_multiple_of(M), "f(1) must be invertible");
    let inv_f1 = pow_mod(f.small[1], M - 2, M);
    let mut t = Pairs::new(n);
    let r = t.r as usize;
    for i in 0..=r {
        t.small[i][0] = f.small[i] as u32;
        t.large[i][0] = f.large[i] as u32;
    }
    let df = differences::<M>(&f.small);
    let mut dg = vec![0u32; r + 1];
    let mut out = Quotients::new(n, 0u64);
    // value x: sum_{i, j: ij <= x} f(i) g(j) = 1, with G(x) the only unknown:
    // G(x) f(1) = 1 - [ sum_{2 <= i <= s} f(i) G(x / i) + sum_{j <= s} g(j) F(x / j) - F(s) G(s) ]
    let solve = |t: &Pairs, dg: &[u32], x: u64, k: usize, gs: &dyn Fn(usize) -> u64| -> u64 {
        let s = isqrt(x) as usize;
        let rest = hyperbola::<M>(t, &df, dg, x, k, true);
        let sub = f.small[s] * gs(s) % M;
        let rest = (rest + M - sub) % M;
        (1 + M - rest) % M * inv_f1 % M
    };
    for v in 1..=r {
        let gv = {
            let small_g = |i: usize| out.small[i];
            solve(&t, &dg, v as u64, 0, &small_g)
        };
        out.small[v] = gv;
        t.small[v][1] = gv as u32;
        dg[v] = if v == 1 { gv as u32 } else { ((gv + M - out.small[v - 1]) % M) as u32 };
    }
    for k in (1..=r).rev() {
        let x = n / k as u64;
        let gx = {
            let small_g = |i: usize| out.small[i];
            solve(&t, &dg, x, k, &small_g)
        };
        out.large[k] = gx;
        t.large[k][1] = gx as u32;
    }
    out
}
