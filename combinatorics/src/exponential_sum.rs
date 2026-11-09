//! Sums of `r^i i^d`.

use crate::stirling::power_table;
use modint::StaticModInt;
use polynomial::fps;

type Z<const M: u32> = StaticModInt<M>;

/// The value at `x` (an integer, taken modulo `M`) of the polynomial of degree `< values.len()`
/// with the given values at `0, 1, ...`; `inv` holds the inverses of `1..` (`inv[i] = 1 / i`).
/// `O(len)`.
fn lagrange<const M: u32>(values: &[Z<M>], inv: &[Z<M>], x: u64) -> Z<M> {
    let len = values.len();
    let x = x % M as u64;
    if (x as usize) < len {
        return values[x as usize];
    }
    let d = len - 1;
    let xz = Z::<M>::new(x);
    // suffix products of (x - k)
    let mut suffix = vec![Z::<M>::one(); len + 1];
    for k in (0..len).rev() {
        suffix[k] = suffix[k + 1] * (xz - Z::new(k as u64));
    }
    // weights (-1)^{d-j} C(d, j) / d!
    let fact_d = (1..=d as u64).fold(Z::<M>::one(), |a, i| a * Z::new(i));
    let inv_fact_d = fact_d.inv();
    let (mut binom, mut prefix, mut total) = (Z::<M>::one(), Z::<M>::one(), Z::<M>::zero());
    for j in 0..len {
        let mut term = values[j] * binom * prefix * suffix[j + 1];
        if (d - j) % 2 == 1 {
            term = -term;
        }
        total += term;
        prefix *= xz - Z::new(j as u64);
        if j < d {
            binom = binom * Z::new((d - j) as u64) * inv[j + 1];
        }
    }
    total * inv_fact_d
}

/// `S(j) = sum_{i < j} r^i i^d` for `j = 0..=d+1` (with `0^0 = 1`).
fn prefix_sums<const M: u32>(r: Z<M>, d: usize) -> Vec<Z<M>> {
    let pw = power_table::<M>(d + 1, d as u64);
    let mut s = Vec::with_capacity(d + 2);
    s.push(Z::<M>::zero());
    let mut rp = Z::<M>::one();
    for &p in pw.iter() {
        let last = *s.last().expect("non-empty");
        s.push(last + rp * p);
        rp *= r;
    }
    s
}

/// The constant `c = F(0)` where `F` is the polynomial of degree `d` with
/// `sum_{i<n} r^i i^d = r^n F(n) - c`, for `r != 0, 1`.
fn solve_constant<const M: u32>(s: &[Z<M>], inv: &[Z<M>], r: Z<M>, d: usize) -> Z<M> {
    let ir = r.inv();
    // F(j) = (S(j) + c) ir^j has a vanishing (d+1)-th difference:
    // sum_j (-1)^{d+1-j} C(d+1, j) F(j) = 0
    let (mut binom, mut irp, mut num) = (Z::<M>::one(), Z::<M>::one(), Z::<M>::zero());
    for j in 0..=d + 1 {
        let mut term = binom * s[j] * irp;
        if (d + 1 - j) % 2 == 1 {
            term = -term;
        }
        num += term;
        irp *= ir;
        if j <= d {
            binom = binom * Z::new((d + 1 - j) as u64) * inv[j + 1];
        }
    }
    let den = (ir - Z::one()).pow(d as u64 + 1);
    -num / den
}

fn check<const M: u32>(d: usize) {
    assert!((M as usize) > d + 3, "the modulus must exceed d + 3");
}

/// `sum_{i=0}^{n-1} r^i i^d mod M` (`0^0 = 1`) for a prime `M > d + 3`, in `O(d)` time and memory
/// plus `O(d)` modular operations for the powers `i^d` (a linear sieve).
pub fn sum_exponential_times_polynomial<const M: u32>(r: u32, d: usize, n: u64) -> u32 {
    check::<M>(d);
    if n == 0 {
        return 0;
    }
    let r = Z::<M>::new(r as u64);
    if r == Z::zero() {
        return (d == 0) as u32; // only the term i = 0 survives
    }
    let inv = fps::inverses::<M>(d + 3);
    let s = prefix_sums::<M>(r, d);
    if r == Z::one() {
        // a polynomial of degree d + 1 in n
        return lagrange(&s, &inv, n).val();
    }
    let c = solve_constant::<M>(&s, &inv, r, d);
    let ir = r.inv();
    let mut f = Vec::with_capacity(d + 1);
    let mut irp = Z::<M>::one();
    for &sj in s.iter().take(d + 1) {
        f.push((sj + c) * irp);
        irp *= ir;
    }
    (r.pow(n) * lagrange(&f, &inv, n) - c).val()
}

/// `sum_{i=0}^{inf} r~^i i^d mod M` where `r~` is the rational with `|r~| < 1` congruent to `r`
/// (`r != 1`); the value is `-F(0)`, the limit of the finite sum as `r^n -> 0`.
pub fn sum_exponential_times_polynomial_limit<const M: u32>(r: u32, d: usize) -> u32 {
    check::<M>(d);
    let r = Z::<M>::new(r as u64);
    assert!(r != Z::one(), "r must not be 1");
    if r == Z::zero() {
        return (d == 0) as u32;
    }
    let inv = fps::inverses::<M>(d + 3);
    let s = prefix_sums::<M>(r, d);
    (-solve_constant::<M>(&s, &inv, r, d)).val()
}
