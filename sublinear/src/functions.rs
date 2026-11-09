//! Classical sums computed with [`sum_by_division`](super::sum_by_division): the totient and
//! Mertens functions and the number of squarefree integers.

use super::division::sum_by_division;
use super::quotient::isqrt;
use modint::StaticModInt;
use prime_sieve::{mobius_table, totient_table};

/// `sum_{i <= n} phi(i) mod M` in `O(n^{2/3})`.
pub fn totient_sum<const M: u32>(n: u64) -> u64 {
    let bound = n.min(((n as f64).powf(2.0 / 3.0) * 2.0) as u64 + 16);
    let phi = totient_table(bound as usize);
    let mut prefix = vec![StaticModInt::<M>::new(0); bound as usize + 1];
    for i in 1..=bound as usize {
        prefix[i] = prefix[i - 1] + StaticModInt::new(phi[i] as u64);
    }
    // phi * 1 = id
    let g = |v: u64| StaticModInt::<M>::new(v);
    let h = |v: u64| StaticModInt::<M>::new((v as u128 * (v as u128 + 1) / 2 % M as u128) as u64);
    sum_by_division(n, bound, |v| prefix[v as usize], g, h).val() as u64
}

/// The Mertens function `sum_{i <= x} mu(i)` for `x` beyond the table, where `table[i]` holds
/// the prefix sums up to `i = table.len() - 1`.
pub fn mertens(x: u64, table: &[i32]) -> i64 {
    let bound = table.len() as u64 - 1;
    // mu * 1 = [n == 1]
    sum_by_division::<i64>(x, bound, |v| table[v as usize] as i64, |v| v as i64, |v| (v >= 1) as i64)
}

/// The number of squarefree integers in `1..=n`: `sum_{d <= sqrt n} mu(d) floor(n / d^2)`, where
/// the terms with large `d` are regrouped as `sum_k (M(sqrt(n / k)) - M(L))` over Mertens values.
/// About `0.3` s for `n = 10^{16}` and `2` s for `10^{18}`.
pub fn squarefree_count(n: u64) -> u64 {
    assert!(n >= 1, "n must be positive");
    let s = isqrt(n);
    let l = s.min(((n as f64) * 4e5).cbrt() as u64).max(1);
    let mu = mobius_table(l as usize);
    // Mertens values up to `b` for the large queries
    let b = l.min(40_000_000) as usize;
    let mut table = vec![0i32; b + 1];
    for i in 1..=b {
        table[i] = table[i - 1] + mu[i] as i32;
    }
    let mut total = 0i64;
    let mut mertens_l = 0i64;
    for d in 1..=l {
        let m = mu[d as usize];
        mertens_l += m as i64;
        if m != 0 {
            total += m as i64 * (n / (d * d)) as i64;
        }
    }
    // d > l: sum over k >= 1 with sqrt(n / k) > l of (M(isqrt(n / k)) - M(l))
    let mut k = 1u64;
    loop {
        let x = isqrt(n / k);
        if x <= l {
            break;
        }
        let mx = if x as usize <= b { table[x as usize] as i64 } else { mertens(x, &table) };
        total += mx - mertens_l;
        k += 1;
    }
    total as u64
}
