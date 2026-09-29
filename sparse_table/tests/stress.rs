//! Stress test: `SparseTable` for idempotent operations against a naive fold over `[l, r)`.

use sparse_table::{gcd_sparse_table, max_sparse_table, min_sparse_table, SparseTable};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn gcd(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn naive(a: &[i64], l: usize, r: usize, op: fn(i64, i64) -> i64) -> Option<i64> {
    if l >= r || r > a.len() {
        return None;
    }
    // op(x, x) first: idempotent ops return x, but gcd(-x, -x) = |x|, which a
    // one-element range must also report.
    a[l + 1..r].iter().fold(Some(op(a[l], a[l])), |acc, &x| acc.map(|y| op(y, x)))
}

fn and(a: i64, b: i64) -> i64 {
    a & b
}
fn or(a: i64, b: i64) -> i64 {
    a | b
}
fn min(a: i64, b: i64) -> i64 {
    a.min(b)
}
fn max(a: i64, b: i64) -> i64 {
    a.max(b)
}

#[test]
fn matches_naive_fold() {
    let mut st = Stress::new("sparse_table/naive");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(0, 70) as usize;
        let a: Vec<i64> = (0..n).map(|_| rng.range(-100, 100)).collect();
        let nonneg: Vec<i64> = a.iter().map(|x| x.abs()).collect();
        let by_min = min_sparse_table(&a);
        let by_max = max_sparse_table(&a);
        let by_gcd = gcd_sparse_table(&a);
        let by_and = SparseTable::new(&nonneg, and as fn(i64, i64) -> i64);
        let by_or = SparseTable::new(&nonneg, or as fn(i64, i64) -> i64);
        st.case();
        for _ in 0..120 {
            // Includes empty, inverted and out-of-bounds ranges, which must give None.
            let l = rng.below(n + 3);
            let r = rng.below(n + 3);
            assert_eq!(by_min.query(l, r), naive(&a, l, r, min), "min [{l},{r}) case {case} n={n}");
            assert_eq!(by_max.query(l, r), naive(&a, l, r, max), "max [{l},{r}) case {case} n={n}");
            assert_eq!(by_gcd.query(l, r), naive(&a, l, r, gcd), "gcd [{l},{r}) case {case} n={n}");
            assert_eq!(by_and.query(l, r), naive(&nonneg, l, r, and), "and [{l},{r}) case {case} n={n}");
            assert_eq!(by_or.query(l, r), naive(&nonneg, l, r, or), "or [{l},{r}) case {case} n={n}");
            st.ops(5);
        }
    }
    st.done();
}

#[test]
fn large_random_ranges() {
    let mut st = Stress::new("sparse_table/large");
    let n = 200_000;
    let mut rng = Rng::new(5);
    let a: Vec<i64> = (0..n).map(|_| rng.range(-1_000_000_000, 1_000_000_000)).collect();
    let tmin = min_sparse_table(&a);
    let tmax = max_sparse_table(&a);
    st.case();
    for _ in 0..100_000 * scale() {
        let l = rng.below(n);
        let r = (l + 1 + rng.below(400)).min(n);
        assert_eq!(tmin.query(l, r), naive(&a, l, r, min));
        assert_eq!(tmax.query(l, r), naive(&a, l, r, max));
        st.ops(2);
    }
    assert_eq!(tmin.query(0, n), a.iter().copied().min());
    assert_eq!(tmax.query(0, n), a.iter().copied().max());
    st.done();
}
