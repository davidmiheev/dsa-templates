//! Stress test: `SparseTable` for idempotent operations against a naive fold over `[l, r)`.
//! Array sizes cross several block boundaries (the table uses blocks of 32), so single-block,
//! two-block and many-block queries are all exercised.

use sparse_table::{gcd_sparse_table, max_sparse_table, min_sparse_table, SparseTable, BLOCK};
use std::panic::{catch_unwind, AssertUnwindSafe};

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

/// op(x, x) first: idempotent ops return x, but gcd(-x, -x) = |x|, which a
/// one-element range must also report.
fn naive(a: &[i64], l: usize, r: usize, op: fn(i64, i64) -> i64) -> i64 {
    a[l + 1..r].iter().fold(op(a[l], a[l]), |acc, &x| op(acc, x))
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
        // Mix of tiny arrays and ones spanning up to ~12 blocks.
        let n = if case % 4 == 0 { rng.range(1, 70) } else { rng.range(1, 12 * BLOCK as i64) } as usize;
        let a: Vec<i64> = (0..n).map(|_| rng.range(-100, 100)).collect();
        let nonneg: Vec<i64> = a.iter().map(|x| x.abs()).collect();
        let by_min = min_sparse_table(&a);
        let by_max = max_sparse_table(&a);
        let by_gcd = gcd_sparse_table(&a);
        let by_and = SparseTable::new(&nonneg, and as fn(i64, i64) -> i64);
        let by_or = SparseTable::new(&nonneg, or as fn(i64, i64) -> i64);
        st.case();
        assert_eq!((by_min.len(), by_min.is_empty()), (n, false));
        for _ in 0..120 {
            let l = rng.below(n);
            // Mostly arbitrary ranges; some short ones to hit the single-block path.
            let r = if rng.chance(1, 3) { (l + 1 + rng.below(BLOCK)).min(n) } else { l + 1 + rng.below(n - l) };
            let ctx = format!("[{l},{r}) case {case} n={n}");
            assert_eq!(by_min.query(l, r), naive(&a, l, r, min), "min {ctx}");
            assert_eq!(by_max.query(l, r), naive(&a, l, r, max), "max {ctx}");
            assert_eq!(by_gcd.query(l, r), naive(&a, l, r, gcd), "gcd {ctx}");
            assert_eq!(by_and.query(l, r), naive(&nonneg, l, r, and), "and {ctx}");
            assert_eq!(by_or.query(l, r), naive(&nonneg, l, r, or), "or {ctx}");
            st.ops(5);
        }
    }
    st.done();
}

#[test]
fn works_for_other_types() {
    let words = ["pear", "apple", "fig", "kiwi", "plum", "date", "lime"];
    let st = min_sparse_table(&words);
    assert_eq!(st.query(0, 7), "apple");
    assert_eq!(st.query(2, 4), "fig");
    let mx = max_sparse_table(&[1u8, 9, 3]);
    assert_eq!(mx.query(0, 3), 9);
}

#[test]
fn invalid_ranges_panic() {
    for n in [0usize, 1, 5, BLOCK, BLOCK + 1, 3 * BLOCK] {
        let st = min_sparse_table(&vec![1i64; n]);
        assert_eq!(st.is_empty(), n == 0);
        for (l, r) in [(0, 0), (n, n), (1, 1), (3, 2), (0, n + 1), (n, n + 1)] {
            let res = catch_unwind(AssertUnwindSafe(|| st.query(l, r)));
            assert!(res.is_err(), "[{l},{r}) must panic for n={n}");
        }
    }
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
    for _ in 0..1_000 {
        let l = rng.below(n);
        let r = l + 1 + rng.below(n - l);
        assert_eq!(tmin.query(l, r), a[l..r].iter().copied().min().unwrap());
        assert_eq!(tmax.query(l, r), a[l..r].iter().copied().max().unwrap());
        st.ops(2);
    }
    st.done();
}
