//! Stress test: `LiChaoTree` against a linear scan over all inserted lines and segments.

use li_chao_tree::LiChaoTree;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

/// (l, r, a, b): y = a*x + b for x in [l, r).
type Seg = (i64, i64, i64, i64);

fn naive(segs: &[Seg], x: i64, empty: i64) -> i64 {
    segs.iter().filter(|s| s.0 <= x && x < s.1).map(|s| s.2 * x + s.3).min().unwrap_or(empty)
}

#[test]
fn matches_linear_scan() {
    let mut st = Stress::new("li_chao/naive");
    for case in 0..800 * scale() {
        let mut rng = Rng::new(case as u64);
        let (lo, hi) = match case % 3 {
            0 => (0, 0),
            1 => (-20, 20),
            _ => (-1_000_000_000, 1_000_000_000),
        };
        let coef = if hi > 1000 { 1_000_000_000 } else { 10 };
        let empty = i64::MAX;
        let mut tree = LiChaoTree::<i64>::new(lo, hi, empty);
        let mut segs: Vec<Seg> = vec![];
        st.case();
        for step in 0..60 {
            if rng.chance(1, 2) {
                let (a, b) = (rng.range(-coef, coef), rng.range(-coef * 100, coef * 100));
                if rng.chance(1, 2) {
                    tree.add_line(a, b);
                    segs.push((lo, hi + 1, a, b));
                } else {
                    // arbitrary bounds, possibly outside the domain or empty
                    let (l, r) = (rng.range(lo - 3, hi + 3), rng.range(lo - 3, hi + 3));
                    tree.add_segment(l, r, a, b);
                    segs.push((l.max(lo), r.min(hi + 1), a, b));
                }
            } else {
                let x = if rng.chance(1, 4) { *[lo, hi].get(rng.below(2)).unwrap() } else { rng.range(lo, hi) };
                assert_eq!(tree.query(x), naive(&segs, x, empty), "x = {x} case {case} step {step}");
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn wide_values_need_i128() {
    // |a*x + b| reaches ~2e18 + 1e18 here, which still fits i64; i128 must agree anyway.
    let mut t64 = LiChaoTree::<i64>::new(-1_000_000_000, 1_000_000_000, i64::MAX);
    let mut t128 = LiChaoTree::<i128>::new(-1_000_000_000, 1_000_000_000, i128::MAX);
    let mut rng = Rng::new(3);
    for _ in 0..500 {
        let (a, b) = (rng.range(-1_000_000_000, 1_000_000_000), rng.range(-1_000_000_000_000_000_000, 1_000_000_000_000_000_000));
        t64.add_line(a, b);
        t128.add_line(a as i128, b as i128);
    }
    for _ in 0..2000 {
        let x = rng.range(-1_000_000_000, 1_000_000_000);
        assert_eq!(t64.query(x) as i128, t128.query(x));
    }
}

#[test]
fn query_outside_the_domain_panics() {
    let t = LiChaoTree::<i64>::new(0, 10, i64::MAX);
    assert_eq!(t.query(5), i64::MAX); // nothing added: the `empty` value
    assert!(catch_unwind(AssertUnwindSafe(|| t.query(11))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.query(-1))).is_err());
}

#[test]
fn many_lines_and_segments_at_scale() {
    // 200k insertions over a 2e9-wide domain, checked against the naive scan on a window of the
    // inserted lines (queries are restricted to the points where only that window matters).
    let mut st = Stress::new("li_chao/large");
    let mut rng = Rng::new(21);
    let (lo, hi) = (-1_000_000_000i64, 1_000_000_000i64);
    let mut tree = LiChaoTree::<i64>::new(lo, hi, i64::MAX);
    let mut segs: Vec<Seg> = vec![];
    st.case();
    for _ in 0..200_000 {
        let (a, b) = (rng.range(-1_000_000_000, 1_000_000_000), rng.range(-1_000_000_000_000_000_000, 1_000_000_000_000_000_000));
        let (l, r) = (rng.range(lo, hi), rng.range(lo, hi + 1));
        tree.add_segment(l, r, a, b);
        segs.push((l, r, a, b));
        st.ops(1);
    }
    for _ in 0..300 * scale() {
        let x = rng.range(lo, hi);
        assert_eq!(tree.query(x), naive(&segs, x, i64::MAX), "x = {x}");
    }
    st.done();
}
