//! Stress test: `FastSet` against a `BTreeSet`, for universe sizes around every level boundary.

use fast_set::{FastSet, NONE};
use std::collections::BTreeSet;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn check_all(set: &FastSet, model: &BTreeSet<usize>, n: usize, ctx: &str) {
    for k in 0..n {
        assert_eq!(set.contains(k), model.contains(&k), "contains({k}) {ctx}");
        assert_eq!(set.next(k), model.range(k..).next().copied().unwrap_or(NONE), "next({k}) {ctx}");
        assert_eq!(set.prev(k), model.range(..=k).next_back().copied().unwrap_or(NONE), "prev({k}) {ctx}");
    }
}

#[test]
fn matches_btreeset() {
    let mut st = Stress::new("fast_set/naive");
    // Sizes straddling the 64 and 64^2 boundaries (and tiny ones).
    let sizes = [1usize, 2, 63, 64, 65, 127, 128, 129, 1000, 4095, 4096, 4097, 5000];
    for case in 0..300 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = sizes[case % sizes.len()];
        let mut set = FastSet::new(n);
        let mut model = BTreeSet::new();
        st.case();
        assert_eq!(set.capacity(), n);
        let density = *[1usize, 8, 64].get(rng.below(3)).unwrap(); // few / many elements
        for step in 0..150 {
            let k = rng.below(n);
            let ctx = format!("n={n} case {case} step {step} k={k}");
            match rng.below(4) {
                0 | 1 => {
                    set.insert(k);
                    model.insert(k);
                }
                2 => {
                    if rng.below(density.max(1)) < 2 {
                        set.remove(k);
                        model.remove(&k);
                    }
                }
                _ => {
                    assert_eq!(set.next(k), model.range(k..).next().copied().unwrap_or(NONE), "next {ctx}");
                    assert_eq!(set.prev(k), model.range(..=k).next_back().copied().unwrap_or(NONE), "prev {ctx}");
                }
            }
            st.ops(1);
        }
        if n <= 5000 {
            check_all(&set, &model, n.min(300), &format!("final n={n} case {case}"));
        }
    }
    st.done();
}

#[test]
fn bounds_are_checked() {
    let mut s = FastSet::new(100);
    assert!(catch_unwind(AssertUnwindSafe(|| s.insert(100))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| s.prev(100))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| s.contains(100))).is_err());
    assert_eq!(s.next(100), NONE); // asking past the end is simply "no successor"
    assert_eq!(FastSet::new(0).next(0), NONE);
}

#[test]
fn ten_million_elements() {
    let mut st = Stress::new("fast_set/large");
    let n = 10_000_000;
    let mut rng = Rng::new(8);
    let mut set = FastSet::new(n);
    let mut sorted: Vec<usize> = vec![];
    st.case();
    for _ in 0..200_000 {
        let k = rng.below(n);
        set.insert(k);
        sorted.push(k);
    }
    sorted.sort_unstable();
    sorted.dedup();
    for _ in 0..200_000 * scale() {
        let k = rng.below(n);
        let i = sorted.partition_point(|&v| v < k);
        assert_eq!(set.next(k), sorted.get(i).copied().unwrap_or(NONE));
        let j = sorted.partition_point(|&v| v <= k);
        assert_eq!(set.prev(k), if j == 0 { NONE } else { sorted[j - 1] });
        st.ops(2);
    }
    st.done();
}
