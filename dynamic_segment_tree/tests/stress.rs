//! Stress test: dynamic segment trees against naive models, on small arrays (full `Vec`) and on
//! huge arrays (a few distinct coordinates, so the model works on elementary intervals).

use algebra::{MapMonoid, Monoid};
use dynamic_segment_tree::{DynamicLazySegmentTree, DynamicSegmentTree};
use std::collections::BTreeMap;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const P: u64 = 998_244_353;

/// Composition of affine maps `x -> a*x + b`; `op(f, g)` applies `f` first (non-commutative).
enum Compose {}
impl Monoid for Compose {
    type S = (u64, u64);
    fn identity() -> (u64, u64) {
        (1, 0)
    }
    fn op(f: &(u64, u64), g: &(u64, u64)) -> (u64, u64) {
        (f.0 * g.0 % P, (g.0 * f.1 + g.1) % P)
    }
}

/// (sum, length) under range-affine actions.
enum SumLen {}
impl Monoid for SumLen {
    type S = (u64, u64);
    fn identity() -> (u64, u64) {
        (0, 0)
    }
    fn op(a: &(u64, u64), b: &(u64, u64)) -> (u64, u64) {
        ((a.0 + b.0) % P, a.1 + b.1)
    }
}
enum Affine {}
impl MapMonoid for Affine {
    type M = SumLen;
    type F = (u64, u64);
    fn id_map() -> (u64, u64) {
        (1, 0)
    }
    fn apply(f: &(u64, u64), x: &(u64, u64)) -> (u64, u64) {
        ((f.0 * x.0 + f.1 * x.1) % P, x.1)
    }
    fn compose(f: &(u64, u64), g: &(u64, u64)) -> (u64, u64) {
        (f.0 * g.0 % P, (f.0 * g.1 + f.1) % P)
    }
}

fn range(rng: &mut Rng, n: usize) -> (usize, usize) {
    let l = rng.below(n + 1);
    (l, l + rng.below(n + 1 - l))
}

#[test]
fn point_tree_matches_vec() {
    let mut st = Stress::new("dynamic_segment_tree/point-naive");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 70) as usize;
        let mut t = DynamicSegmentTree::<Compose>::new(n, (1, 0));
        let mut model = vec![(1u64, 0u64); n];
        st.case();
        for step in 0..100 {
            let ctx = format!("case {case} step {step} n={n}");
            if rng.chance(1, 2) {
                let i = rng.below(n);
                let f = (rng.below(P as usize) as u64, rng.below(P as usize) as u64);
                t.set(i, f);
                model[i] = f;
            } else {
                let (l, r) = range(&mut rng, n);
                let want = model[l..r].iter().fold((1, 0), |acc, g| Compose::op(&acc, g));
                assert_eq!(t.prod(l, r), want, "prod({l}, {r}) {ctx}");
                let i = rng.below(n);
                assert_eq!(t.get(i), model[i], "get({i}) {ctx}");
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn point_tree_on_a_huge_array() {
    let mut st = Stress::new("dynamic_segment_tree/point-huge");
    let n = 1_000_000_000usize;
    for case in 0..100 * scale() {
        let mut rng = Rng::new(case as u64 + 7);
        let mut t = DynamicSegmentTree::<Compose>::new(n, (1, 0));
        let mut model: BTreeMap<usize, (u64, u64)> = BTreeMap::new();
        let pool: Vec<usize> = (0..15).map(|_| rng.below(n)).collect();
        st.case();
        for step in 0..150 {
            let ctx = format!("case {case} step {step}");
            if rng.chance(1, 2) {
                let i = pool[rng.below(pool.len())];
                let f = (rng.below(P as usize) as u64, rng.below(P as usize) as u64);
                t.set(i, f);
                model.insert(i, f);
            } else {
                let (mut l, mut r) = (rng.below(n + 1), rng.below(n + 1));
                if rng.chance(1, 2) {
                    l = pool[rng.below(pool.len())];
                    r = (pool[rng.below(pool.len())] + 1).min(n);
                }
                if l > r {
                    std::mem::swap(&mut l, &mut r);
                }
                let want = model.range(l..r).fold((1, 0), |acc, (_, g)| Compose::op(&acc, g));
                assert_eq!(t.prod(l, r), want, "prod({l}, {r}) {ctx}");
            }
            st.ops(1);
        }
        assert!(t.node_count() <= 150 * 32 + 1, "nodes: {}", t.node_count());
    }
    st.done();
}

#[test]
fn lazy_tree_matches_vec() {
    let mut st = Stress::new("dynamic_lazy_segment_tree/naive");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64 + 100);
        let n = rng.range(1, 70) as usize;
        let init = rng.below(P as usize) as u64;
        let mut t = DynamicLazySegmentTree::<Affine>::new(n, (init, 1));
        let mut model = vec![init; n];
        st.case();
        for step in 0..100 {
            let ctx = format!("case {case} step {step} n={n}");
            match rng.below(4) {
                0 | 1 => {
                    let (l, r) = range(&mut rng, n);
                    let (a, b) = (rng.below(P as usize) as u64, rng.below(P as usize) as u64);
                    t.apply_range(l, r, (a, b));
                    for x in &mut model[l..r] {
                        *x = (a * *x + b) % P;
                    }
                }
                2 => {
                    let (i, x) = (rng.below(n), rng.below(P as usize) as u64);
                    t.set(i, (x, 1));
                    model[i] = x;
                }
                _ => {
                    let (l, r) = range(&mut rng, n);
                    let want = (model[l..r].iter().sum::<u64>() % P, (r - l) as u64);
                    assert_eq!(t.prod(l, r), want, "prod({l}, {r}) {ctx}");
                    let i = rng.below(n);
                    assert_eq!(t.get(i), (model[i], 1), "get({i}) {ctx}");
                }
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn lazy_tree_on_a_huge_array() {
    let mut st = Stress::new("dynamic_lazy_segment_tree/huge");
    let n = 1_000_000_000usize;
    for case in 0..100 * scale() {
        let mut rng = Rng::new(case as u64 + 9);
        // Boundaries: 0, n, and a few random points; points `p` also contribute `p + 1`, so that
        // a one-element interval exists for `set`.
        let mut cuts = vec![0, n];
        for _ in 0..8 {
            let p = rng.below(n);
            cuts.push(p);
            cuts.push(p + 1);
        }
        cuts.sort();
        cuts.dedup();
        let m = cuts.len() - 1;
        let mut t = DynamicLazySegmentTree::<Affine>::new(n, (0, 1));
        let mut value = vec![0u64; m]; // every element of interval j holds value[j]
        let len = |j: usize| (cuts[j + 1] - cuts[j]) as u64;
        st.case();
        for step in 0..120 {
            let ctx = format!("case {case} step {step}");
            let (mut a, mut b) = (rng.below(m + 1), rng.below(m + 1));
            if a > b {
                std::mem::swap(&mut a, &mut b);
            }
            let (l, r) = (cuts[a], cuts[b]);
            match rng.below(4) {
                0 | 1 => {
                    let (x, y) = (rng.below(P as usize) as u64, rng.below(P as usize) as u64);
                    t.apply_range(l, r, (x, y));
                    for v in &mut value[a..b] {
                        *v = (x * *v + y) % P;
                    }
                }
                2 => {
                    let ones: Vec<usize> = (0..m).filter(|&j| len(j) == 1).collect();
                    if let Some(&j) = ones.get(rng.below(ones.len().max(1))) {
                        let x = rng.below(P as usize) as u64;
                        t.set(cuts[j], (x, 1));
                        value[j] = x;
                    }
                }
                _ => {
                    let want = (a..b).fold(0, |s, j| (s + value[j] * (len(j) % P)) % P);
                    assert_eq!(t.prod(l, r), (want, (r - l) as u64), "prod({l}, {r}) {ctx}");
                }
            }
            st.ops(1);
        }
        assert!(t.node_count() <= 120 * 70 + 100, "nodes: {}", t.node_count());
    }
    st.done();
}

#[test]
fn panics_on_bad_arguments() {
    let mut t = DynamicSegmentTree::<Compose>::new(10, (1, 0));
    assert!(catch_unwind(AssertUnwindSafe(|| t.set(10, (1, 0)))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.prod(3, 11))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.prod(4, 3))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| DynamicSegmentTree::<Compose>::new(0, (1, 0)))).is_err());
    assert_eq!(t.prod(5, 5), (1, 0));
    let mut l = DynamicLazySegmentTree::<Affine>::new(10, (0, 1));
    assert!(catch_unwind(AssertUnwindSafe(|| l.apply_range(0, 11, (1, 1)))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| l.get(10))).is_err());
    assert_eq!(l.prod(7, 7), (0, 0));
    assert_eq!(l.prod(0, 10), (0, 10));
    assert_eq!(l.node_count(), 0, "reads do not allocate");
}
