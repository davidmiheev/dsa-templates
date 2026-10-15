//! Stress test: `KdTree` against a plain list of weighted points under rectangle-affine,
//! rectangle-sum, set and get; including inactive points that are activated later.

use algebra::{MapMonoid, Monoid};
use kd_tree::KdTree;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const P: u64 = 998_244_353;

/// (sum of weights of active points, number of active points).
enum SumCount {}
impl Monoid for SumCount {
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
    type M = SumCount;
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

#[test]
fn matches_naive_with_activation() {
    let mut st = Stress::new("kd_tree/naive");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(0, if case % 4 == 0 { 400 } else { 40 }) as usize;
        let c = *[3i64, 10, 1_000_000_000].get(case % 3).unwrap();
        let points: Vec<(i64, i64)> = (0..n).map(|_| (rng.range(0, c), rng.range(0, c))).collect();
        // Weight `None` = inactive point.
        let mut w: Vec<Option<u64>> = (0..n).map(|_| if rng.chance(1, 3) { None } else { Some(rng.below(P as usize) as u64) }).collect();
        let init: Vec<(u64, u64)> = w.iter().map(|x| x.map_or((0, 0), |v| (v, 1))).collect();
        let mut t = KdTree::<Affine, i64, i64>::new(&points, &init);
        assert_eq!(t.len(), n);
        st.case();
        for step in 0..100 {
            let ctx = format!("case {case} step {step} n={n}");
            let (x0, x1) = (rng.range(0, c + 1), rng.range(0, c + 1));
            let (y0, y1) = (rng.range(0, c + 1), rng.range(0, c + 1));
            let (xl, xh, yl, yh) = (x0.min(x1), x0.max(x1), y0.min(y1), y0.max(y1));
            let inside = |p: &(i64, i64)| xl <= p.0 && p.0 < xh && yl <= p.1 && p.1 < yh;
            match rng.below(4) {
                0 => {
                    let f = (rng.below(P as usize) as u64, rng.below(P as usize) as u64);
                    t.apply_rect(xl, xh, yl, yh, f);
                    for (i, p) in points.iter().enumerate() {
                        if inside(p) {
                            w[i] = w[i].map(|v| (f.0 * v + f.1) % P);
                        }
                    }
                }
                1 if n > 0 => {
                    let i = rng.below(n);
                    let v = rng.below(P as usize) as u64;
                    t.set(i, (v, 1));
                    w[i] = Some(v);
                }
                2 if n > 0 => {
                    let i = rng.below(n);
                    assert_eq!(t.get(i), w[i].map_or((0, 0), |v| (v, 1)), "get({i}) {ctx}");
                }
                _ => {
                    let want = points.iter().zip(&w).filter(|(p, _)| inside(p)).fold((0, 0), |acc, (_, x)| match x {
                        Some(v) => ((acc.0 + v) % P, acc.1 + 1),
                        None => acc,
                    });
                    assert_eq!(t.prod_rect(xl, xh, yl, yh), want, "prod_rect({xl}, {xh}, {yl}, {yh}) {ctx}");
                }
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn large_workload() {
    let mut st = Stress::new("kd_tree/large");
    let n = 200_000;
    let mut rng = Rng::new(5);
    let points: Vec<(i64, i64)> = (0..n).map(|_| (rng.range(0, 1_000_000_000), rng.range(0, 1_000_000_000))).collect();
    let init: Vec<(u64, u64)> = (0..n).map(|_| (rng.below(P as usize) as u64, 1)).collect();
    let mut t = KdTree::<Affine, i64, i64>::new(&points, &init);
    let total_before = t.prod_rect(0, 1_000_000_001, 0, 1_000_000_001);
    assert_eq!(total_before.1, n as u64);
    st.case();
    for _ in 0..100_000 * scale() {
        let (x0, x1) = (rng.range(0, 1_000_000_000), rng.range(0, 1_000_000_000));
        let (y0, y1) = (rng.range(0, 1_000_000_000), rng.range(0, 1_000_000_000));
        let (xl, xh, yl, yh) = (x0.min(x1), x0.max(x1) + 1, y0.min(y1), y0.max(y1) + 1);
        if rng.chance(1, 2) {
            t.apply_rect(xl, xh, yl, yh, (rng.below(P as usize) as u64, rng.below(P as usize) as u64));
        } else {
            let _ = t.prod_rect(xl, xh, yl, yh);
        }
        st.ops(1);
    }
    st.done();
}

#[test]
fn panics_on_bad_arguments() {
    let mut t = KdTree::<Affine, i32, i32>::new(&[(1, 1), (2, 2)], &[(1, 1), (2, 1)]);
    assert!(catch_unwind(AssertUnwindSafe(|| t.set(2, (0, 1)))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.get(2))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| KdTree::<Affine, i32, i32>::new(&[(1, 1)], &[]))).is_err());
    assert_eq!(t.prod_rect(5, 2, 0, 9), (0, 0), "inverted rectangle is empty");
    let mut e = KdTree::<Affine, i32, i32>::new(&[], &[]);
    assert_eq!(e.prod_rect(0, 9, 0, 9), (0, 0));
    e.apply_rect(0, 9, 0, 9, (2, 2));
}
