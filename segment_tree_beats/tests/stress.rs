//! Stress test: `SegmentTreeBeats` against a plain `Vec<i64>` under random chmin / chmax /
//! add / sum / min / max.

use segment_tree_beats::SegmentTreeBeats;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn range(rng: &mut Rng, n: usize) -> (usize, usize) {
    let l = rng.below(n + 1);
    (l, l + rng.below(n + 1 - l))
}

fn run_case(rng: &mut Rng, n: usize, vmax: i64, steps: usize, ctx: &str) {
    let mut model: Vec<i64> = (0..n).map(|_| rng.range(-vmax, vmax)).collect();
    let mut t = SegmentTreeBeats::new(&model);
    for step in 0..steps {
        let ctx = format!("{ctx} step {step}");
        let (l, r) = range(rng, n);
        match rng.below(6) {
            0 => {
                let x = rng.range(-vmax, vmax);
                t.chmin(l, r, x);
                model[l..r].iter_mut().for_each(|v| *v = (*v).min(x));
            }
            1 => {
                let x = rng.range(-vmax, vmax);
                t.chmax(l, r, x);
                model[l..r].iter_mut().for_each(|v| *v = (*v).max(x));
            }
            2 => {
                let x = rng.range(-vmax, vmax);
                t.add(l, r, x);
                model[l..r].iter_mut().for_each(|v| *v += x);
            }
            3 => assert_eq!(t.sum(l, r), model[l..r].iter().sum::<i64>(), "sum({l}, {r}) {ctx}"),
            4 if l < r => {
                assert_eq!(t.min(l, r), *model[l..r].iter().min().unwrap(), "min({l}, {r}) {ctx}");
                assert_eq!(t.max(l, r), *model[l..r].iter().max().unwrap(), "max({l}, {r}) {ctx}");
            }
            _ => {
                let i = rng.below(n);
                assert_eq!(t.get(i), model[i], "get({i}) {ctx}");
            }
        }
    }
    for (i, &v) in model.iter().enumerate() {
        assert_eq!(t.get(i), v, "final get({i}) {ctx}");
    }
}

#[test]
fn matches_vec_model() {
    let mut st = Stress::new("segment_tree_beats/naive");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 40) as usize;
        let vmax = *[2i64, 10, 1000, 1_000_000_000_000].get(case % 4).unwrap();
        run_case(&mut rng, n, vmax, 120, &format!("case {case} n={n}"));
        st.case();
        st.ops(120);
    }
    st.done();
}

#[test]
fn large_random_workload() {
    let mut st = Stress::new("segment_tree_beats/large");
    let (n, q) = (200_000, 200_000 * scale());
    let mut rng = Rng::new(17);
    let a: Vec<i64> = (0..n).map(|_| rng.range(-1_000_000_000, 1_000_000_000)).collect();
    let mut t = SegmentTreeBeats::new(&a);
    let mut check_sum = 0i64;
    st.case();
    for _ in 0..q {
        let (l, r) = range(&mut rng, n);
        let x = rng.range(-1_000_000_000, 1_000_000_000);
        match rng.below(4) {
            0 => t.chmin(l, r, x),
            1 => t.chmax(l, r, x),
            2 => t.add(l, r, x % 1000),
            _ => check_sum = check_sum.wrapping_add(t.sum(l, r)),
        }
        st.ops(1);
    }
    let _ = check_sum;
    st.done();
}

#[test]
fn panics_on_bad_arguments() {
    let mut t = SegmentTreeBeats::new(&[1, 2, 3]);
    assert!(catch_unwind(AssertUnwindSafe(|| t.chmin(0, 4, 1))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.add(2, 1, 1))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.sum(0, 4))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.min(1, 1))).is_err(), "min of an empty range");
    assert!(catch_unwind(AssertUnwindSafe(|| t.get(3))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| SegmentTreeBeats::new(&[]))).is_err());
    assert_eq!(t.sum(2, 2), 0);
    t.chmax(1, 1, 100); // empty range: no-op
    assert_eq!(t.sum(0, 3), 6);
}
