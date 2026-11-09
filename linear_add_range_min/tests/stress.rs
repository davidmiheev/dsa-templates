//! Stress test: `LinearAddRangeMin` against a plain `Vec<i64>`.

use linear_add_range_min::LinearAddRangeMin;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_vec_model() {
    let mut st = Stress::new("linear_add_range_min/naive");
    for case in 0..1000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 50) as usize;
        let block = rng.range(1, 12) as usize;
        let amp = *[3i64, 100, 1_000_000_000].get(case % 3).unwrap();
        let mut model: Vec<i64> = (0..n).map(|_| rng.range(-amp, amp)).collect();
        let mut t = if case % 2 == 0 { LinearAddRangeMin::with_block_size(&model, block) } else { LinearAddRangeMin::new(&model) };
        st.case();
        for step in 0..120 {
            let ctx = format!("case {case} step {step} n={n} block={block}");
            let l = rng.below(n + 1);
            let r = l + rng.below(n + 1 - l);
            if rng.chance(1, 2) {
                let (b, c) = (rng.range(-100, 100), rng.range(-1000, 1000));
                t.add_linear(l, r, b, c);
                for i in l..r {
                    model[i] += b * i as i64 + c;
                }
            } else if l < r {
                assert_eq!(t.min(l, r), *model[l..r].iter().min().unwrap(), "min({l}, {r}) {ctx}");
            } else {
                let i = rng.below(n);
                assert_eq!(t.get(i), model[i], "get({i}) {ctx}");
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn large_workload_with_constraint_sized_values() {
    let mut st = Stress::new("linear_add_range_min/large");
    let (n, q) = (100_000, 100_000 * scale());
    let mut rng = Rng::new(13);
    let a: Vec<i64> = (0..n).map(|_| rng.range(-1_000_000_000, 1_000_000_000)).collect();
    let mut t = LinearAddRangeMin::new(&a);
    let mut model_spot = a.clone();
    st.case();
    for k in 0..q {
        let l = rng.below(n);
        let r = l + 1 + rng.below(n - l);
        if rng.chance(1, 2) {
            let (b, c) = (rng.range(-99, 99), rng.range(-9_999_999, 9_999_999));
            t.add_linear(l, r, b, c);
            if k < 3000 {
                for i in l..r {
                    model_spot[i] += b * i as i64 + c;
                }
            }
        } else {
            let m = t.min(l, r);
            if k < 3000 {
                assert_eq!(m, *model_spot[l..r].iter().min().unwrap(), "min({l}, {r}) at op {k}");
            }
        }
        st.ops(1);
    }
    st.done();
}

#[test]
fn panics_on_bad_arguments() {
    let mut t = LinearAddRangeMin::new(&[1, 2, 3]);
    assert!(catch_unwind(AssertUnwindSafe(|| t.add_linear(0, 4, 1, 1))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.min(1, 1))).is_err(), "empty minimum");
    assert!(catch_unwind(AssertUnwindSafe(|| t.min(0, 4))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.get(3))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| LinearAddRangeMin::new(&[]))).is_err());
    t.add_linear(2, 2, 100, 100); // empty range: no-op
    assert_eq!(t.min(0, 3), 1);
}
