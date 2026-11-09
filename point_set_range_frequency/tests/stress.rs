//! Stress test: `PointSetRangeFrequency` against a plain array and linear counting.

use point_set_range_frequency::PointSetRangeFrequency;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn range(rng: &mut Rng, n: usize) -> (usize, usize) {
    let l = rng.below(n + 1);
    (l, l + rng.below(n + 1 - l))
}

#[test]
fn matches_naive_counting() {
    let mut st = Stress::new("point_set_range_frequency/naive");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 40) as usize;
        let alphabet = rng.range(1, 6);
        let mut a: Vec<i64> = (0..n).map(|_| rng.range(-1, alphabet)).collect();
        let updates: Vec<(usize, i64)> = (0..rng.range(0, 60)).map(|_| (rng.below(n), rng.range(-3, alphabet + 3))).collect();
        let mut f = PointSetRangeFrequency::new(&a, &updates);
        st.case();
        let mut next_update = 0;
        for step in 0..100 {
            let ctx = format!("case {case} step {step}");
            if next_update < updates.len() && rng.chance(1, 2) {
                let (p, v) = updates[next_update];
                next_update += 1;
                f.set(p, v);
                a[p] = v;
                assert_eq!(f.get(p), v, "get {ctx}");
            } else {
                let (l, r) = range(&mut rng, n);
                let x = rng.range(-4, alphabet + 4);
                let want = a[l..r].iter().filter(|&&y| y == x).count();
                assert_eq!(f.count(l, r, x), want, "count({l}, {r}, {x}) {ctx}");
            }
            st.ops(1);
        }
        // Re-assigning the current value (declared or not) is a no-op.
        let cur = a[0];
        f.set(0, cur);
        assert_eq!(f.count(0, n, cur), a.iter().filter(|&&y| y == cur).count());
    }
    st.done();
}

#[test]
fn large_scale() {
    let mut st = Stress::new("point_set_range_frequency/large");
    let n = 200_000;
    let q = 200_000 * scale();
    let mut rng = Rng::new(3);
    let mut a: Vec<u32> = (0..n).map(|_| rng.below(1000) as u32).collect();
    let updates: Vec<(usize, u32)> = (0..q).map(|_| (rng.below(n), rng.below(1000) as u32)).collect();
    let mut f = PointSetRangeFrequency::new(&a, &updates);
    st.case();
    for (i, &(p, v)) in updates.iter().enumerate() {
        f.set(p, v);
        a[p] = v;
        if i % 1000 == 0 {
            let (l, r) = range(&mut rng, n);
            let x = rng.below(1000) as u32;
            assert_eq!(f.count(l, r, x), a[l..r].iter().filter(|&&y| y == x).count());
        }
        st.ops(1);
    }
    st.done();
}

#[test]
fn panics_on_bad_arguments() {
    let mut f = PointSetRangeFrequency::new(&[1, 2, 3], &[(0, 9)]);
    assert!(catch_unwind(AssertUnwindSafe(|| f.set(1, 9))).is_err(), "undeclared pair");
    assert!(catch_unwind(AssertUnwindSafe(|| f.set(3, 1))).is_err(), "position");
    assert!(catch_unwind(AssertUnwindSafe(|| f.count(2, 4, 1))).is_err(), "range");
    assert!(catch_unwind(AssertUnwindSafe(|| f.count(2, 1, 1))).is_err(), "range");
    assert!(catch_unwind(AssertUnwindSafe(|| PointSetRangeFrequency::new(&[1], &[(1, 1)]))).is_err());
    let e = PointSetRangeFrequency::<u8>::new(&[], &[]);
    assert_eq!((e.len(), e.count(0, 0, 3)), (0, 0));
}
