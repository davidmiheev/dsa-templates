//! Stress test: `PersistentQueue` against one full `Vec` per version.

use persistent_queue::PersistentQueue;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_vec_per_version() {
    let mut st = Stress::new("persistent_queue/naive");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let mut q = PersistentQueue::new();
        let mut model: Vec<Vec<i64>> = vec![vec![]]; // version 0 = empty
        st.case();
        for step in 0..150 {
            let v = rng.below(model.len());
            let ctx = format!("case {case} step {step} version {v}");
            if model[v].is_empty() || rng.chance(1, 2) {
                let x = rng.range(-1000, 1000);
                let nv = q.push_back(v, x);
                let mut m = model[v].clone();
                m.push(x);
                assert_eq!(nv, model.len(), "ids are consecutive {ctx}");
                model.push(m);
            } else {
                let (front, nv) = q.pop_front(v);
                assert_eq!(front, model[v][0], "front {ctx}");
                assert_eq!(nv, model.len());
                model.push(model[v][1..].to_vec());
            }
            // Every version stays intact: check a couple of random ones, and the newest.
            for &w in &[model.len() - 1, rng.below(model.len())] {
                assert_eq!(q.len(w), model[w].len(), "len(v{w}) {ctx}");
                assert_eq!(q.is_empty(w), model[w].is_empty());
                if !model[w].is_empty() {
                    assert_eq!(q.front(w), model[w][0], "front(v{w}) {ctx}");
                }
            }
            st.ops(1);
        }
        assert_eq!(q.num_versions(), model.len());
    }
    st.done();
}

#[test]
fn empty_version_has_no_front() {
    let mut q = PersistentQueue::<u8>::new();
    assert!(q.is_empty(PersistentQueue::<u8>::EMPTY));
    assert!(catch_unwind(AssertUnwindSafe(|| q.front(0))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| q.pop_front(0))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| q.len(5))).is_err());
    let v = q.push_back(0, 7);
    assert_eq!((q.pop_front(v).0, q.len(v)), (7, 1)); // popping does not change `v`
}

#[test]
fn deep_histories_stay_logarithmic() {
    // A chain of 500k pushes, then pops from the middle versions: each `front` must be
    // O(log n), so the whole loop is fast; correctness is checked against the position.
    let mut st = Stress::new("persistent_queue/large");
    let n = 500_000;
    let mut q = PersistentQueue::new();
    let mut v = PersistentQueue::<usize>::EMPTY;
    let mut versions = vec![v];
    st.case();
    for i in 0..n {
        v = q.push_back(v, i);
        versions.push(v);
    }
    let mut rng = Rng::new(4);
    for _ in 0..200_000 * scale() {
        let ver = versions[rng.below(versions.len())];
        if q.len(ver) > 0 {
            let (front, nv) = q.pop_front(ver); // not recorded: all recorded versions stay unpopped
            assert_eq!(front, 0, "the queue was never popped: the front is the first push");
            assert_eq!(q.len(nv), q.len(ver) - 1);
        }
        st.ops(1);
    }
    // After many pops chained on one branch the fronts walk forward one by one.
    let mut cur = *versions.get(n).unwrap();
    for expect in 0..1000 {
        let (f, next) = q.pop_front(cur);
        assert_eq!(f, expect);
        cur = next;
    }
    st.done();
}
