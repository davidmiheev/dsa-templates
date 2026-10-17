//! Stress test: `QueueAggregation` and `DequeAggregation` against a `VecDeque` and a fresh
//! left-to-right fold, with a non-commutative monoid (affine maps mod p).

use algebra::Monoid;
use queue_aggregation::{DequeAggregation, QueueAggregation};
use std::collections::VecDeque;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const P: u64 = 998_244_353;
/// `x -> a*x + b`; `op(f, g)` = apply `f` first, then `g` (not commutative).
enum Affine {}
impl Monoid for Affine {
    type S = (u64, u64);
    fn identity() -> (u64, u64) {
        (1, 0)
    }
    fn op(f: &(u64, u64), g: &(u64, u64)) -> (u64, u64) {
        (g.0 * f.0 % P, (g.0 * f.1 + g.1) % P)
    }
}

fn naive_fold(items: &VecDeque<(u64, u64)>) -> (u64, u64) {
    items.iter().fold(Affine::identity(), |acc, f| Affine::op(&acc, f))
}

fn rand_fn(rng: &mut Rng) -> (u64, u64) {
    (rng.below(P as usize - 1) as u64 + 1, rng.below(P as usize) as u64)
}

#[test]
fn queue_matches_naive() {
    let mut st = Stress::new("queue_agg/queue");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let mut q = QueueAggregation::<Affine>::new();
        let mut model: VecDeque<(u64, u64)> = VecDeque::new();
        st.case();
        for step in 0..200 {
            if model.is_empty() || rng.chance(3, 5) {
                let f = rand_fn(&mut rng);
                q.push_back(f);
                model.push_back(f);
            } else {
                assert_eq!(q.pop_front(), model.pop_front().unwrap(), "pop_front case {case} step {step}");
            }
            assert_eq!(q.fold(), naive_fold(&model), "fold case {case} step {step}");
            assert_eq!((q.len(), q.is_empty()), (model.len(), model.is_empty()));
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn deque_matches_naive() {
    let mut st = Stress::new("queue_agg/deque");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(10_000 + case as u64);
        let mut d = DequeAggregation::<Affine>::new();
        let mut model: VecDeque<(u64, u64)> = VecDeque::new();
        st.case();
        // Phases favouring one end make the rebalancing paths (pop from an empty stack) common.
        for step in 0..300 {
            let phase = (step / 50) % 3;
            let op = if model.is_empty() { rng.below(2) } else { rng.below(4) };
            let ctx = format!("case {case} step {step}");
            match (op, phase) {
                (0, _) => {
                    let f = rand_fn(&mut rng);
                    d.push_front(f);
                    model.push_front(f);
                }
                (1, _) => {
                    let f = rand_fn(&mut rng);
                    d.push_back(f);
                    model.push_back(f);
                }
                (2, 1) | (3, 2) | (2, 0) => assert_eq!(d.pop_front(), model.pop_front().unwrap(), "pop_front {ctx}"),
                _ => assert_eq!(d.pop_back(), model.pop_back().unwrap(), "pop_back {ctx}"),
            }
            assert_eq!(d.fold(), naive_fold(&model), "fold {ctx}");
            assert_eq!((d.len(), d.is_empty()), (model.len(), model.is_empty()), "len {ctx}");
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn popping_from_empty_panics() {
    let mut q = QueueAggregation::<Affine>::new();
    let mut d = DequeAggregation::<Affine>::new();
    assert_eq!((q.fold(), d.fold()), (Affine::identity(), Affine::identity()));
    assert!(catch_unwind(AssertUnwindSafe(|| q.pop_front())).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| d.pop_front())).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| d.pop_back())).is_err());
}

#[test]
fn large_amortized_workloads() {
    // 1M operations each; a pop-heavy deque alternating ends exercises the halving rebuilds.
    let mut st = Stress::new("queue_agg/large");
    let mut rng = Rng::new(2);
    let mut d = DequeAggregation::<Affine>::new();
    let mut model: VecDeque<(u64, u64)> = VecDeque::new();
    st.case();
    for i in 0..1_000_000usize {
        if model.len() < 2 || rng.chance(1, 2) {
            let f = rand_fn(&mut rng);
            if rng.chance(1, 2) { d.push_front(f); model.push_front(f) } else { d.push_back(f); model.push_back(f) }
        } else if rng.chance(1, 2) {
            assert_eq!(d.pop_front(), model.pop_front().unwrap());
        } else {
            assert_eq!(d.pop_back(), model.pop_back().unwrap());
        }
        if i % 5000 == 0 {
            assert_eq!(d.fold(), naive_fold(&model));
        }
        st.ops(1);
    }
    st.done();
}
