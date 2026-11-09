//! Stress test: `ImplicitTreap` against a plain `Vec` under insert / remove / reverse /
//! range-affine / range-sum, plus a non-commutative-action check without reversal.

use algebra::{MapMonoid, Monoid};
use implicit_treap::ImplicitTreap;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const P: u64 = 998_244_353;

/// (sum, length) under range-affine actions `x -> a*x + b`.
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

fn elem(x: u64) -> (u64, u64) {
    (x, 1)
}

fn range(rng: &mut Rng, n: usize) -> (usize, usize) {
    let l = rng.below(n + 1);
    (l, l + rng.below(n + 1 - l))
}

#[test]
fn matches_vec_model() {
    let mut st = Stress::new("implicit_treap/naive");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(0, 40) as usize;
        let mut model: Vec<u64> = (0..n).map(|_| rng.below(P as usize) as u64).collect();
        let init: Vec<(u64, u64)> = model.iter().map(|&x| elem(x)).collect();
        let mut t = if case % 2 == 0 { ImplicitTreap::<Affine>::from_slice(&init) } else {
            let mut t = ImplicitTreap::<Affine>::new();
            for (i, e) in init.iter().enumerate() { t.insert(i, *e); }
            t
        };
        st.case();
        for step in 0..120 {
            let ctx = format!("case {case} step {step} len={}", model.len());
            let len = model.len();
            match rng.below(7) {
                0 | 1 => {
                    let (i, x) = (rng.below(len + 1), rng.below(P as usize) as u64);
                    t.insert(i, elem(x));
                    model.insert(i, x);
                }
                2 if len > 0 => {
                    let i = rng.below(len);
                    assert_eq!(t.remove(i), elem(model.remove(i)), "remove({i}) {ctx}");
                }
                3 => {
                    let (l, r) = range(&mut rng, len);
                    t.reverse(l, r);
                    model[l..r].reverse();
                }
                4 => {
                    let (l, r) = range(&mut rng, len);
                    let (a, b) = (rng.below(P as usize) as u64, rng.below(P as usize) as u64);
                    t.apply(l, r, (a, b));
                    model[l..r].iter_mut().for_each(|x| *x = (a * *x + b) % P);
                }
                5 => {
                    let (l, r) = range(&mut rng, len);
                    let sum = model[l..r].iter().fold(0, |s, &x| (s + x) % P);
                    assert_eq!(t.prod(l, r), (sum, (r - l) as u64), "prod({l},{r}) {ctx}");
                }
                _ if len > 0 => {
                    let i = rng.below(len);
                    assert_eq!(t.get(i), elem(model[i]), "get({i}) {ctx}");
                }
                _ => {}
            }
            assert_eq!((t.len(), t.is_empty()), (model.len(), model.is_empty()), "len {ctx}");
            st.ops(1);
        }
        assert_eq!(t.to_vec(), model.iter().map(|&x| elem(x)).collect::<Vec<_>>(), "final contents case {case}");
        st.ops(model.len());
    }
    st.done();
}

/// Without reversal the monoid may be non-commutative: here function composition.
enum Compose {}
impl Monoid for Compose {
    type S = (u64, u64); // x -> a*x + b; op = apply left first, then right
    fn identity() -> (u64, u64) {
        (1, 0)
    }
    fn op(f: &(u64, u64), g: &(u64, u64)) -> (u64, u64) {
        (g.0 * f.0 % P, (g.0 * f.1 + g.1) % P)
    }
}
/// No-op action, so only the order-sensitive aggregate is exercised.
enum NoAction {}
impl MapMonoid for NoAction {
    type M = Compose;
    type F = ();
    fn id_map() {}
    fn apply(_: &(), x: &(u64, u64)) -> (u64, u64) {
        *x
    }
    fn compose(_: &(), _: &()) {}
}

#[test]
fn order_sensitive_aggregate_without_reverse() {
    let mut st = Stress::new("implicit_treap/compose");
    for case in 0..300 * scale() {
        let mut rng = Rng::new(7_000 + case as u64);
        let mut model: Vec<(u64, u64)> = vec![];
        let mut t = ImplicitTreap::<NoAction>::new();
        st.case();
        for step in 0..100 {
            let len = model.len();
            if len == 0 || rng.chance(1, 2) {
                let (i, f) = (rng.below(len + 1), (rng.below(P as usize - 1) as u64 + 1, rng.below(P as usize) as u64));
                t.insert(i, f);
                model.insert(i, f);
            } else if rng.chance(1, 4) {
                let i = rng.below(len);
                assert_eq!(t.remove(i), model.remove(i));
            } else {
                let (l, r) = range(&mut rng, len);
                let want = model[l..r].iter().fold(Compose::identity(), |a, f| Compose::op(&a, f));
                assert_eq!(t.prod(l, r), want, "prod({l},{r}) case {case} step {step}");
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn invalid_arguments_panic() {
    let mut t = ImplicitTreap::<Affine>::from_slice(&[elem(1), elem(2)]);
    assert!(catch_unwind(AssertUnwindSafe(|| t.insert(3, elem(0)))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.remove(2))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.get(2))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.reverse(1, 3))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.prod(2, 1))).is_err());
    assert_eq!(t.len(), 2); // the failed calls left the sequence intact
    assert_eq!(t.prod(0, 0), (0, 0));
}

#[test]
fn large_sequence_operations() {
    let mut st = Stress::new("implicit_treap/large");
    let mut rng = Rng::new(31);
    let n = 100_000;
    let mut model: Vec<u64> = (0..n).map(|_| rng.below(P as usize) as u64).collect();
    let mut t = ImplicitTreap::<Affine>::from_slice(&model.iter().map(|&x| elem(x)).collect::<Vec<_>>());
    st.case();
    for _ in 0..3_000 * scale() {
        let len = model.len();
        match rng.below(5) {
            0 => {
                let (i, x) = (rng.below(len + 1), rng.below(P as usize) as u64);
                t.insert(i, elem(x));
                model.insert(i, x);
            }
            1 => {
                let i = rng.below(len);
                assert_eq!(t.remove(i), elem(model.remove(i)));
            }
            2 => {
                let (l, r) = range(&mut rng, len);
                t.reverse(l, r);
                model[l..r].reverse();
            }
            3 => {
                let (l, r) = range(&mut rng, len);
                let (a, b) = (rng.below(P as usize) as u64, rng.below(P as usize) as u64);
                t.apply(l, r, (a, b));
                model[l..r].iter_mut().for_each(|x| *x = (a * *x + b) % P);
            }
            _ => {
                let (l, r) = range(&mut rng, len);
                let sum = model[l..r].iter().fold(0, |s, &x| (s + x) % P);
                assert_eq!(t.prod(l, r), (sum, (r - l) as u64));
            }
        }
        st.ops(1);
    }
    assert_eq!(t.to_vec(), model.iter().map(|&x| elem(x)).collect::<Vec<_>>());
    st.done();
}
