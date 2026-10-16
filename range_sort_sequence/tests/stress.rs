//! Stress test: `RangeSortSequence` against a plain `Vec<(key, value)>` with real sorting,
//! using the non-commutative affine-composition monoid.

use algebra::Monoid;
use range_sort_sequence::RangeSortSequence;
use std::collections::HashSet;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const P: u64 = 998_244_353;

/// Composition `op(f, g) = g after f` of maps `x -> a*x + b`.
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

fn fresh_key(rng: &mut Rng, used: &mut HashSet<u64>, bits: u32) -> u64 {
    loop {
        let k = rng.next_u64() >> (64 - bits);
        if used.insert(k) {
            return k;
        }
    }
}

fn run_case(rng: &mut Rng, n: usize, bits: u32, steps: usize, ctx: &str) {
    let mut used = HashSet::new();
    let mut model: Vec<(u64, (u64, u64))> = (0..n)
        .map(|_| (fresh_key(rng, &mut used, bits), (rng.below(P as usize) as u64, rng.below(P as usize) as u64)))
        .collect();
    let mut s = RangeSortSequence::<Compose>::new(bits, &model);
    for step in 0..steps {
        let ctx = format!("{ctx} step {step}");
        let l = rng.below(n + 1);
        let r = l + rng.below(n + 1 - l);
        match rng.below(5) {
            0 => {
                s.sort_ascending(l, r);
                model[l..r].sort_by_key(|it| it.0);
            }
            1 => {
                s.sort_descending(l, r);
                model[l..r].sort_by_key(|it| std::cmp::Reverse(it.0));
            }
            2 if n > 0 => {
                let i = rng.below(n);
                used.remove(&model[i].0);
                let key = fresh_key(rng, &mut used, bits);
                let value = (rng.below(P as usize) as u64, rng.below(P as usize) as u64);
                s.set(i, key, value);
                model[i] = (key, value);
            }
            _ => {
                let want = model[l..r].iter().fold((1, 0), |acc, it| Compose::op(&acc, &it.1));
                assert_eq!(s.prod(l, r), want, "prod({l}, {r}) {ctx}");
            }
        }
    }
    let want = model.iter().fold((1, 0), |acc, it| Compose::op(&acc, &it.1));
    assert_eq!(s.prod(0, n), want, "final fold {ctx}");
    for i in 0..n {
        assert_eq!(s.prod(i, i + 1), model[i].1, "final item {i} {ctx}");
    }
}

#[test]
fn matches_vec_model() {
    let mut st = Stress::new("range_sort_sequence/naive");
    for case in 0..1000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(0, 30) as usize;
        // Few bits make keys dense (many splits land on shared tree nodes).
        let bits = (rng.range(5, 12) as u32).max((usize::BITS - n.leading_zeros()) + 1);
        run_case(&mut rng, n, bits, 100, &format!("case {case} n={n} bits={bits}"));
        st.case();
        st.ops(100);
    }
    st.done();
}

#[test]
fn large_workload() {
    let mut st = Stress::new("range_sort_sequence/large");
    let (n, q) = (100_000, 100_000 * scale());
    let mut rng = Rng::new(31);
    let mut used = HashSet::new();
    let items: Vec<(u64, (u64, u64))> = (0..n)
        .map(|_| (fresh_key(&mut rng, &mut used, 30), (rng.below(P as usize) as u64, rng.below(P as usize) as u64)))
        .collect();
    let mut s = RangeSortSequence::<Compose>::new(30, &items);
    st.case();
    let mut acc = 0u64;
    for _ in 0..q {
        let l = rng.below(n);
        let r = l + 1 + rng.below(n - l);
        match rng.below(4) {
            0 => s.sort_ascending(l, r),
            1 => s.sort_descending(l, r),
            2 => {
                let i = rng.below(n);
                let key = fresh_key(&mut rng, &mut used, 30);
                s.set(i, key, (rng.below(P as usize) as u64, rng.below(P as usize) as u64));
            }
            _ => acc ^= s.prod(l, r).0,
        }
        st.ops(1);
    }
    let _ = acc;
    st.done();
}

#[test]
fn panics_on_bad_arguments() {
    let items = [(1u64, (2u64, 3u64)), (2, (4, 5))];
    let mut s = RangeSortSequence::<Compose>::new(4, &items);
    assert!(catch_unwind(AssertUnwindSafe(|| s.sort_ascending(0, 3))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| s.prod(2, 1))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| s.set(2, 3, (1, 1)))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| s.set(0, 16, (1, 1)))).is_err(), "key does not fit");
    assert!(catch_unwind(AssertUnwindSafe(|| RangeSortSequence::<Compose>::new(4, &[(3, (1, 0)), (3, (1, 0))]))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| RangeSortSequence::<Compose>::new(0, &items))).is_err());
    assert_eq!(s.prod(1, 1), (1, 0));
    // Setting an item to the key of another item and merging must be caught.
    let mut t = RangeSortSequence::<Compose>::new(4, &items);
    t.set(0, 2, (1, 1));
    assert!(catch_unwind(AssertUnwindSafe(|| t.sort_ascending(0, 2))).is_err());
}
