//! Stress test: `SegmentTree` for several associative operations, including a
//! non-commutative one (leftmost non-zero), against a left-to-right fold.

use segment_tree::SegmentTree;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

type Op = fn(i64, i64) -> i64;

fn sum(a: i64, b: i64) -> i64 {
    a + b
}
fn min(a: i64, b: i64) -> i64 {
    a.min(b)
}
fn max(a: i64, b: i64) -> i64 {
    a.max(b)
}
fn xor(a: i64, b: i64) -> i64 {
    a ^ b
}
/// Associative but not commutative: exposes swapped left/right merges.
fn first_nonzero(a: i64, b: i64) -> i64 {
    if a != 0 {
        a
    } else {
        b
    }
}
/// (name, op, identity)
const OPS: [(&str, Op, i64); 5] = [
    ("sum", sum, 0),
    ("min", min, i64::MAX),
    ("max", max, i64::MIN),
    ("xor", xor, 0),
    ("first_nonzero", first_nonzero, 0),
];

fn fold(model: &[i64], l: usize, r: usize, op: Op, id: i64) -> i64 {
    model[l..=r].iter().fold(id, |acc, &x| op(acc, x))
}

#[test]
fn matches_naive_fold() {
    let mut st = Stress::new("segment_tree/naive");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(case as u64);
        let (name, op, id) = OPS[case % OPS.len()];
        let n = rng.range(1, 70) as usize;
        let small = rng.chance(1, 2); // small values make zeros/duplicates common
        let val = |rng: &mut Rng| if small { rng.range(-3, 3) } else { rng.range(-1_000_000, 1_000_000) };
        let mut model: Vec<i64> = (0..n).map(|_| val(&mut rng)).collect();
        let mut tree = SegmentTree::new(&model, id, op);
        st.case();
        for step in 0..150 {
            if rng.chance(1, 3) {
                let (i, v) = (rng.below(n), val(&mut rng));
                model[i] = v;
                tree.update(i, v);
            } else {
                let l = rng.below(n);
                let r = l + rng.below(n - l);
                assert_eq!(
                    tree.query(l, r),
                    fold(&model, l, r, op, id),
                    "{name} query({l},{r}) case {case} step {step} n={n}"
                );
            }
            st.ops(1);
        }
        // Whole-array and single-element queries.
        assert_eq!(tree.query(0, n - 1), fold(&model, 0, n - 1, op, id), "{name} full case {case}");
        for i in 0..n {
            assert_eq!(tree.query(i, i), model[i], "{name} point case {case}");
        }
        st.ops(n + 1);
    }
    st.done();
}

#[test]
fn large_random_ranges() {
    let mut st = Stress::new("segment_tree/large");
    let n = 200_000;
    let mut rng = Rng::new(9);
    let mut model: Vec<i64> = (0..n).map(|_| rng.range(-1_000_000, 1_000_000)).collect();
    let mut tree = SegmentTree::new(&model, i64::MIN, max);
    let mut sum_tree = SegmentTree::new(&model, 0, sum);
    let mut total: i64 = model.iter().sum();
    st.case();
    for _ in 0..100_000 * scale() {
        if rng.chance(1, 2) {
            let (i, v) = (rng.below(n), rng.range(-1_000_000, 1_000_000));
            total += v - model[i];
            model[i] = v;
            tree.update(i, v);
            sum_tree.update(i, v);
        } else {
            let l = rng.below(n);
            let r = (l + rng.below(300)).min(n - 1); // short ranges keep the naive fold cheap
            assert_eq!(tree.query(l, r), fold(&model, l, r, max, i64::MIN));
            assert_eq!(sum_tree.query(l, r), fold(&model, l, r, sum, 0));
            assert_eq!(sum_tree.query(0, n - 1), total);
        }
        st.ops(1);
    }
    st.done();
}
