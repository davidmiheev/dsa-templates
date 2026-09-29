//! Stress test: `IntervalTree` against a linear scan using the documented overlap rule
//! `lo < qhi && qlo < hi`, plus a counting cross-check at scale.

use interval_tree::{Interval, IntervalTree};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

type Item = (i64, i64, usize);

fn key(it: &Interval<usize>) -> Item {
    (it.lo, it.hi, it.value)
}

fn naive(items: &[Item], qlo: i64, qhi: i64) -> Vec<Item> {
    let mut v: Vec<Item> = items.iter().copied().filter(|&(lo, hi, _)| lo < qhi && qlo < hi).collect();
    v.sort_unstable();
    v
}

fn build(rng: &mut Rng, n: usize, span: i64, max_len: i64) -> (Vec<Item>, IntervalTree<usize>) {
    let items: Vec<Item> = (0..n)
        .map(|id| {
            let lo = rng.range(-span, span);
            (lo, lo + rng.range(0, max_len), id) // max_len >= 0; length 0 gives empty intervals
        })
        .collect();
    let ivs: Vec<Interval<usize>> = items.iter().map(|&(lo, hi, id)| Interval::new(lo, hi, id)).collect();
    (items, IntervalTree::new(&ivs))
}

#[test]
fn matches_linear_scan() {
    let mut st = Stress::new("interval_tree/naive");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.below(60);
        let span = *[5i64, 30, 1000].get(case % 3).unwrap();
        let (items, tree) = build(&mut rng, n, span, span);
        st.case();
        for _ in 0..40 {
            let qlo = rng.range(-span - 2, span + 2);
            let qhi = qlo + rng.range(0, span); // includes empty queries qlo == qhi
            let mut got: Vec<Item> = tree.query(qlo, qhi).iter().map(key).collect();
            got.sort_unstable();
            assert_eq!(got, naive(&items, qlo, qhi), "query [{qlo},{qhi}) case {case} n={n}");
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn large_counts_match_sorted_endpoint_formula() {
    // For a query with qlo < qhi, an interval fails to overlap iff lo >= qhi or hi <= qlo
    // (never both), so #overlaps = n - #(lo >= qhi) - #(hi <= qlo).
    let mut st = Stress::new("interval_tree/large");
    let mut rng = Rng::new(13);
    let n = if cfg!(debug_assertions) { 20_000 } else { 100_000 }; // smaller in debug builds
    let (items, tree) = build(&mut rng, n, 1_000_000_000, 50_000_000);
    let mut los: Vec<i64> = items.iter().map(|i| i.0).collect();
    let mut his: Vec<i64> = items.iter().map(|i| i.1).collect();
    los.sort_unstable();
    his.sort_unstable();
    st.case();
    for q in 0..if cfg!(debug_assertions) { 10_000 } else { 50_000 } * scale() {
        let qlo = rng.range(-1_000_000_000, 1_000_000_000);
        let qhi = qlo + rng.range(1, 100_000_000);
        let hits = tree.query(qlo, qhi);
        let ge_qhi = n - los.partition_point(|&x| x < qhi);
        let le_qlo = his.partition_point(|&x| x <= qlo);
        assert_eq!(hits.len(), n - ge_qhi - le_qlo, "count for [{qlo},{qhi})");
        if q % 250 == 0 {
            let mut got: Vec<Item> = hits.iter().map(key).collect();
            got.sort_unstable();
            assert_eq!(got, naive(&items, qlo, qhi), "exact set for [{qlo},{qhi})");
        }
        st.ops(1);
    }
    st.done();
}

#[test]
fn empty_tree_returns_nothing() {
    let tree: IntervalTree<usize> = IntervalTree::new(&[]);
    assert!(tree.query(-5, 5).is_empty());
}
