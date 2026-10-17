//! Stress test: `MultiSet` against a plain count array.

use multiset::MultiSet;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const V: usize = 12;

fn check(ms: &MultiSet<i32>, counts: &[usize; V], ctx: &str) {
    assert_eq!(ms.len(), counts.iter().sum::<usize>(), "len {ctx}");
    assert_eq!(ms.is_empty(), counts.iter().all(|&c| c == 0), "is_empty {ctx}");
    assert_eq!(ms.first(), (0..V).find(|&v| counts[v] > 0).map(|v| v as i32), "first {ctx}");
    assert_eq!(ms.last(), (0..V).rev().find(|&v| counts[v] > 0).map(|v| v as i32), "last {ctx}");
    for v in 0..V {
        assert_eq!(ms.count(&(v as i32)), counts[v], "count({v}) {ctx}");
        assert_eq!(ms.contains(&(v as i32)), counts[v] > 0, "contains({v}) {ctx}");
    }
    let listed: Vec<(i32, usize)> = ms.iter().map(|(k, c)| (*k, *c)).collect();
    let expect: Vec<(i32, usize)> = (0..V).filter(|&v| counts[v] > 0).map(|v| (v as i32, counts[v])).collect();
    assert_eq!(listed, expect, "iter order/contents {ctx}");
}

#[test]
fn matches_count_array() {
    let mut st = Stress::new("multiset/naive");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(case as u64);
        let mut ms = MultiSet::new();
        let mut counts = [0usize; V];
        st.case();
        for step in 0..150 {
            let v = rng.below(V);
            let ctx = format!("case {case} step {step} v {v}");
            match rng.below(7) {
                0..=2 => {
                    ms.insert(v as i32);
                    counts[v] += 1;
                }
                3 => {
                    let had = counts[v] > 0;
                    assert_eq!(ms.remove(&(v as i32)), had, "remove {ctx}");
                    counts[v] = counts[v].saturating_sub(1);
                }
                4 => {
                    let had = counts[v] > 0;
                    assert_eq!(ms.remove_all(&(v as i32)), had, "remove_all {ctx}");
                    counts[v] = 0;
                }
                5 => {
                    if rng.chance(1, 20) {
                        ms.clear();
                        counts = [0; V];
                    }
                }
                _ => {
                    for (_, c) in ms.iter_mut() {
                        *c += 1; // exposed mutable counts must be honoured
                    }
                    for c in counts.iter_mut().filter(|c| **c > 0) {
                        *c += 1;
                    }
                }
            }
            check(&ms, &counts, &ctx);
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn from_iterator_and_into_iter_round_trip() {
    let mut st = Stress::new("multiset/from_iter");
    for case in 0..200 * scale() {
        let mut rng = Rng::new(9000 + case as u64);
        let items: Vec<i32> = (0..rng.below(60)).map(|_| rng.below(V) as i32).collect();
        let ms: MultiSet<i32> = items.iter().copied().collect();
        let mut counts = [0usize; V];
        for &x in &items {
            counts[x as usize] += 1;
        }
        check(&ms, &counts, &format!("from_iter case {case}"));
        let drained: Vec<(i32, usize)> = ms.into_iter().collect();
        let expect: Vec<(i32, usize)> = (0..V).filter(|&v| counts[v] > 0).map(|v| (v as i32, counts[v])).collect();
        assert_eq!(drained, expect, "into_iter case {case}");
        st.case();
        st.ops(items.len() + V);
    }
    st.done();
}
