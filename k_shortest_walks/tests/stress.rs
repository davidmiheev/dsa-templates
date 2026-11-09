//! Stress test: best-first enumeration of all walks as the reference.

use k_shortest_walks::k_shortest_walks;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

/// Lengths of the k shortest walks by expanding walks in order of length. Requires that no
/// cycle of weight 0 exists (otherwise the search may not terminate).
fn brute(n: usize, edges: &[(usize, usize, u64)], s: usize, t: usize, k: usize) -> Vec<u64> {
    // only walk through vertices that can still reach t
    let mut can = vec![false; n];
    can[t] = true;
    for _ in 0..n {
        for &(a, b, _) in edges {
            if can[b] {
                can[a] = true;
            }
        }
    }
    let mut pq = BinaryHeap::new();
    if can[s] {
        pq.push(Reverse((0u64, s)));
    }
    let mut out = Vec::new();
    let mut pops = 0;
    while let Some(Reverse((c, v))) = pq.pop() {
        pops += 1;
        if pops > 200_000 {
            break;
        }
        if v == t {
            out.push(c);
            if out.len() == k {
                break;
            }
        }
        for &(a, b, w) in edges {
            if a == v && can[b] {
                pq.push(Reverse((c + w, b)));
            }
        }
    }
    out
}

#[test]
fn matches_best_first_enumeration() {
    let mut st = Stress::new("k_shortest_walks");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 6) as usize;
        let m = rng.range(0, 10) as usize;
        let edges: Vec<(usize, usize, u64)> = (0..m)
            .map(|_| {
                let (u, v) = (rng.below(n), rng.below(n));
                // zero weights only on forward edges, so no zero-weight cycle exists
                let w = if u < v { rng.range(0, 4) } else { rng.range(1, 4) };
                (u, v, w as u64)
            })
            .collect();
        let (s, t) = (rng.below(n), rng.below(n));
        let k = rng.range(1, 12) as usize;
        let got = k_shortest_walks(n, &edges, s, t, k);
        let expect = brute(n, &edges, s, t, k);
        assert_eq!(got, expect, "case {case}: n={n} {edges:?} s={s} t={t} k={k}");
        st.case();
        st.ops(got.len());
    }
    st.done();
}
