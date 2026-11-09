//! Stress test: Dreyfus–Wagner against enumeration of all edge subsets.

use steiner_tree::steiner_tree;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_edge_subset_enumeration() {
    let mut st = Stress::new("steiner_tree");
    for case in 0..2500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 7) as usize;
        // a random spanning tree plus extra edges, so the graph is connected
        let mut edges: Vec<(usize, usize, u64)> = (1..n).map(|v| (rng.below(v), v, rng.range(0, 9) as u64)).collect();
        for _ in 0..rng.range(0, 6) {
            edges.push((rng.below(n), rng.below(n), rng.range(0, 9) as u64));
        }
        let k = rng.range(1, n as i64) as usize;
        let mut terminals: Vec<usize> = (0..n).collect();
        for i in 0..k {
            let j = i + rng.below(n - i);
            terminals.swap(i, j);
        }
        terminals.truncate(k);
        terminals.sort();
        let m = edges.len();
        let mut best = u64::MAX;
        for mask in 0u32..1 << m {
            let mut label: Vec<usize> = (0..n).collect();
            let mut cost = 0;
            for (i, &(u, v, w)) in edges.iter().enumerate() {
                if mask >> i & 1 == 1 {
                    cost += w;
                    let (a, b) = (label[u], label[v]);
                    for l in label.iter_mut() {
                        if *l == b {
                            *l = a;
                        }
                    }
                }
            }
            if terminals.iter().all(|&t| label[t] == label[terminals[0]]) {
                best = best.min(cost);
            }
        }
        let (cost, chosen) = steiner_tree(n, &edges, &terminals);
        assert_eq!(cost, best, "case {case}: {edges:?} {terminals:?}");
        // the chosen edges connect the terminals and weigh no more than the optimum
        let mut label: Vec<usize> = (0..n).collect();
        let mut w = 0;
        for &e in &chosen {
            w += edges[e].2;
            let (a, b) = (label[edges[e].0], label[edges[e].1]);
            for l in label.iter_mut() {
                if *l == b {
                    *l = a;
                }
            }
        }
        assert!(terminals.iter().all(|&t| label[t] == label[terminals[0]]), "case {case}: not connected");
        assert_eq!(w, best, "case {case}: chosen edges weigh {w}, optimum {best}");
        st.case();
    }
    st.done();
}
