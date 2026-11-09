//! Stress test: the dropped edge against recomputing the MSF with Kruskal after every insertion.

use incremental_msf::{IncrementalMsf, NONE};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn kruskal(n: usize, edges: &[(usize, usize, u32)]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..edges.len()).collect();
    order.sort_by_key(|&i| edges[i].2);
    let mut comp: Vec<usize> = (0..n).collect();
    let mut chosen = Vec::new();
    for i in order {
        let (u, v, _) = edges[i];
        let (a, b) = (comp[u], comp[v]);
        if a != b {
            for c in comp.iter_mut() {
                if *c == b {
                    *c = a;
                }
            }
            chosen.push(i);
        }
    }
    chosen.sort();
    chosen
}

#[test]
fn matches_kruskal_recomputation() {
    let mut st = Stress::new("incremental_msf");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 9) as usize;
        let m = rng.range(1, 25) as usize;
        // distinct weights
        let mut weights: Vec<u32> = (1..=m as u32).collect();
        for i in (1..m).rev() {
            weights.swap(i, rng.below(i + 1));
        }
        let mut msf = IncrementalMsf::new(n, m);
        let mut edges: Vec<(usize, usize, u32)> = Vec::new();
        let mut prev: Vec<usize> = Vec::new();
        for i in 0..m {
            let (u, v) = (rng.below(n), rng.below(n));
            edges.push((u, v, weights[i]));
            let dropped = msf.add_edge(u, v, weights[i]);
            let now: Vec<usize> = kruskal(n, &edges)
                .into_iter()
                .filter(|&e| edges[e].0 != edges[e].1)
                .collect();
            let mut with_new = prev.clone();
            with_new.push(i);
            with_new.sort();
            if dropped == NONE {
                assert_eq!(now, with_new, "case {case} step {i}");
            } else {
                with_new.retain(|&e| e != dropped);
                assert_eq!(now, with_new, "case {case} step {i}: dropped {dropped}");
            }
            prev = now;
            st.ops(1);
        }
        st.case();
    }
    st.done();
}
