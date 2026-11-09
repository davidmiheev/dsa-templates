//! Stress test: components of the complement against a plain DFS on the explicit complement.

use complement_graph::complement_components;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_explicit_complement() {
    let mut st = Stress::new("complement_graph");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 12) as usize;
        let density = rng.range(0, 10) as u64;
        let mut adj = vec![vec![false; n]; n];
        let mut edges = Vec::new();
        for u in 0..n {
            for v in u + 1..n {
                if rng.chance(density, 10) {
                    adj[u][v] = true;
                    adj[v][u] = true;
                    edges.push((u, v));
                }
            }
        }
        let mut label = vec![usize::MAX; n];
        let mut k = 0;
        for s in 0..n {
            if label[s] != usize::MAX {
                continue;
            }
            label[s] = k;
            let mut stack = vec![s];
            while let Some(v) = stack.pop() {
                for w in 0..n {
                    if w != v && !adj[v][w] && label[w] == usize::MAX {
                        label[w] = k;
                        stack.push(w);
                    }
                }
            }
            k += 1;
        }
        let comps = complement_components(n, &edges);
        assert_eq!(comps.len(), k, "case {case}");
        let mut seen = vec![false; n];
        for c in &comps {
            for &v in c {
                assert!(!seen[v]);
                seen[v] = true;
                assert_eq!(label[v], label[c[0]], "case {case}: {edges:?}");
            }
        }
        assert!(seen.iter().all(|&b| b));
        st.case();
        st.ops(n);
    }
    st.done();
}
