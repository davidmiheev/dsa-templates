//! Stress test: maximum independent set / clique against enumeration of all vertex subsets.

use independent_set::{maximum_clique, maximum_independent_set};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_subset_enumeration() {
    let mut st = Stress::new("independent_set");
    for case in 0..2500 * scale() {
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
        let (mut best_i, mut best_c) = (0, 0);
        for mask in 0u32..1 << n {
            let vs: Vec<usize> = (0..n).filter(|&v| mask >> v & 1 == 1).collect();
            let pairs = || vs.iter().flat_map(|&a| vs.iter().map(move |&b| (a, b))).filter(|&(a, b)| a < b);
            if pairs().all(|(a, b)| !adj[a][b]) {
                best_i = best_i.max(vs.len());
            }
            if pairs().all(|(a, b)| adj[a][b]) {
                best_c = best_c.max(vs.len());
            }
        }
        let is = maximum_independent_set(n, &edges);
        assert_eq!(is.len(), best_i, "case {case}: {edges:?}");
        assert!(is.iter().all(|&a| is.iter().all(|&b| a == b || !adj[a][b])));
        let cl = maximum_clique(n, &edges);
        assert_eq!(cl.len(), best_c, "case {case}: {edges:?}");
        assert!(cl.iter().all(|&a| cl.iter().all(|&b| a == b || adj[a][b])));
        st.case();
    }
    st.done();
}
