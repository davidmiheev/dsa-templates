//! Stress test: `Scc` against mutual reachability (Floyd-Warshall) and the topological order
//! contract, on small random graphs and a few large ones.

use scc::Scc;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_mutual_reachability() {
    let mut st = Stress::new("scc/naive");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 14) as usize;
        let m = rng.range(0, 3 * n as i64) as usize;
        let edges: Vec<(usize, usize)> = (0..m).map(|_| (rng.below(n), rng.below(n))).collect();
        let mut reach = vec![vec![false; n]; n];
        for i in 0..n {
            reach[i][i] = true;
        }
        for &(u, v) in &edges {
            reach[u][v] = true;
        }
        for k in 0..n {
            for i in 0..n {
                for j in 0..n {
                    if reach[i][k] && reach[k][j] {
                        reach[i][j] = true;
                    }
                }
            }
        }
        let s = Scc::from_edges(n, &edges);
        for u in 0..n {
            for v in 0..n {
                assert_eq!(s.comp[u] == s.comp[v], reach[u][v] && reach[v][u], "case {case} ({u}, {v})");
            }
        }
        for &(u, v) in &edges {
            assert!(s.comp[u] <= s.comp[v], "topological order, case {case}");
        }
        assert_eq!(s.count, s.comp.iter().max().unwrap() + 1);
        assert_eq!(s.groups().iter().map(Vec::len).sum::<usize>(), n);
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn large_graphs() {
    let mut st = Stress::new("scc/large");
    let n = 1_000_000;
    st.case();
    // One huge cycle, a long path, and a random graph.
    let cycle: Vec<(usize, usize)> = (0..n).map(|i| (i, (i + 1) % n)).collect();
    assert_eq!(Scc::from_edges(n, &cycle).count, 1);
    let path: Vec<(usize, usize)> = (1..n).map(|i| (i - 1, i)).collect();
    let s = Scc::from_edges(n, &path);
    assert_eq!(s.count, n);
    assert!(s.comp.windows(2).all(|w| w[0] < w[1]));
    let mut rng = Rng::new(9);
    let random: Vec<(usize, usize)> = (0..2 * n).map(|_| (rng.below(n), rng.below(n))).collect();
    let s = Scc::from_edges(n, &random);
    assert!(random.iter().all(|&(u, v)| s.comp[u] <= s.comp[v]));
    st.ops(n);
    st.done();
}
