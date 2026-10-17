//! Stress test: BFS against Bellman-Ford-style relaxation (distances) and a
//! layer-by-layer reference (visit order), on directed and undirected graphs.

use bfs::{bfs_distances, bfs_order};
use std::collections::HashMap;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn random_graph(rng: &mut Rng, n: usize, m: usize, undirected: bool) -> Vec<Vec<usize>> {
    let mut g = vec![Vec::new(); n];
    for _ in 0..m {
        let (u, v) = (rng.below(n), rng.below(n)); // self-loops and parallel edges included
        g[u].push(v);
        if undirected {
            g[v].push(u);
        }
    }
    g
}

fn ref_dist(g: &[Vec<usize>], s: usize) -> Vec<Option<usize>> {
    let mut d = vec![None; g.len()];
    d[s] = Some(0);
    for _ in 0..g.len() {
        for u in 0..g.len() {
            if let Some(du) = d[u] {
                for &v in &g[u] {
                    if d[v].map_or(true, |dv| du + 1 < dv) {
                        d[v] = Some(du + 1);
                    }
                }
            }
        }
    }
    d
}

/// Same visit order as a queue BFS, computed frontier by frontier.
fn ref_order(g: &[Vec<usize>], s: usize) -> Vec<usize> {
    let mut seen = vec![false; g.len()];
    seen[s] = true;
    let (mut order, mut frontier) = (vec![s], vec![s]);
    while !frontier.is_empty() {
        let mut next = vec![];
        for &u in &frontier {
            for &v in &g[u] {
                if !seen[v] {
                    seen[v] = true;
                    next.push(v);
                }
            }
        }
        order.extend(&next);
        frontier = next;
    }
    order
}

#[test]
fn matches_references() {
    let mut st = Stress::new("bfs/naive");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 25) as usize;
        let m = rng.below(3 * n + 1);
        let g = random_graph(&mut rng, n, m, case % 2 == 0);
        st.case();
        for s in 0..n {
            let want = ref_dist(&g, s);
            let got = bfs_distances(&g, s);
            let want_map: HashMap<usize, usize> = want.iter().enumerate().filter_map(|(v, d)| d.map(|d| (v, d))).collect();
            assert_eq!(got, want_map, "distances from {s}, case {case}");

            let order = bfs_order(&g, s);
            assert_eq!(order, ref_order(&g, s), "order from {s}, case {case}");
            // Independent validity checks on the order itself.
            assert_eq!(order[0], s);
            assert_eq!(order.len(), want_map.len(), "exactly the reachable vertices, case {case}");
            assert!(order.windows(2).all(|w| want[w[0]] <= want[w[1]]), "non-decreasing distance, case {case}");
            st.ops(3);
        }
    }
    st.done();
}

#[test]
fn source_out_of_range_panics() {
    let g = vec![vec![1], vec![0]];
    assert!(std::panic::catch_unwind(|| bfs_distances(&g, 2)).is_err());
    assert!(std::panic::catch_unwind(|| bfs_order(&g, 2)).is_err());
}

#[test]
fn large_sparse_graph_distances_are_consistent() {
    // n = 200k: distances must satisfy dist[v] <= dist[u] + 1 on every edge, and every
    // reached vertex other than the source needs a predecessor exactly one step closer.
    let mut st = Stress::new("bfs/large");
    let n = 200_000;
    let mut rng = Rng::new(8);
    let g = random_graph(&mut rng, n, 300_000, true);
    st.case();
    let d = bfs_distances(&g, 0);
    let mut has_pred = vec![false; n];
    has_pred[0] = true;
    for u in 0..n {
        if let Some(&du) = d.get(&u) {
            for &v in &g[u] {
                let dv = *d.get(&v).expect("neighbour of a reached vertex is reached");
                assert!(dv <= du + 1);
                has_pred[v] |= dv == du + 1;
            }
        }
    }
    assert!(d.keys().all(|&v| has_pred[v]));
    assert_eq!(bfs_order(&g, 0).len(), d.len());
    st.ops(2 * n + 600_000);
    st.done();
}
