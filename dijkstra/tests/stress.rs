//! Stress test: Dijkstra against Bellman-Ford on random directed / undirected multigraphs,
//! with zero weights, self-loops and parallel edges; paths are checked edge by edge.

use csr_graph::Csr;
use dijkstra::dijkstra;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const INF: u64 = u64::MAX / 4;

fn bellman_ford(n: usize, arcs: &[(usize, usize, u64)], s: usize) -> Vec<u64> {
    let mut d = vec![INF; n];
    d[s] = 0;
    for _ in 0..n {
        for &(u, v, w) in arcs {
            if d[u] < INF && d[u] + w < d[v] {
                d[v] = d[u] + w;
            }
        }
    }
    d
}

#[test]
fn matches_bellman_ford() {
    let mut st = Stress::new("dijkstra/naive");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 12) as usize;
        let m = rng.range(0, 30) as usize;
        let directed = case % 2 == 0;
        let wmax = *[1usize, 5, 1000].get(case % 3).unwrap();
        let edges: Vec<(usize, usize)> = (0..m).map(|_| (rng.below(n), rng.below(n))).collect();
        let w: Vec<u64> = (0..m).map(|_| rng.below(wmax + 1) as u64).collect();
        let g = if directed { Csr::directed(n, &edges) } else { Csr::undirected(n, &edges) };
        let mut arcs: Vec<(usize, usize, u64)> = edges.iter().zip(&w).map(|(&(u, v), &w)| (u, v, w)).collect();
        if !directed {
            arcs.extend(edges.iter().zip(&w).map(|(&(u, v), &w)| (v, u, w)));
        }
        let s = rng.below(n);
        let sp = dijkstra(&g, &w, s);
        let want = bellman_ford(n, &arcs, s);
        for v in 0..n {
            assert_eq!(sp.reachable(v), want[v] < INF, "reachable({v}) case {case}");
            if want[v] < INF {
                assert_eq!(sp.dist[v], want[v], "dist[{v}] case {case}");
                // The returned path is a real path of that length.
                let vs = sp.path_vertices(v);
                let es = sp.path_edges(v);
                assert_eq!((vs[0], *vs.last().unwrap(), vs.len()), (s, v, es.len() + 1));
                let mut total = 0;
                for (i, &e) in es.iter().enumerate() {
                    let (a, b) = edges[e];
                    assert!((a == vs[i] && b == vs[i + 1]) || (!directed && b == vs[i] && a == vs[i + 1]), "edge {e} joins path vertices, case {case}");
                    total += w[e];
                }
                assert_eq!(total, want[v], "path length case {case}");
            }
        }
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn large_graph() {
    let mut st = Stress::new("dijkstra/large");
    let (n, m) = (500_000, 1_000_000);
    let mut rng = Rng::new(3);
    let edges: Vec<(usize, usize)> = (0..m).map(|_| (rng.below(n), rng.below(n))).collect();
    let w: Vec<u64> = (0..m).map(|_| rng.below(1_000_000_000) as u64).collect();
    st.case();
    let sp = dijkstra(&Csr::directed(n, &edges), &w, 0);
    assert!(sp.reachable(0) && sp.dist[0] == 0);
    // Triangle inequality over every edge from a reachable tail.
    for (e, &(u, v)) in edges.iter().enumerate() {
        if sp.reachable(u) {
            assert!(sp.reachable(v) && sp.dist[v] <= sp.dist[u] + w[e]);
        }
    }
    st.ops(m);
    st.done();
}
