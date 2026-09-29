//! Stress test: Dijkstra against Bellman-Ford on random directed graphs (zero weights,
//! parallel edges and self-loops included), plus an optimality certificate at scale.

use dijkstra::dijkstra;
use std::collections::HashMap;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

type Graph = Vec<Vec<(usize, u64)>>;

fn random_graph(rng: &mut Rng, n: usize, m: usize, max_w: u64, min_w: u64) -> Graph {
    let mut g = vec![Vec::new(); n];
    for _ in 0..m {
        let (u, v) = (rng.below(n), rng.below(n));
        g[u].push((v, min_w + rng.next_u64() % (max_w - min_w + 1)));
    }
    g
}

fn bellman_ford(g: &Graph, s: usize) -> HashMap<usize, u64> {
    let mut d: Vec<Option<u64>> = vec![None; g.len()];
    d[s] = Some(0);
    for _ in 0..g.len() {
        let mut changed = false;
        for u in 0..g.len() {
            if let Some(du) = d[u] {
                for &(v, w) in &g[u] {
                    if d[v].map_or(true, |dv| du + w < dv) {
                        d[v] = Some(du + w);
                        changed = true;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
    d.iter().enumerate().filter_map(|(v, x)| x.map(|x| (v, x))).collect()
}

#[test]
fn matches_bellman_ford() {
    let mut st = Stress::new("dijkstra/naive");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 25) as usize;
        let max_w = *[0u64, 3, 100, 1_000_000_000_000].get(case % 4).unwrap();
        let m = rng.below(4 * n + 1);
        let g = random_graph(&mut rng, n, m, max_w, 0);
        st.case();
        for s in 0..n {
            assert_eq!(dijkstra(&g, s), bellman_ford(&g, s), "from {s}, case {case}, max_w {max_w}");
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn source_out_of_range_panics() {
    let g: Graph = vec![vec![(0, 1)]];
    assert!(std::panic::catch_unwind(|| dijkstra(&g, 1)).is_err());
}

#[test]
fn large_graph_satisfies_optimality_certificate() {
    // With strictly positive weights, distances are optimal iff (1) no edge can be relaxed and
    // (2) every reached vertex other than the source has a "tight" incoming edge.
    let mut st = Stress::new("dijkstra/large");
    let n = 200_000;
    let mut rng = Rng::new(6);
    let g = random_graph(&mut rng, n, 600_000, 1_000_000, 1);
    st.case();
    let d = dijkstra(&g, 0);
    let mut tight = vec![false; n];
    tight[0] = true;
    assert_eq!(d[&0], 0);
    for u in 0..n {
        if let Some(&du) = d.get(&u) {
            for &(v, w) in &g[u] {
                let dv = *d.get(&v).expect("successor of a reached vertex is reached");
                assert!(dv <= du + w, "edge {u}->{v} can still be relaxed");
                tight[v] |= dv == du + w;
            }
        }
    }
    assert!(d.keys().all(|&v| tight[v]), "every reached vertex needs a tight predecessor");
    st.ops(n + 600_000);
    st.done();
}
