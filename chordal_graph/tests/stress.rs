//! Stress test: recognition against brute force, and validation of both certificates.

use chordal_graph::{recognize, Chordality};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn is_peo(n: usize, adj: &[Vec<bool>], peo: &[usize]) -> bool {
    let mut pos = vec![usize::MAX; n];
    for (i, &v) in peo.iter().enumerate() {
        assert_eq!(pos[v], usize::MAX, "repeated vertex");
        pos[v] = i;
    }
    for &v in peo {
        let later: Vec<usize> = (0..n).filter(|&w| adj[v][w] && pos[w] > pos[v]).collect();
        for &a in &later {
            for &b in &later {
                if a != b && !adj[a][b] {
                    return false;
                }
            }
        }
    }
    true
}

/// Chordal iff some elimination (repeatedly remove a simplicial vertex) clears the graph.
fn brute_chordal(n: usize, adj: &[Vec<bool>]) -> bool {
    let mut alive = vec![true; n];
    for _ in 0..n {
        let mut found = false;
        for v in 0..n {
            if !alive[v] {
                continue;
            }
            let nb: Vec<usize> = (0..n).filter(|&w| alive[w] && adj[v][w]).collect();
            if nb.iter().all(|&a| nb.iter().all(|&b| a == b || adj[a][b])) {
                alive[v] = false;
                found = true;
                break;
            }
        }
        if !found {
            return false;
        }
    }
    true
}

#[test]
fn recognition_and_certificates() {
    let mut st = Stress::new("chordal_graph");
    let (mut yes, mut no) = (0, 0);
    for case in 0..6000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 10) as usize;
        let density = rng.range(1, 9) as u64;
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
        match recognize(n, &edges) {
            Chordality::Chordal { peo } => {
                yes += 1;
                assert_eq!(peo.len(), n);
                assert!(is_peo(n, &adj, &peo), "case {case}: {edges:?} bad PEO {peo:?}");
                assert!(brute_chordal(n, &adj));
            }
            Chordality::Hole { cycle } => {
                no += 1;
                assert!(!brute_chordal(n, &adj), "case {case}: {edges:?}");
                let k = cycle.len();
                assert!(k >= 4, "case {case}: short cycle {cycle:?}");
                for i in 0..k {
                    for j in i + 1..k {
                        assert_ne!(cycle[i], cycle[j]);
                        let consecutive = j == i + 1 || (i == 0 && j == k - 1);
                        assert_eq!(adj[cycle[i]][cycle[j]], consecutive, "case {case}: {edges:?} cycle {cycle:?}");
                    }
                }
            }
        }
        st.case();
        st.ops(n);
    }
    assert!(yes > 100 && no > 100, "generator should cover both outcomes ({yes}/{no})");
    st.done();
}
