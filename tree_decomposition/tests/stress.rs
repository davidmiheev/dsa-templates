//! Stress test: validity of every decomposition, and agreement with a brute-force K4-minor-free
//! test (treewidth <= 2) on small graphs.

use tree_decomposition::tree_decomposition_width2;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

/// Treewidth <= 2 iff repeated removal of degree <= 2 vertices (with fill-in) empties the graph;
/// verified independently by trying every elimination order on tiny graphs.
fn brute_tw_le_2(n: usize, adj: &[Vec<bool>]) -> bool {
    fn rec(alive: u32, adj: &mut Vec<Vec<bool>>, n: usize) -> bool {
        if alive == 0 {
            return true;
        }
        for v in 0..n {
            if alive >> v & 1 == 0 {
                continue;
            }
            let nb: Vec<usize> = (0..n).filter(|&w| w != v && alive >> w & 1 == 1 && adj[v][w]).collect();
            if nb.len() > 2 {
                continue;
            }
            let saved = adj.clone();
            for &a in &nb {
                for &b in &nb {
                    if a != b {
                        adj[a][b] = true;
                    }
                }
            }
            let ok = rec(alive & !(1 << v), adj, n);
            *adj = saved;
            if ok {
                return true;
            }
        }
        false
    }
    let mut a = adj.to_vec();
    rec((1u32 << n) - 1, &mut a, n)
}

#[test]
fn decompositions_are_valid() {
    let mut st = Stress::new("tree_decomposition");
    let (mut yes, mut no) = (0, 0);
    for case in 0..4000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 8) as usize;
        let density = rng.range(1, 6) as u64;
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
        let d = tree_decomposition_width2(n, &edges);
        let expect = brute_tw_le_2(n, &adj);
        assert_eq!(!d.bags.is_empty(), expect, "case {case}: {edges:?}");
        if expect {
            yes += 1;
            let k = d.bags.len();
            assert!(k <= 10 * n + edges.len() + 10);
            assert!(d.bags.iter().all(|b| b.len() <= 3));
            assert_eq!(d.edges.len(), k - 1);
            // the bag graph is a tree
            let mut label: Vec<usize> = (0..k).collect();
            for &(a, b) in &d.edges {
                let (x, y) = (label[a], label[b]);
                assert_ne!(x, y, "case {case}: cycle in bag tree");
                for l in label.iter_mut() {
                    if *l == y {
                        *l = x;
                    }
                }
            }
            // every edge is inside some bag
            for &(u, v) in &edges {
                assert!(d.bags.iter().any(|b| b.contains(&u) && b.contains(&v)), "case {case}: edge {u}-{v} uncovered");
            }
            // bags containing a vertex are connected in the tree
            for v in 0..n {
                let with: Vec<usize> = (0..k).filter(|&i| d.bags[i].contains(&v)).collect();
                assert!(!with.is_empty(), "case {case}: vertex {v} in no bag");
                let mut comp: Vec<usize> = (0..k).collect();
                for &(a, b) in &d.edges {
                    if d.bags[a].contains(&v) && d.bags[b].contains(&v) {
                        let (x, y) = (comp[a], comp[b]);
                        for c in comp.iter_mut() {
                            if *c == y {
                                *c = x;
                            }
                        }
                    }
                }
                assert!(with.iter().all(|&i| comp[i] == comp[with[0]]), "case {case}: vertex {v} bags disconnected");
            }
        } else {
            no += 1;
        }
        st.case();
    }
    assert!(yes > 200 && no > 200, "{yes}/{no}");
    st.done();
}
