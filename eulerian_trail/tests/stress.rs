//! Stress test: Eulerian trails against exhaustive search over edge orders (tiny graphs) and
//! validation of every returned trail (larger graphs).

use eulerian_trail::{eulerian_trail_directed, eulerian_trail_undirected, Trail};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn valid(edges: &[(usize, usize)], t: &Trail, directed: bool) -> bool {
    if t.vertices.len() != t.edges.len() + 1 || t.edges.len() != edges.len() {
        return false;
    }
    let mut seen = vec![false; edges.len()];
    for (i, &e) in t.edges.iter().enumerate() {
        if e >= edges.len() || seen[e] {
            return false;
        }
        seen[e] = true;
        let (a, b) = edges[e];
        let (x, y) = (t.vertices[i], t.vertices[i + 1]);
        if !((a == x && b == y) || (!directed && a == y && b == x)) {
            return false;
        }
    }
    true
}

/// Is there any order of all edges forming a trail?
fn exists(n: usize, edges: &[(usize, usize)], directed: bool) -> bool {
    fn rec(edges: &[(usize, usize)], used: &mut Vec<bool>, at: usize, left: usize, directed: bool) -> bool {
        if left == 0 {
            return true;
        }
        for (e, &(a, b)) in edges.iter().enumerate() {
            if used[e] {
                continue;
            }
            for (x, y) in if directed { vec![(a, b)] } else { vec![(a, b), (b, a)] } {
                if x == at {
                    used[e] = true;
                    if rec(edges, used, y, left - 1, directed) {
                        used[e] = false;
                        return true;
                    }
                    used[e] = false;
                }
            }
        }
        false
    }
    if edges.is_empty() {
        return true;
    }
    (0..n).any(|s| rec(edges, &mut vec![false; edges.len()], s, edges.len(), directed))
}

#[test]
fn matches_exhaustive_search() {
    let mut st = Stress::new("eulerian_trail/naive");
    for case in 0..2000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 6) as usize;
        let edges: Vec<(usize, usize)> = (0..rng.range(0, 8)).map(|_| (rng.below(n), rng.below(n))).collect();
        for directed in [true, false] {
            let t = if directed { eulerian_trail_directed(n, &edges) } else { eulerian_trail_undirected(n, &edges) };
            assert_eq!(t.exists(), exists(n, &edges, directed), "directed={directed} case {case}: {edges:?}");
            if t.exists() {
                assert!(valid(&edges, &t, directed), "invalid trail {t:?} for {edges:?} (directed={directed})");
            }
        }
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn large_circuits() {
    let mut st = Stress::new("eulerian_trail/large");
    let n = 200_000;
    let mut rng = Rng::new(5);
    // A union of random cycles is Eulerian when connected: chain them through vertex 0.
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for _ in 0..1000 {
        let len = 2 + rng.below(400);
        let mut cur = 0;
        for _ in 1..len {
            let nxt = rng.below(n);
            edges.push((cur, nxt));
            cur = nxt;
        }
        edges.push((cur, 0));
    }
    st.case();
    let t = eulerian_trail_directed(n, &edges);
    assert!(t.exists() && valid(&edges, &t, true));
    assert_eq!(t.vertices.first(), t.vertices.last());
    let t = eulerian_trail_undirected(n, &edges);
    assert!(t.exists() && valid(&edges, &t, false));
    st.ops(2 * edges.len());
    st.done();
}

#[test]
fn edgeless_graph_has_a_trivial_trail() {
    assert_eq!(eulerian_trail_directed(3, &[]).vertices, vec![0]);
    assert_eq!(eulerian_trail_undirected(3, &[]).vertices, vec![0]);
}
