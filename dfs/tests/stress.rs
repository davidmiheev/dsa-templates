//! Stress test: `dfs_order` against a recursive DFS, and `has_cycle_undirected`
//! against the edge-count criterion (a simple graph has a cycle iff m > n - components).

use dfs::{dfs_order, has_cycle_undirected};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn recursive_preorder(g: &[Vec<usize>], v: usize, seen: &mut Vec<bool>, out: &mut Vec<usize>) {
    seen[v] = true;
    out.push(v);
    for &u in &g[v] {
        if !seen[u] {
            recursive_preorder(g, u, seen, out);
        }
    }
}

#[test]
fn dfs_order_matches_recursive_dfs() {
    let mut st = Stress::new("dfs/order");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 30) as usize;
        let mut g = vec![Vec::new(); n];
        for _ in 0..rng.below(3 * n + 1) {
            let (u, v) = (rng.below(n), rng.below(n)); // directed, with loops and parallel edges
            g[u].push(v);
            if case % 2 == 0 {
                g[v].push(u);
            }
        }
        st.case();
        for s in 0..n {
            let mut out = vec![];
            recursive_preorder(&g, s, &mut vec![false; n], &mut out);
            assert_eq!(dfs_order(&g, s), out, "from {s}, case {case}");
            st.ops(1);
        }
    }
    st.done();
}

fn components(n: usize, edges: &[(usize, usize)]) -> usize {
    let mut p: Vec<usize> = (0..n).collect();
    fn find(p: &mut Vec<usize>, x: usize) -> usize {
        if p[x] != x {
            let r = find(p, p[x]);
            p[x] = r;
        }
        p[x]
    }
    let mut c = n;
    for &(u, v) in edges {
        let (a, b) = (find(&mut p, u), find(&mut p, v));
        if a != b {
            p[a] = b;
            c -= 1;
        }
    }
    c
}

#[test]
fn cycle_detection_matches_edge_count_on_simple_graphs() {
    let mut st = Stress::new("dfs/cycle");
    let (mut with, mut without) = (0, 0);
    for case in 0..1000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 20) as usize;
        // Start from a random forest (acyclic), then add a few random extra edges.
        let mut edges: Vec<(usize, usize)> = vec![];
        for v in 1..n {
            if rng.chance(3, 4) {
                edges.push((rng.below(v), v));
            }
        }
        for _ in 0..rng.below(6) {
            let (u, v) = (rng.below(n), rng.below(n));
            let dup = edges.iter().any(|&(a, b)| (a, b) == (u, v) || (a, b) == (v, u));
            if u != v && !dup && rng.chance(2, 3) {
                edges.push((u, v));
            }
        }
        let mut g = vec![Vec::new(); n];
        for &(u, v) in &edges {
            g[u].push(v);
            g[v].push(u);
        }
        let expect = edges.len() > n - components(n, &edges);
        if expect { with += 1 } else { without += 1 }
        assert_eq!(has_cycle_undirected(&g), expect, "case {case} n={n} edges={edges:?}");
        st.case();
        st.ops(1);
    }
    assert!(with > 100 && without > 100, "generator must produce both outcomes ({with} cyclic, {without} acyclic)");
    st.done();
}

#[test]
fn self_loop_is_a_cycle() {
    assert!(has_cycle_undirected(&[vec![0]]));
    assert!(!has_cycle_undirected(&[vec![]]));
    assert!(!has_cycle_undirected(&[]));
}

#[test]
fn deep_path_does_not_overflow_the_stack() {
    // 1,000,000-vertex path (200,000 in debug builds): iterative implementations must cope.
    let n = if cfg!(debug_assertions) { 200_000 } else { 1_000_000 };
    let mut g = vec![Vec::new(); n];
    for v in 1..n {
        g[v - 1].push(v);
        g[v].push(v - 1);
    }
    let mut st = Stress::new("dfs/deep_path");
    st.case();
    assert_eq!(dfs_order(&g, 0), (0..n).collect::<Vec<_>>());
    assert!(!has_cycle_undirected(&g));
    g[n - 1].push(0);
    g[0].push(n - 1);
    assert!(has_cycle_undirected(&g));
    st.ops(3 * n);
    st.done();
}

#[test]
fn multigraph_cycles_match_edge_count() {
    // Parallel edges and self-loops are cycles too: a multigraph is acyclic iff m == n - components.
    let mut st = Stress::new("dfs/multigraph");
    let (mut with, mut without) = (0, 0);
    for case in 0..1000 * scale() {
        let mut rng = Rng::new(50_000 + case as u64);
        let n = rng.range(1, 15) as usize;
        let mut edges: Vec<(usize, usize)> = vec![];
        for v in 1..n {
            if rng.chance(3, 4) {
                edges.push((rng.below(v), v));
            }
        }
        for _ in 0..rng.below(3) {
            if rng.chance(1, 2) {
                edges.push(match rng.below(3) {
                    0 => (rng.below(n), rng.below(n)),                       // any pair (maybe a loop)
                    1 if !edges.is_empty() => edges[rng.below(edges.len())], // duplicate an existing edge
                    _ => {
                        let v = rng.below(n);
                        (v, v) // explicit self-loop
                    }
                });
            }
        }
        let mut g = vec![Vec::new(); n];
        for &(u, v) in &edges {
            g[u].push(v);
            if u != v {
                g[v].push(u); // a self-loop is listed once
            }
        }
        let expect = edges.len() > n - components(n, &edges);
        if expect { with += 1 } else { without += 1 }
        assert_eq!(has_cycle_undirected(&g), expect, "case {case} n={n} edges={edges:?}");
        st.case();
        st.ops(1);
    }
    assert!(with > 100 && without > 100, "generator must produce both outcomes ({with} cyclic, {without} acyclic)");
    st.done();
}
