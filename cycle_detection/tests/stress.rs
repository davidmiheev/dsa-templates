//! Stress test: cycle detection against reachability / union-find brute force, and validation
//! of every returned cycle.

use cycle_detection::{cycle_vertices_directed, cycle_vertices_undirected, find_cycle_directed, find_cycle_undirected};
use std::collections::HashSet;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn directed_has_cycle(n: usize, edges: &[(usize, usize)]) -> bool {
    let mut reach = vec![vec![false; n]; n];
    for &(u, v) in edges {
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
    (0..n).any(|v| reach[v][v])
}

fn undirected_has_cycle(n: usize, edges: &[(usize, usize)]) -> bool {
    let mut p: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        if p[x] != x {
            let r = find(p, p[x]);
            p[x] = r;
        }
        p[x]
    }
    for &(u, v) in edges {
        let (a, b) = (find(&mut p, u), find(&mut p, v));
        if a == b {
            return true;
        }
        p[a] = b;
    }
    false
}

#[test]
fn directed_matches_reachability() {
    let mut st = Stress::new("cycle_detection/directed");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 10) as usize;
        let edges: Vec<(usize, usize)> = (0..rng.range(0, 14)).map(|_| (rng.below(n), rng.below(n))).collect();
        let cycle = find_cycle_directed(n, &edges);
        assert_eq!(!cycle.is_empty(), directed_has_cycle(n, &edges), "case {case}");
        if !cycle.is_empty() {
            let distinct: HashSet<_> = cycle.iter().collect();
            assert_eq!(distinct.len(), cycle.len(), "edges distinct");
            for i in 0..cycle.len() {
                assert_eq!(edges[cycle[i]].1, edges[cycle[(i + 1) % cycle.len()]].0, "consecutive, case {case}");
            }
            let vs = cycle_vertices_directed(&edges, &cycle);
            assert_eq!(vs.iter().collect::<HashSet<_>>().len(), vs.len(), "simple cycle, case {case}");
        }
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn undirected_matches_union_find() {
    let mut st = Stress::new("cycle_detection/undirected");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64 + 50);
        let n = rng.range(1, 10) as usize;
        let edges: Vec<(usize, usize)> = (0..rng.range(0, 10)).map(|_| (rng.below(n), rng.below(n))).collect();
        let cycle = find_cycle_undirected(n, &edges);
        assert_eq!(!cycle.is_empty(), undirected_has_cycle(n, &edges), "case {case}: {edges:?}");
        if !cycle.is_empty() {
            assert_eq!(cycle.iter().collect::<HashSet<_>>().len(), cycle.len());
            let vs = cycle_vertices_undirected(&edges, &cycle);
            assert_eq!(vs.len(), cycle.len());
            assert_eq!(vs.iter().collect::<HashSet<_>>().len(), vs.len(), "simple, case {case}");
            for i in 0..cycle.len() {
                let (a, b) = edges[cycle[i]];
                let (x, y) = (vs[i], vs[(i + 1) % vs.len()]);
                assert!((a, b) == (x, y) || (a, b) == (y, x), "edge joins consecutive vertices, case {case}");
            }
        }
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn large_graphs() {
    let mut st = Stress::new("cycle_detection/large");
    let n = 1_000_000;
    st.case();
    let mut path: Vec<(usize, usize)> = (1..n).map(|i| (i - 1, i)).collect();
    assert!(find_cycle_directed(n, &path).is_empty());
    assert!(find_cycle_undirected(n, &path).is_empty());
    path.push((n - 1, 0));
    assert_eq!(find_cycle_directed(n, &path).len(), n);
    assert_eq!(find_cycle_undirected(n, &path).len(), n);
    st.ops(n);
    st.done();
}
