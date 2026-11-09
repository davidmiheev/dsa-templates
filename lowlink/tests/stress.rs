//! Stress test: `Lowlink` against brute force (remove an edge / vertex and recount components;
//! blocks as classes of edges that share a simple cycle).

use csr_graph::Csr;
use lowlink::Lowlink;
use std::collections::BTreeSet;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn components(n: usize, edges: &[(usize, usize)], skip_edge: Option<usize>, skip_vertex: Option<usize>) -> Vec<usize> {
    let mut comp = vec![usize::MAX; n];
    let mut c = 0;
    for s in 0..n {
        if comp[s] != usize::MAX || Some(s) == skip_vertex {
            continue;
        }
        comp[s] = c;
        let mut stack = vec![s];
        while let Some(x) = stack.pop() {
            for (e, &(u, v)) in edges.iter().enumerate() {
                if Some(e) == skip_edge {
                    continue;
                }
                for (a, b) in [(u, v), (v, u)] {
                    if a == x && Some(b) != skip_vertex && comp[b] == usize::MAX {
                        comp[b] = c;
                        stack.push(b);
                    }
                }
            }
        }
        c += 1;
    }
    comp
}

fn count(comp: &[usize]) -> usize {
    comp.iter().filter(|&&c| c != usize::MAX).map(|&c| c).collect::<BTreeSet<_>>().len()
}

/// All simple cycles of a loop-free multigraph, as sets of edge indices.
fn simple_cycles(n: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut out = vec![];
    fn dfs(n: usize, edges: &[(usize, usize)], start: usize, v: usize, visited: &mut Vec<bool>, path: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        for (e, &(a, b)) in edges.iter().enumerate() {
            let w = if a == v { b } else if b == v { a } else { continue };
            if path.contains(&e) {
                continue;
            }
            if w == start {
                if path.len() >= 1 {
                    let mut cyc = path.clone();
                    cyc.push(e);
                    out.push(cyc);
                }
            } else if !visited[w] && w > start {
                visited[w] = true;
                path.push(e);
                dfs(n, edges, start, w, visited, path, out);
                path.pop();
                visited[w] = false;
            }
        }
    }
    for s in 0..n {
        let mut visited = vec![false; n];
        visited[s] = true;
        dfs(n, edges, s, s, &mut visited, &mut vec![], &mut out);
    }
    out
}

#[test]
fn matches_brute_force() {
    let mut st = Stress::new("lowlink/naive");
    for case in 0..1200 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 8) as usize;
        let loops = case % 3 == 0;
        let m = if n == 1 && !loops { 0 } else { rng.range(0, 10) as usize };
        let edges: Vec<(usize, usize)> = (0..m)
            .map(|_| loop {
                let (a, b) = (rng.below(n), rng.below(n));
                if loops || a != b {
                    break (a, b);
                }
            })
            .collect();
        let g = Csr::undirected(n, &edges);
        let ll = Lowlink::new(&g);
        let base = components(n, &edges, None, None);
        let base_count = count(&base);

        // Bridges.
        let want_bridges: BTreeSet<usize> = (0..m).filter(|&e| count(&components(n, &edges, Some(e), None)) > base_count).collect();
        assert_eq!(ll.bridges().into_iter().collect::<BTreeSet<_>>(), want_bridges, "bridges case {case}");
        assert_eq!(ll.bridges().len(), want_bridges.len());

        // Articulation points: removing v leaves more components than before, ignoring the
        // component that v alone formed.
        let want_art: Vec<usize> = (0..n)
            .filter(|&v| {
                let alone = !edges.iter().any(|&(a, b)| (a == v || b == v) && a != b);
                if alone {
                    return false;
                }
                let comp = components(n, &edges, None, Some(v));
                count(&comp) > base_count
            })
            .collect();
        assert_eq!(ll.articulation_points(), want_art, "articulation points case {case}");

        // 2-edge-connected components: components after removing all bridges.
        let kept: Vec<(usize, usize)> = (0..m).filter(|e| !want_bridges.contains(e)).map(|e| edges[e]).collect();
        let want = components(n, &kept, None, None);
        let (k, comp) = ll.two_edge_connected_components();
        assert_eq!(k, count(&want), "2ecc count case {case}");
        for u in 0..n {
            for v in 0..n {
                assert_eq!(comp[u] == comp[v], want[u] == want[v], "2ecc case {case} ({u}, {v})");
            }
        }

        // Blocks (loop-free graphs only): edges sharing a simple cycle are in one block.
        if !loops {
            let mut class: Vec<usize> = (0..m).collect();
            for cyc in simple_cycles(n, &edges) {
                let target = class[cyc[0]];
                let merged: Vec<usize> = cyc.iter().map(|&e| class[e]).collect();
                for c in class.iter_mut() {
                    if merged.contains(c) {
                        *c = target;
                    }
                }
            }
            let mut want_sets: BTreeSet<Vec<usize>> = BTreeSet::new();
            for id in class.iter().copied().collect::<BTreeSet<_>>() {
                let mut s: BTreeSet<usize> = BTreeSet::new();
                for e in 0..m {
                    if class[e] == id {
                        s.insert(edges[e].0);
                        s.insert(edges[e].1);
                    }
                }
                want_sets.insert(s.into_iter().collect());
            }
            for v in 0..n {
                if !edges.iter().any(|&(a, b)| a == v || b == v) {
                    want_sets.insert(vec![v]);
                }
            }
            let got: BTreeSet<Vec<usize>> = ll.biconnected_components().into_iter().collect();
            assert_eq!(got, want_sets, "blocks case {case}: {edges:?}");
        }
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn large_graphs_do_not_overflow_the_stack() {
    let mut st = Stress::new("lowlink/large");
    let n = 1_000_000;
    st.case();
    let path: Vec<(usize, usize)> = (1..n).map(|i| (i - 1, i)).collect();
    let ll = Lowlink::new(&Csr::undirected(n, &path));
    assert_eq!(ll.bridges().len(), n - 1);
    assert_eq!(ll.articulation_points().len(), n - 2);
    assert_eq!(ll.two_edge_connected_components().0, n);
    let cycle: Vec<(usize, usize)> = (0..n).map(|i| (i, (i + 1) % n)).collect();
    let ll = Lowlink::new(&Csr::undirected(n, &cycle));
    assert!(ll.bridges().is_empty() && ll.articulation_points().is_empty());
    assert_eq!(ll.biconnected_components().len(), 1);
    st.ops(2 * n);
    st.done();
}
