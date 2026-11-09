//! Stress test: Lengauer–Tarjan against "remove a vertex and re-run reachability".

use dominator_tree::{DominatorTree, NONE};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn reach(n: usize, adj: &[Vec<usize>], root: usize, banned: usize) -> Vec<bool> {
    let mut seen = vec![false; n];
    if root == banned {
        return seen;
    }
    seen[root] = true;
    let mut st = vec![root];
    while let Some(v) = st.pop() {
        for &w in &adj[v] {
            if w != banned && !seen[w] {
                seen[w] = true;
                st.push(w);
            }
        }
    }
    seen
}

#[test]
fn matches_removal_brute_force() {
    let mut st = Stress::new("dominator_tree");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 9) as usize;
        let m = rng.range(0, 20) as usize;
        let edges: Vec<(usize, usize)> = (0..m).map(|_| (rng.below(n), rng.below(n))).collect();
        let root = rng.below(n);
        let mut adj = vec![Vec::new(); n];
        for &(a, b) in &edges {
            adj[a].push(b);
        }
        let d = DominatorTree::from_edges(n, &edges, root);
        let base = reach(n, &adj, root, NONE);
        // dom[v] = set of strict dominators
        let mut doms = vec![vec![false; n]; n];
        for x in 0..n {
            let r = reach(n, &adj, root, x);
            for v in 0..n {
                if base[v] && v != x && !r[v] {
                    doms[v][x] = true;
                }
            }
        }
        for v in 0..n {
            if !base[v] {
                assert_eq!(d.idom[v], NONE, "case {case}");
                continue;
            }
            if v == root {
                assert_eq!(d.idom[v], root);
                continue;
            }
            // idom: the strict dominator that has the most strict dominators itself.
            let best = (0..n)
                .filter(|&x| doms[v][x])
                .max_by_key(|&x| doms[x].iter().filter(|&&b| b).count())
                .unwrap_or(root);
            let best = if doms[v].iter().all(|&b| !b) { root } else { best };
            assert_eq!(d.idom[v], best, "case {case}: {edges:?} root {root} v {v}");
            for x in 0..n {
                if base[x] {
                    assert_eq!(d.dominates(x, v), x == v || doms[v][x], "case {case}");
                }
            }
            st.ops(1);
        }
        st.case();
    }
    st.done();
}
