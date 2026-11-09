//! Stress test: directed MST against exhaustive search over every choice of incoming edges.

use directed_mst::{directed_mst, NONE};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

/// Whether `parent` (for non-roots) forms a tree rooted at `root`.
fn is_arborescence(n: usize, root: usize, parent: &[usize]) -> bool {
    for v in 0..n {
        let (mut x, mut steps) = (v, 0);
        while x != root {
            x = parent[x];
            steps += 1;
            if steps > n {
                return false;
            }
        }
    }
    true
}

fn brute(n: usize, root: usize, edges: &[(usize, usize, i64)]) -> Option<i64> {
    let mut by_head: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (i, &(_, b, _)) in edges.iter().enumerate() {
        by_head[b].push(i);
    }
    let mut best: Option<i64> = None;
    let mut pick = vec![0usize; n];
    fn rec(v: usize, n: usize, root: usize, by_head: &[Vec<usize>], edges: &[(usize, usize, i64)], pick: &mut Vec<usize>, best: &mut Option<i64>) {
        if v == n {
            let parent: Vec<usize> = (0..n).map(|x| if x == root { root } else { edges[pick[x]].0 }).collect();
            if is_arborescence(n, root, &parent) {
                let c: i64 = (0..n).filter(|&x| x != root).map(|x| edges[pick[x]].2).sum();
                if best.map_or(true, |b| c < b) {
                    *best = Some(c);
                }
            }
            return;
        }
        if v == root {
            return rec(v + 1, n, root, by_head, edges, pick, best);
        }
        for &e in &by_head[v] {
            pick[v] = e;
            rec(v + 1, n, root, by_head, edges, pick, best);
        }
    }
    rec(0, n, root, &by_head, edges, &mut pick, &mut best);
    best
}

#[test]
fn matches_exhaustive_search() {
    let mut st = Stress::new("directed_mst");
    for case in 0..4000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 6) as usize;
        let m = rng.range(0, 14) as usize;
        let edges: Vec<(usize, usize, i64)> = (0..m)
            .map(|_| (rng.below(n), rng.below(n), rng.range(0, 9)))
            .filter(|&(a, b, _)| a != b)
            .collect();
        let root = rng.below(n);
        let Some(expect) = brute(n, root, &edges) else { continue };
        let t = directed_mst(n, root, &edges);
        assert_eq!(t.cost, expect, "case {case}: {edges:?} root {root}");
        assert_eq!(t.parent[root], root);
        assert_eq!(t.edge[root], NONE);
        assert!(is_arborescence(n, root, &t.parent), "case {case}: not a tree");
        let mut total = 0;
        for v in (0..n).filter(|&v| v != root) {
            let (a, b, w) = edges[t.edge[v]];
            assert_eq!((a, b), (t.parent[v], v));
            total += w;
        }
        assert_eq!(total, expect, "case {case}: chosen edges do not sum to cost");
        st.case();
        st.ops(n);
    }
    st.done();
}

#[test]
#[should_panic]
fn unreachable_vertex_panics() {
    directed_mst(3, 0, &[(0, 1, 1i64)]);
}
