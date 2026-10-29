//! Stress tests of the LCA structures against naive parent walking.
#![allow(clippy::needless_range_loop)]

use lca::{BinaryLifting, Lca, LevelAncestor, RootedTree};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

/// A random tree on `0..n` with randomly permuted labels: returns (edges, parent array for root `root`).
fn random_tree(rng: &mut Rng, n: usize, shape: usize) -> Vec<(usize, usize)> {
    let mut perm: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        perm.swap(i, rng.below(i + 1));
    }
    (1..n)
        .map(|v| {
            let p = match shape {
                0 => rng.below(v),                       // random recursive tree
                1 => v - 1,                              // path
                2 => 0,                                  // star
                3 => (v - 1) / 2,                        // binary heap shape
                4 => v.saturating_sub(1 + rng.below(3)), // almost a path
                _ => v - 1 - rng.below(v.min(8)),        // bushy near the end
            };
            (perm[v], perm[p])
        })
        .collect()
}

fn naive_parents(n: usize, edges: &[(usize, usize)], root: usize) -> (Vec<usize>, Vec<usize>) {
    let mut adj = vec![Vec::new(); n];
    for &(u, v) in edges {
        adj[u].push(v);
        adj[v].push(u);
    }
    let (mut parent, mut depth) = (vec![usize::MAX; n], vec![0usize; n]);
    let mut stack = vec![root];
    let mut seen = vec![false; n];
    seen[root] = true;
    while let Some(v) = stack.pop() {
        for &u in &adj[v] {
            if !seen[u] {
                seen[u] = true;
                parent[u] = v;
                depth[u] = depth[v] + 1;
                stack.push(u);
            }
        }
    }
    (parent, depth)
}

fn naive_lca(parent: &[usize], depth: &[usize], mut u: usize, mut v: usize) -> usize {
    while depth[u] > depth[v] {
        u = parent[u];
    }
    while depth[v] > depth[u] {
        v = parent[v];
    }
    while u != v {
        u = parent[u];
        v = parent[v];
    }
    u
}

#[test]
fn lca_and_level_ancestors() {
    let mut st = Stress::new("lca/structures");
    for case in 0..2500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = 1 + rng.below(if case % 25 == 0 { 400 } else { 40 });
        let shape = rng.below(6);
        let edges = random_tree(&mut rng, n, shape);
        let root = rng.below(n);
        let (parent, depth) = naive_parents(n, &edges, root);
        let t = RootedTree::from_edges(n, &edges, root);
        let (l, la, bl) = (Lca::new(&t), LevelAncestor::new(&t), BinaryLifting::new(&t));
        // parents given as an array build the same tree
        let mut pa = parent.clone();
        pa[root] = root;
        let t2 = RootedTree::from_parents(&pa, root);
        for v in 0..n {
            assert_eq!(t.depth(v), depth[v]);
            assert_eq!(t2.depth(v), depth[v]);
            assert_eq!(
                t.parent(v),
                if v == root {
                    RootedTree::NONE
                } else {
                    parent[v]
                }
            );
        }
        // preorder: subtree ranges are contiguous and contain exactly the descendants
        for v in 0..n {
            let (lo, hi) = l.subtree(v);
            assert_eq!(
                hi - lo,
                (0..n)
                    .filter(|&u| naive_lca(&parent, &depth, u, v) == v)
                    .count()
            );
            for u in 0..n {
                let desc = naive_lca(&parent, &depth, u, v) == v;
                assert_eq!(l.is_ancestor(v, u), desc);
                assert_eq!((lo..hi).contains(&l.tin(u)), desc);
            }
            assert_eq!(l.vertex_at(l.tin(v)), v);
        }
        for u in 0..n {
            for v in 0..n {
                let want = naive_lca(&parent, &depth, u, v);
                assert_eq!(l.lca(u, v), want, "lca({u}, {v}) n={n}");
                assert_eq!(bl.lca(u, v), want);
                assert_eq!(l.dist(u, v), depth[u] + depth[v] - 2 * depth[want]);
            }
            let mut a = u;
            for k in 0..=depth[u] {
                assert_eq!(la.kth_ancestor(u, k), a);
                assert_eq!(bl.kth_ancestor(u, k), a);
                if k < depth[u] {
                    a = parent[a];
                }
            }
        }
        // k-th vertex on a path: walk both halves
        for _ in 0..10 {
            let (u, v) = (rng.below(n), rng.below(n));
            let w = naive_lca(&parent, &depth, u, v);
            let mut path = Vec::new();
            let mut a = u;
            while a != w {
                path.push(a);
                a = parent[a];
            }
            path.push(w);
            let mut tail = Vec::new();
            let mut b = v;
            while b != w {
                tail.push(b);
                b = parent[b];
            }
            path.extend(tail.iter().rev());
            for (k, &x) in path.iter().enumerate() {
                assert_eq!(la.path_vertex(&l, u, v, k), x);
            }
        }
        st.case();
        st.ops(n * n);
    }
    // parent[v] < v input, a deep path and a large random tree
    let mut rng = Rng::new(99);
    for shape in [0usize, 1, 2, 3, 4] {
        let n = 20_000;
        let mut parent = vec![0usize; n];
        for v in 1..n {
            parent[v] = match shape {
                0 => rng.below(v),
                1 => v - 1,
                2 => 0,
                3 => (v - 1) / 2,
                _ => v.saturating_sub(1 + rng.below(3)),
            };
        }
        let t = RootedTree::from_increasing_parents(&parent);
        let (l, la, bl) = (Lca::new(&t), LevelAncestor::new(&t), BinaryLifting::new(&t));
        let (par, depth) = {
            let mut d = vec![0usize; n];
            for v in 1..n {
                d[v] = d[parent[v]] + 1;
            }
            (parent.clone(), d)
        };
        for _ in 0..3000 {
            let (u, v) = (rng.below(n), rng.below(n));
            let want = naive_lca(&par, &depth, u, v);
            assert_eq!(l.lca(u, v), want);
            assert_eq!(bl.lca(u, v), want);
            let k = rng.below(depth[u] + 1);
            let mut a = u;
            for _ in 0..k {
                a = par[a];
            }
            assert_eq!(la.kth_ancestor(u, k), a);
        }
    }
    st.done();
}
