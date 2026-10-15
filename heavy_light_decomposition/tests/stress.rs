//! Stress test: `Hld` against naive parent/depth/size/path computations on random trees of
//! several shapes (random, path, star, caterpillar, binary), through all three constructors.

use heavy_light_decomposition::Hld;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

/// A random tree on `n` vertices as a parent array (`parent[i] < i`) with a shape, then
/// relabelled by a random permutation so that vertex numbers carry no structure.
fn random_edges(rng: &mut Rng, n: usize, shape: usize) -> Vec<(usize, usize)> {
    let mut perm: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        perm.swap(i, rng.below(i + 1));
    }
    (1..n)
        .map(|i| {
            let p = match shape {
                0 => rng.below(i),                       // random
                1 => i - 1,                              // path
                2 => 0,                                  // star
                3 => if i % 2 == 1 { i / 2 } else { i - 1 }, // caterpillar-ish
                _ => (i - 1) / 2,                        // binary
            };
            (perm[p], perm[i])
        })
        .collect()
}

fn adjacency(n: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut adj = vec![vec![]; n];
    for &(u, v) in edges {
        adj[u].push(v);
        adj[v].push(u);
    }
    adj
}

fn naive(adj: &[Vec<usize>], root: usize) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let n = adj.len();
    let (mut parent, mut depth, mut order) = (vec![Hld::NONE; n], vec![0; n], vec![root]);
    let mut i = 0;
    while i < order.len() {
        let v = order[i];
        i += 1;
        for &u in &adj[v] {
            if u != parent[v] {
                parent[u] = v;
                depth[u] = depth[v] + 1;
                order.push(u);
            }
        }
    }
    let mut size = vec![1; n];
    for &v in order.iter().skip(1).rev() {
        size[parent[v]] += size[v];
    }
    (parent, depth, size)
}

fn naive_path(parent: &[usize], depth: &[usize], mut u: usize, mut v: usize) -> Vec<usize> {
    let (mut left, mut right) = (vec![], vec![]);
    while u != v {
        if depth[u] >= depth[v] {
            left.push(u);
            u = parent[u];
        } else {
            right.push(v);
            v = parent[v];
        }
    }
    left.push(u);
    right.reverse();
    left.extend(right);
    left
}

fn check(hs: &[Hld], adj: &[Vec<usize>], root: usize, rng: &mut Rng, ctx: &str, queries: usize) {
    let n = adj.len();
    let (parent, depth, size) = naive(adj, root);
    for h in hs {
        assert_eq!((&h.parent, &h.depth, &h.size), (&parent, &depth, &size), "parent/depth/size {ctx}");
        assert_eq!(h.len(), n);
        // pos/order are inverse permutations; subtrees are contiguous and start at pos[v].
        for v in 0..n {
            assert_eq!(h.order[h.pos[v]], v, "order inverts pos {ctx}");
            let (lo, hi) = h.subtree(v);
            assert_eq!(hi - lo, size[v]);
            for p in lo..hi {
                assert!(h.is_ancestor(v, h.order[p]), "subtree block {ctx}");
            }
            // The heavy child directly follows its parent.
            if v != root && h.head[v] == h.head[parent[v]] {
                assert_eq!(h.pos[v], h.pos[parent[v]] + 1, "chains are contiguous {ctx}");
            }
        }
        for _ in 0..queries {
            let (u, v) = (rng.below(n), rng.below(n));
            let path = naive_path(&parent, &depth, u, v);
            assert_eq!(h.dist(u, v), path.len() - 1, "dist({u}, {v}) {ctx}");
            let l = *path.iter().min_by_key(|&&x| depth[x]).unwrap();
            assert_eq!(h.lca(u, v), l, "lca({u}, {v}) {ctx}");
            let k = rng.below(path.len());
            assert_eq!(h.path_vertex(u, v, k), path[k], "path_vertex({u}, {v}, {k}) {ctx}");
            let j = rng.below(depth[u] + 1);
            let mut anc = u;
            for _ in 0..j {
                anc = parent[anc];
            }
            assert_eq!(h.jump(u, j), anc, "jump({u}, {j}) {ctx}");
            // Unordered segments cover exactly the path's vertices.
            let mut covered: Vec<usize> = h.path_segments(u, v).iter().flat_map(|&(a, b)| (a..b).map(|p| h.order[p])).collect();
            covered.sort_unstable();
            let mut want = path.clone();
            want.sort_unstable();
            assert_eq!(covered, want, "path_segments({u}, {v}) {ctx}");
            // Ordered segments, walked in their direction, give the path in order.
            let walked: Vec<usize> = h
                .path_ordered(u, v)
                .iter()
                .flat_map(|&(a, b, rev)| {
                    let r: Vec<usize> = (a..b).map(|p| h.order[p]).collect();
                    if rev { r.into_iter().rev().collect::<Vec<_>>() } else { r }
                })
                .collect();
            assert_eq!(walked, path, "path_ordered({u}, {v}) {ctx}");
        }
    }
}

#[test]
fn matches_naive_on_random_trees() {
    let mut st = Stress::new("hld/naive");
    for case in 0..800 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 40) as usize;
        let edges = random_edges(&mut rng, n, case % 5);
        let adj = adjacency(n, &edges);
        let root = rng.below(n);
        let (parent, _, _) = naive(&adj, root);
        let hs = [Hld::new(&adj, root), Hld::from_edges(n, &edges, root), Hld::from_parents(&parent, root)];
        check(&hs, &adj, root, &mut rng, &format!("case {case} n={n}"), 20);
        st.case();
        st.ops(20);
    }
    st.done();
}

#[test]
fn million_vertex_path_does_not_overflow_the_stack() {
    let mut st = Stress::new("hld/large");
    let n = 1_000_000;
    let edges: Vec<(usize, usize)> = (1..n).map(|i| (i - 1, i)).collect();
    let h = Hld::from_edges(n, &edges, 0);
    st.case();
    assert_eq!(h.lca(n - 1, n / 2), n / 2);
    assert_eq!(h.jump(n - 1, n - 1), 0);
    assert_eq!(h.path_segments(0, n - 1), vec![(0, n)]);
    // A random tree with a million vertices.
    let mut rng = Rng::new(3);
    let edges: Vec<(usize, usize)> = (1..n).map(|i| (rng.below(i), i)).collect();
    let h = Hld::from_edges(n, &edges, 0);
    for _ in 0..100_000 * scale() {
        let (u, v) = (rng.below(n), rng.below(n));
        let l = h.lca(u, v);
        assert!(h.is_ancestor(l, u) && h.is_ancestor(l, v));
        assert_eq!(h.dist(u, v), h.depth[u] + h.depth[v] - 2 * h.depth[l]);
        st.ops(1);
    }
    st.done();
}

#[test]
#[should_panic(expected = "not connected")]
fn disconnected_graph_is_rejected() {
    Hld::new(&[vec![1], vec![0], vec![3], vec![2]], 0);
}
