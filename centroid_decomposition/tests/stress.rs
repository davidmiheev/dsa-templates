//! Stress test: verify that `build_centroid_tree` returns a genuine centroid decomposition.
//!
//! For every vertex `c`, its component is `c` plus its descendants in the centroid tree. That
//! component must be a connected piece of the original tree, `c` must be one of its centroids
//! (every piece left after removing `c` has at most half the vertices), and `c`'s children in
//! the centroid tree must be exactly one vertex per remaining piece.

use centroid_decomposition::{CentroidTree, ContourAdd, ContourSum};
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn verify(adj: &[Vec<usize>], parent_raw: &[usize], depth: &[usize], ctx: &str) {
    let parent: Vec<Option<usize>> = parent_raw.iter().map(|&p| if p == CentroidTree::NONE { None } else { Some(p) }).collect();
    let n = adj.len();
    assert_eq!((parent.len(), depth.len()), (n, n));
    let roots: Vec<usize> = (0..n).filter(|&v| parent[v].is_none()).collect();
    assert_eq!(roots.len(), 1, "exactly one centroid-tree root {ctx}");
    assert_eq!(depth[roots[0]], 0);
    for v in 0..n {
        if let Some(p) = parent[v] {
            assert_eq!(depth[v], depth[p] + 1, "depth of {v} {ctx}");
        }
    }
    // children lists and descendant sets (n is small in this check)
    let mut kids = vec![vec![]; n];
    for v in 0..n {
        if let Some(p) = parent[v] {
            kids[p].push(v);
        }
    }
    let mut comp: Vec<Vec<usize>> = vec![vec![]; n];
    for c in 0..n {
        let mut stack = vec![c];
        while let Some(x) = stack.pop() {
            comp[c].push(x);
            stack.extend(&kids[x]);
        }
    }
    assert_eq!(comp[roots[0]].len(), n, "the root's component is the whole tree {ctx}");
    for c in 0..n {
        let inside: Vec<bool> = (0..n).map(|v| comp[c].contains(&v)).collect();
        // pieces of the component after removing c
        let mut piece_of = vec![usize::MAX; n];
        let mut pieces: Vec<usize> = vec![];
        for &start in &comp[c] {
            if start == c || piece_of[start] != usize::MAX {
                continue;
            }
            let id = pieces.len();
            let (mut stack, mut size) = (vec![start], 0);
            piece_of[start] = id;
            while let Some(x) = stack.pop() {
                size += 1;
                for &u in &adj[x] {
                    if inside[u] && u != c && piece_of[u] == usize::MAX {
                        piece_of[u] = id;
                        stack.push(u);
                    }
                }
            }
            assert!(size <= comp[c].len() / 2, "piece of size {size} after removing {c} from a component of {} {ctx}", comp[c].len());
            pieces.push(size);
        }
        // the component (including c) must be connected: c touches every piece
        for id in 0..pieces.len() {
            assert!(adj[c].iter().any(|&u| inside[u] && piece_of[u] == id), "component of {c} must be connected {ctx}");
        }
        // one child per piece, each inside its own piece
        let mut child_pieces: Vec<usize> = kids[c].iter().map(|&k| piece_of[k]).collect();
        child_pieces.sort_unstable();
        assert_eq!(child_pieces, (0..pieces.len()).collect::<Vec<_>>(), "children of {c} <-> pieces {ctx}");
    }
    let bound = (usize::BITS - n.leading_zeros()) as usize; // floor(log2 n) + 1
    assert!(depth.iter().all(|&d| d < bound), "depth must be O(log n) {ctx}");
}

fn bfs_dist(adj: &[Vec<usize>], s: usize) -> Vec<usize> {
    let mut d = vec![usize::MAX; adj.len()];
    d[s] = 0;
    let mut q = vec![s];
    let mut i = 0;
    while i < q.len() {
        let v = q[i];
        i += 1;
        for &u in &adj[v] {
            if d[u] == usize::MAX {
                d[u] = d[v] + 1;
                q.push(u);
            }
        }
    }
    d
}

#[test]
fn produces_valid_decompositions_and_distances() {
    let mut st = Stress::new("centroid/naive");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 45) as usize;
        let adj = rng.tree(n, case % 4);
        let t = CentroidTree::new(&adj);
        verify(&adj, &t.parent, &t.depth, &format!("case {case} n={n} shape={}", case % 4));
        // Distances to every centroid ancestor, radius, and the iterator.
        let all: Vec<Vec<usize>> = (0..n).map(|s| bfs_dist(&adj, s)).collect();
        for v in 0..n {
            let mut chain = vec![];
            let mut c = v;
            while c != CentroidTree::NONE {
                chain.push(c);
                c = t.parent[c];
            }
            let got: Vec<(usize, usize)> = t.ancestors(v).collect();
            let want: Vec<(usize, usize)> = chain.iter().map(|&c| (c, all[v][c])).collect();
            assert_eq!(got, want, "ancestors of {v}, case {case}");
            for (k, &c) in chain.iter().enumerate() {
                assert_eq!(t.dist_to_ancestor(v, t.depth[v] - k), all[v][c]);
            }
        }
        // radius[c] = farthest vertex in c's component.
        for c in 0..n {
            let comp_far = (0..n).filter(|&u| { let mut x = u; loop { if x == c { break true; } if x == CentroidTree::NONE { break false; } x = t.parent[x]; } }).map(|u| all[c][u]).max().unwrap();
            assert_eq!(t.radius[c], comp_far, "radius of {c}, case {case}");
        }
        assert!(catch_unwind(AssertUnwindSafe(|| t.dist_to_ancestor(t.root, 1))).is_err());
        st.case();
        st.ops(n);
    }
    st.done();
}

#[test]
fn contour_sum_matches_bfs() {
    let mut st = Stress::new("centroid/contour_sum");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64 + 100);
        let n = rng.range(1, 40) as usize;
        let adj = rng.tree(n, case % 4);
        let all: Vec<Vec<usize>> = (0..n).map(|s| bfs_dist(&adj, s)).collect();
        let mut c = ContourSum::<i64>::new(CentroidTree::new(&adj));
        let mut val = vec![0i64; n];
        st.case();
        for step in 0..100 {
            if rng.chance(1, 2) {
                let (v, x) = (rng.below(n), rng.range(-1000, 1000));
                c.add(v, x);
                val[v] += x;
            } else {
                let (p, l) = (rng.below(n), rng.below(n + 2));
                let r = l + rng.below(n + 2 - l);
                let want: i64 = (0..n).filter(|&u| l <= all[p][u] && all[p][u] < r).map(|u| val[u]).sum();
                assert_eq!(c.sum(p, l, r), want, "sum({p}, {l}, {r}) case {case} step {step}");
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn contour_add_matches_bfs() {
    let mut st = Stress::new("centroid/contour_add");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64 + 300);
        let n = rng.range(1, 40) as usize;
        let adj = rng.tree(n, case % 4);
        let all: Vec<Vec<usize>> = (0..n).map(|s| bfs_dist(&adj, s)).collect();
        let mut c = ContourAdd::<i64>::new(CentroidTree::new(&adj));
        let mut val = vec![0i64; n];
        st.case();
        for step in 0..100 {
            if rng.chance(1, 2) {
                let (p, l) = (rng.below(n), rng.below(n + 2));
                let r = l + rng.below(n + 2 - l);
                let x = rng.range(-1000, 1000);
                c.add(p, l, r, x);
                for u in 0..n {
                    if l <= all[p][u] && all[p][u] < r {
                        val[u] += x;
                    }
                }
            } else {
                let v = rng.below(n);
                assert_eq!(c.get(v), val[v], "get({v}) case {case} step {step}");
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn large_trees_have_logarithmic_depth_and_work() {
    let mut st = Stress::new("centroid/large");
    for shape in 0..4 {
        let mut rng = Rng::new(200 + shape as u64);
        let n = 200_000;
        let adj = rng.tree(n, shape);
        let t = CentroidTree::new(&adj);
        st.case();
        let bound = (usize::BITS - n.leading_zeros()) as usize;
        assert_eq!(t.parent.iter().filter(|&&p| p == CentroidTree::NONE).count(), 1, "one root, shape {shape}");
        assert!(t.depth.iter().all(|&d| d < bound), "depth < {bound}, shape {shape}");
        let mut sums = ContourSum::<i64>::new(t);
        for i in 0..20_000 * scale() {
            let v = rng.below(n);
            sums.add(v, 1);
            let _ = sums.sum(rng.below(n), 0, 1 + i % 50);
        }
        assert_eq!(sums.sum(0, 0, n + 1), 20_000 * scale() as i64);
        st.ops(n);
    }
    st.done();
}

#[test]
fn bad_inputs_panic() {
    assert!(catch_unwind(|| CentroidTree::new(&[])).is_err());
    assert!(catch_unwind(|| CentroidTree::from_edges(3, &[(0, 1), (0, 1)])).is_err(), "disconnected / cyclic");
}

#[test]
fn pair_counts_by_distance_match_all_pairs() {
    let naive = |a: &[u64], b: &[u64]| {
        let mut c = vec![0u64; a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                c[i + j] += x * y;
            }
        }
        c
    };
    let mut st = Stress::new("centroid/pair_counts");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64 + 900);
        let n = rng.range(1, 40) as usize;
        let adj = rng.tree(n, case % 4);
        let all: Vec<Vec<usize>> = (0..n).map(|s| bfs_dist(&adj, s)).collect();
        let mut want = vec![0u64; n];
        for u in 0..n {
            for v in u + 1..n {
                want[all[u][v]] += 1;
            }
        }
        let t = CentroidTree::new(&adj);
        assert_eq!(t.pair_counts_by_distance(naive), want, "case {case} n={n}");
        st.case();
        st.ops(n * n);
    }
    st.done();
}
