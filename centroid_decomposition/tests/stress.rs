//! Stress test: verify that `build_centroid_tree` returns a genuine centroid decomposition.
//!
//! For every vertex `c`, its component is `c` plus its descendants in the centroid tree. That
//! component must be a connected piece of the original tree, `c` must be one of its centroids
//! (every piece left after removing `c` has at most half the vertices), and `c`'s children in
//! the centroid tree must be exactly one vertex per remaining piece.

use centroid_decomposition::build_centroid_tree;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn verify(adj: &[Vec<usize>], parent: &[Option<usize>], depth: &[usize], ctx: &str) {
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

#[test]
fn produces_valid_decompositions() {
    let mut st = Stress::new("centroid/naive");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 45) as usize;
        let adj = rng.tree(n, case % 4);
        let (parent, depth) = build_centroid_tree(&adj);
        verify(&adj, &parent, &depth, &format!("case {case} n={n} shape={}", case % 4));
        st.case();
        st.ops(n);
    }
    st.done();
}

#[test]
fn large_trees_have_logarithmic_depth() {
    // Recursion depth equals tree depth in `compute_sizes`, so run on a big stack.
    std::thread::Builder::new()
        .stack_size(512 << 20)
        .spawn(|| {
            let mut st = Stress::new("centroid/large");
            for shape in 0..4 {
                let mut rng = Rng::new(200 + shape as u64);
                let n = 200_000;
                let adj = rng.tree(n, shape);
                let (parent, depth) = build_centroid_tree(&adj);
                st.case();
                let bound = (usize::BITS - n.leading_zeros()) as usize;
                assert_eq!(parent.iter().filter(|p| p.is_none()).count(), 1, "one root, shape {shape}");
                assert!(depth.iter().all(|&d| d < bound), "depth < {bound}, shape {shape}");
                for v in 0..n {
                    if let Some(p) = parent[v] {
                        assert_eq!(depth[v], depth[p] + 1);
                    }
                }
                st.ops(n);
            }
            st.done();
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn empty_tree_panics() {
    assert!(std::panic::catch_unwind(|| build_centroid_tree(&[])).is_err());
}
