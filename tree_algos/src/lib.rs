//! Small, self-contained tree algorithms.
//!
//! * [`cartesian_tree`]: the min-rooted Cartesian tree of a sequence in `O(n)`.
//! * [`tree_diameter`]: weighted diameter and one longest path in `O(n)`.
//! * [`rooted_isomorphism_classes`]: classes of isomorphic rooted subtrees in `O(n log n)`.
//!
//! Parent arrays use [`NONE`] as the parent of a root.
//!
//! # Example
//!
//! ```
//! use tree_algos::{cartesian_tree, tree_diameter, NONE};
//!
//! assert_eq!(cartesian_tree(&[3, 1, 2]), vec![1, NONE, 1]);
//! let (len, path) = tree_diameter(4, &[(0, 1, 5u64), (1, 2, 1), (1, 3, 7)]);
//! assert_eq!((len, path.len()), (12, 3)); // 0 - 1 - 3
//! ```

use std::collections::HashMap;
use std::ops::Add;

/// The parent of a root.
pub const NONE: usize = usize::MAX;

/// The Cartesian tree of `a` as a parent array: the root holds the smallest element (the
/// leftmost one among equals), and the in-order traversal is `0..n`.
pub fn cartesian_tree<T: Ord>(a: &[T]) -> Vec<usize> {
    let n = a.len();
    let mut parent = vec![NONE; n];
    let mut stack: Vec<usize> = Vec::new();
    for i in 0..n {
        let mut last = NONE;
        while let Some(&top) = stack.last() {
            if a[top] > a[i] {
                last = top;
                stack.pop();
            } else {
                break;
            }
        }
        if last != NONE {
            parent[last] = i;
        }
        if let Some(&top) = stack.last() {
            parent[i] = top;
        }
        stack.push(i);
    }
    parent
}

/// The diameter of a weighted tree on `n` vertices (non-negative weights): its length and the
/// vertices of one longest path, in order. Returns `(0, [0])` for a single vertex.
///
/// # Panics
/// Panics if `n == 0`, an endpoint is out of range, or the edges do not form a tree.
pub fn tree_diameter<W: Copy + Ord + Default + Add<Output = W>>(n: usize, edges: &[(usize, usize, W)]) -> (W, Vec<usize>) {
    assert!(n > 0, "the tree must have at least one vertex");
    assert_eq!(edges.len() + 1, n, "a tree on {n} vertices has {} edges", n - 1);
    let mut offsets = vec![0usize; n + 2];
    for &(u, v, _) in edges {
        assert!(u < n && v < n, "edge endpoint out of range");
        offsets[u + 2] += 1;
        offsets[v + 2] += 1;
    }
    for i in 2..n + 2 {
        offsets[i] += offsets[i - 1];
    }
    let mut adj = vec![(0usize, W::default()); 2 * edges.len()];
    for &(u, v, w) in edges {
        adj[offsets[u + 1]] = (v, w);
        offsets[u + 1] += 1;
        adj[offsets[v + 1]] = (u, w);
        offsets[v + 1] += 1;
    }
    // Farthest vertex from `s`, with parents for path recovery.
    let sweep = |s: usize| -> (usize, W, Vec<usize>) {
        let mut parent = vec![NONE; n];
        let mut dist = vec![W::default(); n];
        let mut order = vec![s];
        let mut i = 0;
        while i < order.len() {
            let v = order[i];
            i += 1;
            for &(u, w) in &adj[offsets[v]..offsets[v + 1]] {
                if u != parent[v] {
                    parent[u] = v;
                    dist[u] = dist[v] + w;
                    order.push(u);
                }
            }
        }
        assert_eq!(order.len(), n, "the graph is not connected");
        let far = *order.iter().max_by_key(|&&v| dist[v]).expect("non-empty");
        (far, dist[far], parent)
    };
    let (a, _, _) = sweep(0);
    let (b, len, parent) = sweep(a);
    let mut path = vec![b];
    while *path.last().unwrap() != a {
        path.push(parent[*path.last().unwrap()]);
    }
    path.reverse();
    (len, path)
}

/// Classify the subtrees of a rooted tree by isomorphism. Returns the number of classes `k` and
/// a class id in `0..k` for the subtree rooted at every vertex; two vertices get the same id
/// exactly when their subtrees are isomorphic as rooted trees.
///
/// # Panics
/// Panics if the parent array does not describe a tree with the single root `root`.
pub fn rooted_isomorphism_classes(parent: &[usize], root: usize) -> (usize, Vec<usize>) {
    let n = parent.len();
    assert!(root < n, "root out of range");
    let mut offsets = vec![0usize; n + 2];
    for (v, &p) in parent.iter().enumerate() {
        if v != root {
            assert!(p < n, "parent out of range");
            offsets[p + 2] += 1;
        }
    }
    for i in 2..n + 2 {
        offsets[i] += offsets[i - 1];
    }
    let mut children = vec![0usize; n - 1];
    for (v, &p) in parent.iter().enumerate() {
        if v != root {
            children[offsets[p + 1]] = v;
            offsets[p + 1] += 1;
        }
    }
    offsets.pop();
    let mut order = vec![root];
    let mut i = 0;
    while i < order.len() {
        let v = order[i];
        i += 1;
        order.extend_from_slice(&children[offsets[v]..offsets[v + 1]]);
    }
    assert_eq!(order.len(), n, "the parent array does not describe a tree rooted at {root}");
    let mut class = vec![0usize; n];
    let mut ids: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut key: Vec<u32> = Vec::new();
    for &v in order.iter().rev() {
        key.clear();
        key.extend(children[offsets[v]..offsets[v + 1]].iter().map(|&c| class[c] as u32));
        key.sort_unstable();
        let next = ids.len();
        class[v] = match ids.get(&key) {
            Some(&id) => id,
            None => {
                ids.insert(key.clone(), next);
                next
            }
        };
    }
    (ids.len(), class)
}
