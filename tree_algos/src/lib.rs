//! Small, self-contained tree algorithms.
//!
//! * [`cartesian_tree`]: the min-rooted Cartesian tree of a sequence in `O(n)`.
//! * [`tree_diameter`]: weighted diameter and one longest path in `O(n)`.
//! * [`rooted_isomorphism_classes`]: classes of isomorphic rooted subtrees in `O(n log n)`.
//! * [`rerooting`]: the answer of a tree DP for *every* choice of root in `O(n)` merges.
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

/// Tree DP for every root at once ("rerooting").
///
/// For a root `r` the DP value of a vertex `v` is
/// `vertex(v, fold over neighbours u (away from r) of edge(dp(u), e, u, v))`, where `fold` is
/// the commutative-or-not monoid `(identity, merge)` (the neighbours of `v` are merged in the
/// order of the edge list, which only matters for non-commutative `merge`), `edge(value, e, u, v)`
/// carries the DP value `value` of the subtree at `u` across the edge `e = edges[e] = {u, v}` to
/// `v`, and `vertex(acc, v)` finishes the value of `v` from the merged contributions.
///
/// Returns `ans[r]` = the DP value of `r` when the tree is rooted at `r`, for every `r`.
/// No inverse of `merge` is needed (prefix/suffix products are used), `O(n)` closure calls.
///
/// ```
/// use tree_algos::rerooting;
///
/// // Farthest vertex distance (eccentricity) of every vertex of the path 0 - 1 - 2 - 3.
/// let edges = [(0, 1), (1, 2), (2, 3)];
/// let ecc = rerooting(4, &edges, 0u32, |a, b| *a.max(b), |&d, _, _, _| d + 1, |&acc, _| acc);
/// assert_eq!(ecc, vec![3, 2, 2, 3]);
/// ```
///
/// # Panics
/// Panics if the edges do not form a tree on `n >= 1` vertices.
pub fn rerooting<V, Merge, Edge, Vertex>(n: usize, edges: &[(usize, usize)], identity: V, merge: Merge, edge: Edge, vertex: Vertex) -> Vec<V>
where
    V: Clone,
    Merge: Fn(&V, &V) -> V,
    Edge: Fn(&V, usize, usize, usize) -> V,
    Vertex: Fn(&V, usize) -> V,
{
    assert!(n > 0, "the tree must have at least one vertex");
    assert_eq!(edges.len() + 1, n, "a tree on {n} vertices has {} edges", n - 1);
    // Adjacency in CSR form: (neighbour, edge index).
    let mut offsets = vec![0usize; n + 2];
    for &(u, v) in edges {
        assert!(u < n && v < n, "edge endpoint out of range");
        offsets[u + 2] += 1;
        offsets[v + 2] += 1;
    }
    for i in 2..n + 2 {
        offsets[i] += offsets[i - 1];
    }
    let mut adj = vec![(0usize, 0usize); 2 * edges.len()];
    for (e, &(u, v)) in edges.iter().enumerate() {
        adj[offsets[u + 1]] = (v, e);
        offsets[u + 1] += 1;
        adj[offsets[v + 1]] = (u, e);
        offsets[v + 1] += 1;
    }
    offsets.pop();
    // Root at 0.
    let mut parent = vec![NONE; n];
    let mut parent_edge = vec![NONE; n];
    let mut order = vec![0usize];
    let mut i = 0;
    while i < order.len() {
        let v = order[i];
        i += 1;
        for &(u, e) in &adj[offsets[v]..offsets[v + 1]] {
            if u != parent[v] {
                parent[u] = v;
                parent_edge[u] = e;
                order.push(u);
            }
        }
    }
    assert_eq!(order.len(), n, "the graph is not connected");
    // Bottom-up: `up_value[c]` = contribution of the subtree of `c` to its parent.
    let mut up_value: Vec<Option<V>> = vec![None; n];
    for &v in order.iter().rev() {
        let mut acc = identity.clone();
        for &(u, _) in &adj[offsets[v]..offsets[v + 1]] {
            if u != parent[v] {
                acc = merge(&acc, up_value[u].as_ref().expect("children come first"));
            }
        }
        let value = vertex(&acc, v);
        if parent[v] != NONE {
            up_value[v] = Some(edge(&value, parent_edge[v], v, parent[v]));
        }
    }
    // Top-down: `from_parent[c]` = contribution of everything outside the subtree of `c`.
    let mut from_parent: Vec<Option<V>> = vec![None; n];
    let mut ans: Vec<Option<V>> = vec![None; n];
    for &v in &order {
        let nbrs = &adj[offsets[v]..offsets[v + 1]];
        // Contributions in adjacency order (the parent's slot holds `from_parent[v]`).
        let contrib: Vec<&V> = nbrs.iter().map(|&(u, _)| if u == parent[v] { from_parent[v].as_ref().expect("parent first") } else { up_value[u].as_ref().expect("computed") }).collect();
        let mut suffix = vec![identity.clone(); contrib.len() + 1];
        for k in (0..contrib.len()).rev() {
            suffix[k] = merge(contrib[k], &suffix[k + 1]);
        }
        let mut prefix = identity.clone();
        let mut sent: Vec<(usize, V)> = Vec::new();
        for (k, &(u, e)) in nbrs.iter().enumerate() {
            if u != parent[v] {
                let rest = vertex(&merge(&prefix, &suffix[k + 1]), v);
                sent.push((u, edge(&rest, e, v, u)));
            }
            prefix = merge(&prefix, contrib[k]);
        }
        ans[v] = Some(vertex(&prefix, v));
        drop(contrib);
        for (u, value) in sent {
            from_parent[u] = Some(value);
        }
    }
    ans.into_iter().map(|a| a.expect("every vertex is reached")).collect()
}
