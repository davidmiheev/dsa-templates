//! Centroid Decomposition of a tree.
//!
//! Centroid decomposition recursively splits a tree at its centroid —
//! the vertex whose removal leaves components of size at most `n/2` —
//! and repeats on each component. The resulting recursion tree has
//! height O(log n), and any two original vertices share an ancestor
//! in the decomposition at depth at most O(log n).
//!
//! This is the foundation for "small-to-large" tricks on trees, for
//! answering distance queries offline (count paths of length ≤ k), and
//! for many divide-and-conquer-on-tree problems.
//!
//! This crate exposes [`build_centroid_tree`], which returns the
//! centroid tree: for each vertex, the centroid chosen at the next
//! level up (its "centroid parent") and the depth of the vertex in the
//! centroid tree.
//!
//! # Example
//!
//! ```
//! use centroid_decomposition::build_centroid_tree;
//!
//! // Path 0 -- 1 -- 2 -- 3 -- 4
//! let adj = vec![
//!     vec![1],
//!     vec![0, 2],
//!     vec![1, 3],
//!     vec![2, 4],
//!     vec![3],
//! ];
//! let (parent, depth) = build_centroid_tree(&adj);
//! // Centroid of a 5-path is the middle vertex (2).
//! assert_eq!(parent[2], None);
//! ```

/// Build the centroid decomposition of the tree `adj`.
///
/// Returns two vectors of length `n`:
/// * `parent[v]` — the centroid at the level above `v` in the
///   decomposition (i.e., the centroid chosen for the component that
///   contained `v` in the previous recursion step). `None` for the root
///   of the centroid tree.
/// * `depth[v]` — depth of `v` in the centroid tree (0 for the root).
///
/// # Panics
/// Panics if `adj` is empty.
pub fn build_centroid_tree(adj: &[Vec<usize>]) -> (Vec<Option<usize>>, Vec<usize>) {
    let n = adj.len();
    assert!(n > 0, "graph must be non-empty");

    let mut parent: Vec<Option<usize>> = vec![None; n];
    let mut depth = vec![0usize; n];
    let mut removed = vec![false; n];
    let mut size = vec![0usize; n];

    decompose(adj, &mut removed, &mut size, 0, None, 0, &mut parent, &mut depth);

    (parent, depth)
}

fn decompose(
    adj: &[Vec<usize>],
    removed: &mut [bool],
    size: &mut [usize],
    entry: usize,
    centroid_parent: Option<usize>,
    centroid_depth: usize,
    parent: &mut [Option<usize>],
    depth: &mut [usize],
) {
    // Compute subtree sizes within the current component, rooted at
    // `entry`.
    let root = entry;
    compute_sizes(adj, removed, size, root, usize::MAX);
    let total = size[root];

    // Find centroid: a vertex `c` in the current component where every
    // child component has size ≤ total/2.
    let c = find_centroid(adj, removed, size, root, usize::MAX, total);

    parent[c] = centroid_parent;
    depth[c] = centroid_depth;
    removed[c] = true;

    // Recurse on each component formed by removing `c`.
    for &u in &adj[c] {
        if !removed[u] {
            decompose(adj, removed, size, u, Some(c), centroid_depth + 1, parent, depth);
        }
    }
    removed[c] = false;
}

fn compute_sizes(adj: &[Vec<usize>], removed: &[bool], size: &mut [usize], v: usize, p: usize) -> usize {
    let mut s = 1usize;
    for &u in &adj[v] {
        if u != p && !removed[u] {
            s += compute_sizes(adj, removed, size, u, v);
        }
    }
    size[v] = s;
    s
}

fn find_centroid(
    adj: &[Vec<usize>],
    removed: &[bool],
    size: &[usize],
    v: usize,
    p: usize,
    total: usize,
) -> usize {
    for &u in &adj[v] {
        if u != p && !removed[u] && size[u] > total / 2 {
            return find_centroid(adj, removed, size, u, v, total);
        }
    }
    v
}
