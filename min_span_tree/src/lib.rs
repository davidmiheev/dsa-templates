//! Kruskal's minimum spanning forest, `O(m log m)`, on top of [`union_find`].
//!
//! # Example
//!
//! ```
//! use min_span_tree::minimum_spanning_forest;
//!
//! // A triangle with weights 1, 2, 3 plus an isolated vertex.
//! let (total, edges, components) = minimum_spanning_forest(4, &[(0, 1, 1u64), (1, 2, 2), (2, 0, 3)]);
//! assert_eq!((total, components), (3, 2));
//! assert_eq!(edges, vec![0, 1]);
//! ```

use std::ops::Add;

/// A minimum spanning forest of the undirected multigraph `edges = (u, v, weight)` on `0..n`
/// (self-loops and parallel edges allowed). Returns `(total weight, chosen edge indices in
/// increasing weight order, number of components)`; the forest has `n - components` edges.
///
/// # Panics
/// Panics if an endpoint is `>= n`.
pub fn minimum_spanning_forest<W>(n: usize, edges: &[(usize, usize, W)]) -> (W, Vec<usize>, usize)
where
    W: Copy + Ord + Default + Add<Output = W>,
{
    let mut order: Vec<usize> = (0..edges.len()).collect();
    order.sort_unstable_by_key(|&e| edges[e].2);
    // `union_find::plain(k)` covers the indices `0..=k`, one more than needed.
    let mut uf = union_find::plain(n);
    let (mut total, mut chosen) = (W::default(), Vec::new());
    for e in order {
        let (u, v, w) = edges[e];
        assert!(u < n && v < n, "edge {e} has an endpoint outside 0..{n}");
        if uf.merge(u, v) {
            total = total + w;
            chosen.push(e);
        }
    }
    let components = n - chosen.len();
    (total, chosen, components)
}
