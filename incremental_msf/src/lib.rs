//! Online minimum spanning forest under edge insertions in `O(log n)` amortised per edge.
//!
//! Every edge is a node of a [`link_cut_tree::LinkCutTree`] sitting between its two endpoints, so
//! "the heaviest edge on the tree path between `u` and `v`" is a path query over a max-monoid.
//! Inserting `{u, v}`: if `u` and `v` are in different trees the edge joins them; otherwise the
//! heaviest edge on the cycle it would close — possibly the new edge itself — is dropped.
//!
//! # Example
//!
//! ```
//! use incremental_msf::{IncrementalMsf, NONE};
//!
//! let mut f = IncrementalMsf::new(3, 4);
//! assert_eq!(f.add_edge(0, 1, 5i64), NONE); // edge 0: joins two trees
//! assert_eq!(f.add_edge(1, 2, 7), NONE); // edge 1
//! assert_eq!(f.add_edge(0, 2, 6), 1); // edge 2 replaces the heaviest cycle edge, edge 1
//! assert_eq!(f.add_edge(0, 1, 9), 3); // edge 3 is itself the heaviest, so it is dropped
//! ```

use algebra::Monoid;
use link_cut_tree::LinkCutTree;
use std::marker::PhantomData;

/// "No edge left the forest" in the result of [`IncrementalMsf::add_edge`].
pub const NONE: usize = usize::MAX;

/// Max over `(weight, edge index)`; vertex nodes carry `None`.
struct MaxEdge<W>(PhantomData<W>);

impl<W: Copy + Ord> Monoid for MaxEdge<W> {
    type S = Option<(W, usize)>;
    fn identity() -> Self::S {
        None
    }
    fn op(a: &Self::S, b: &Self::S) -> Self::S {
        match (a, b) {
            (Some(x), Some(y)) => Some(if y.0 > x.0 { *y } else { *x }),
            (Some(_), None) => *a,
            _ => *b,
        }
    }
}

/// A minimum spanning forest that grows by edge insertions.
pub struct IncrementalMsf<W: Copy + Ord> {
    n: usize,
    tree: LinkCutTree<MaxEdge<W>>,
    ends: Vec<(usize, usize)>,
    weights: Vec<W>,
}

impl<W: Copy + Ord> IncrementalMsf<W> {
    /// A forest of `n` isolated vertices that will receive at most `max_edges` insertions (the
    /// link-cut tree has one node per vertex and per edge, allocated up front).
    pub fn new(n: usize, max_edges: usize) -> Self {
        Self {
            n,
            tree: LinkCutTree::new(&vec![None; n + max_edges]),
            ends: Vec::with_capacity(max_edges),
            weights: Vec::with_capacity(max_edges),
        }
    }

    /// Insert the edge `{u, v}` of weight `w`; it gets the next index (`0, 1, ...`). Returns the
    /// index of the edge that is no longer in the forest afterwards (the new edge itself if it is
    /// the heaviest edge of the cycle it closes), or [`NONE`] if the forest just grew.
    ///
    /// # Panics
    /// Panics if more edges are inserted than `max_edges`, or a vertex is out of range.
    pub fn add_edge(&mut self, u: usize, v: usize, w: W) -> usize {
        assert!(u < self.n && v < self.n, "vertex out of range");
        let id = self.ends.len();
        assert!(self.n + id < self.tree.len(), "more edges than the capacity {}", self.tree.len() - self.n);
        self.ends.push((u, v));
        self.weights.push(w);
        if u == v {
            return id;
        }
        let node = self.n + id;
        if self.tree.connected(u, v) {
            let (heaviest, j) = self.tree.path_prod(u, v).expect("a path between distinct vertices has an edge");
            if heaviest <= w {
                return id;
            }
            let (a, b) = self.ends[j];
            let jn = self.n + j;
            self.tree.cut(a, jn);
            self.tree.cut(jn, b);
            self.attach(node, u, v, w, id);
            j
        } else {
            self.attach(node, u, v, w, id);
            NONE
        }
    }

    fn attach(&mut self, node: usize, u: usize, v: usize, w: W, id: usize) {
        self.tree.set(node, Some((w, id)));
        self.tree.link(u, node);
        self.tree.link(node, v);
    }

    /// Whether `u` and `v` are currently connected.
    pub fn connected(&mut self, u: usize, v: usize) -> bool {
        self.tree.connected(u, v)
    }
}
