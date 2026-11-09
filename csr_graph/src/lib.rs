//! Compressed sparse row (CSR) adjacency: all arcs in two flat arrays, grouped by tail vertex.
//!
//! The vertices are `0..n`; the arcs leaving `v` are the indices `start[v]..start[v + 1]`, with
//! head `to[i]` and the index `id[i]` of the input edge that produced the arc. Compared with
//! `Vec<Vec<usize>>` this is one allocation instead of `n + 1` and walks memory linearly.
//!
//! # Example
//!
//! ```
//! use csr_graph::Csr;
//!
//! let g = Csr::undirected(3, &[(0, 1), (1, 2)]);
//! let around_1: Vec<(usize, usize)> = g.arcs(1).map(|i| (g.to[i], g.id[i])).collect();
//! assert_eq!(around_1, vec![(0, 0), (2, 1)]);
//! assert_eq!(g.degree(1), 2);
//! ```

use std::ops::Range;

/// A directed multigraph in CSR form. See the [module documentation](self).
#[derive(Debug, Clone)]
pub struct Csr {
    /// `start[v]..start[v + 1]` are the arcs leaving `v`; `start.len() == n + 1`.
    pub start: Vec<usize>,
    /// Head of every arc.
    pub to: Vec<usize>,
    /// Index of the originating input edge of every arc.
    pub id: Vec<usize>,
}

impl Csr {
    /// The directed graph with an arc `u -> v` for every edge `(u, v)`. Arcs of a vertex keep the
    /// order of the edge list.
    ///
    /// # Panics
    /// Panics if an endpoint is `>= n`.
    pub fn directed(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut start = vec![0usize; n + 2];
        for &(u, v) in edges {
            assert!(u < n && v < n, "edge ({u}, {v}) has an endpoint outside 0..{n}");
            start[u + 2] += 1;
        }
        for i in 2..n + 2 {
            start[i] += start[i - 1];
        }
        let mut to = vec![0usize; edges.len()];
        let mut id = vec![0usize; edges.len()];
        for (e, &(u, v)) in edges.iter().enumerate() {
            let slot = &mut start[u + 1];
            to[*slot] = v;
            id[*slot] = e;
            *slot += 1;
        }
        start.pop();
        Csr { start, to, id }
    }

    /// The undirected graph: every edge `(u, v)` gives the arcs `u -> v` and `v -> u`, both
    /// carrying the edge index. A self-loop `(v, v)` gives two arcs at `v`.
    ///
    /// # Panics
    /// Panics if an endpoint is `>= n`.
    pub fn undirected(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut start = vec![0usize; n + 2];
        for &(u, v) in edges {
            assert!(u < n && v < n, "edge ({u}, {v}) has an endpoint outside 0..{n}");
            start[u + 2] += 1;
            start[v + 2] += 1;
        }
        for i in 2..n + 2 {
            start[i] += start[i - 1];
        }
        let mut to = vec![0usize; 2 * edges.len()];
        let mut id = vec![0usize; 2 * edges.len()];
        for (e, &(u, v)) in edges.iter().enumerate() {
            let slot = &mut start[u + 1];
            to[*slot] = v;
            id[*slot] = e;
            *slot += 1;
            let slot = &mut start[v + 1];
            to[*slot] = u;
            id[*slot] = e;
            *slot += 1;
        }
        start.pop();
        Csr { start, to, id }
    }

    /// Number of vertices.
    pub fn len(&self) -> usize {
        self.start.len() - 1
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Number of arcs.
    pub fn num_arcs(&self) -> usize {
        self.to.len()
    }

    /// The arc indices leaving `v`.
    pub fn arcs(&self, v: usize) -> Range<usize> {
        self.start[v]..self.start[v + 1]
    }

    /// Number of arcs leaving `v`.
    pub fn degree(&self, v: usize) -> usize {
        self.start[v + 1] - self.start[v]
    }

    /// The heads of the arcs leaving `v`.
    pub fn neighbors(&self, v: usize) -> &[usize] {
        &self.to[self.arcs(v)]
    }

    /// The graph with every arc reversed (edge indices kept).
    pub fn reversed(&self) -> Csr {
        let n = self.len();
        let edges: Vec<(usize, usize)> = (0..n).flat_map(|u| self.arcs(u).map(move |i| (self.to[i], u))).collect();
        let ids: Vec<usize> = (0..n).flat_map(|u| self.arcs(u).map(move |i| self.id[i])).collect();
        let mut g = Csr::directed(n, &edges);
        for slot in g.id.iter_mut() {
            *slot = ids[*slot];
        }
        g
    }
}
