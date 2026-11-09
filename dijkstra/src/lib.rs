//! Dijkstra's single-source shortest paths with non-negative weights, `O((n + m) log n)`.
//!
//! The graph is a [`Csr`](csr_graph::Csr) and the weight of an arc is looked up by its edge
//! index (`weight[g.id[arc]]`), so the same weight array serves directed and undirected graphs.
//! The result stores distances, the shortest-path tree and reachability; asking for the path to
//! an unreachable vertex is a precondition violation (asserted), so check
//! [`ShortestPaths::reachable`] first.
//!
//! # Example
//!
//! ```
//! use csr_graph::Csr;
//! use dijkstra::dijkstra;
//!
//! // 0 -> 1 (4), 0 -> 2 (1), 2 -> 1 (2), 3 isolated.
//! let edges = [(0, 1), (0, 2), (2, 1)];
//! let g = Csr::directed(4, &edges);
//! let sp = dijkstra(&g, &[4u64, 1, 2], 0);
//! assert_eq!(sp.dist[1], 3);
//! assert!(!sp.reachable(3));
//! assert_eq!(sp.path_edges(1), vec![1, 2]);
//! assert_eq!(sp.path_vertices(1), vec![0, 2, 1]);
//! ```

use csr_graph::Csr;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ops::Add;

/// The result of a shortest-path search. See the [module documentation](self).
#[derive(Debug, Clone)]
pub struct ShortestPaths<W> {
    /// Distance from the source; meaningful only where [`reachable`](Self::reachable).
    pub dist: Vec<W>,
    reached: Vec<bool>,
    /// Predecessor on a shortest path ([`ShortestPaths::NONE`] for the source and unreachable
    /// vertices).
    pub parent: Vec<usize>,
    /// The edge used to enter each vertex ([`ShortestPaths::NONE`] likewise).
    pub parent_edge: Vec<usize>,
}

impl<W: Copy> ShortestPaths<W> {
    /// "No parent".
    pub const NONE: usize = usize::MAX;

    /// Whether the source reaches `v`.
    pub fn reachable(&self, v: usize) -> bool {
        self.reached[v]
    }

    /// The edge indices of a shortest path from the source to `t`, in order.
    ///
    /// # Panics
    /// Panics if `t` is not reachable.
    pub fn path_edges(&self, t: usize) -> Vec<usize> {
        assert!(self.reached[t], "vertex {t} is not reachable from the source");
        let mut edges = Vec::new();
        let mut v = t;
        while self.parent[v] != Self::NONE {
            edges.push(self.parent_edge[v]);
            v = self.parent[v];
        }
        edges.reverse();
        edges
    }

    /// The vertices of a shortest path from the source to `t`, in order (both ends included).
    ///
    /// # Panics
    /// Panics if `t` is not reachable.
    pub fn path_vertices(&self, t: usize) -> Vec<usize> {
        assert!(self.reached[t], "vertex {t} is not reachable from the source");
        let mut vertices = vec![t];
        let mut v = t;
        while self.parent[v] != Self::NONE {
            v = self.parent[v];
            vertices.push(v);
        }
        vertices.reverse();
        vertices
    }
}

/// Shortest distances from `source` in the graph `g` with edge weights `weight[edge index]`
/// (non-negative; `W::default()` must be zero).
///
/// # Panics
/// Panics if `source >= g.len()`, or a weight is missing for some edge index.
pub fn dijkstra<W>(g: &Csr, weight: &[W], source: usize) -> ShortestPaths<W>
where
    W: Copy + Ord + Default + Add<Output = W>,
{
    let n = g.len();
    assert!(source < n, "source vertex {source} out of range for {n} vertices");
    let mut dist = vec![W::default(); n];
    let mut reached = vec![false; n];
    let mut done = vec![false; n];
    let mut parent = vec![ShortestPaths::<W>::NONE; n];
    let mut parent_edge = vec![ShortestPaths::<W>::NONE; n];
    let mut heap: BinaryHeap<Reverse<(W, usize)>> = BinaryHeap::new();
    reached[source] = true;
    heap.push(Reverse((W::default(), source)));
    while let Some(Reverse((d, u))) = heap.pop() {
        if done[u] {
            continue;
        }
        done[u] = true;
        for a in g.arcs(u) {
            let v = g.to[a];
            if done[v] {
                continue;
            }
            let nd = d + weight[g.id[a]];
            if !reached[v] || nd < dist[v] {
                reached[v] = true;
                dist[v] = nd;
                parent[v] = u;
                parent_edge[v] = g.id[a];
                heap.push(Reverse((nd, v)));
            }
        }
    }
    ShortestPaths { dist, reached, parent, parent_edge }
}
