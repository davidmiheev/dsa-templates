//! Dijkstra's single-source shortest path on a weighted graph.
//!
//! Given a directed or undirected graph with non-negative edge weights,
//! Dijkstra's algorithm computes the shortest distance from a source
//! vertex to every other reachable vertex in `O((V + E) log V)` time
//! using a binary heap.
//!
//! This module exposes a single function [`dijkstra`] that takes the
//! adjacency list as `&[Vec<(usize, u64)>]` (each edge is `(neighbor,
//! weight)`) and returns a `HashMap<usize, u64>` of minimum distances.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

/// Compute shortest distances from `source` to every reachable vertex.
///
/// `graph[u]` is the list of edges leaving `u`: pairs `(v, w)` meaning
/// there is an edge to `v` with weight `w`. Returns a map from vertex
/// id to minimum distance. Unreachable vertices are absent.
///
/// # Panics
///
/// Panics if `source >= graph.len()`.
pub fn dijkstra(graph: &[Vec<(usize, u64)>], source: usize) -> HashMap<usize, u64> {
    assert!(source < graph.len(), "source vertex out of range");
    let mut dist: HashMap<usize, u64> = HashMap::new();
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new();
    dist.insert(source, 0);
    heap.push(Reverse((0, source)));
    while let Some(Reverse((d, u))) = heap.pop() {
        if d > dist[&u] {
            continue; // stale entry
        }
        for &(v, w) in &graph[u] {
            let nd = d + w;
            let is_better = match dist.get(&v) {
                None => true,
                Some(cur) => nd < *cur,
            };
            if is_better {
                dist.insert(v, nd);
                heap.push(Reverse((nd, v)));
            }
        }
    }
    dist
}
