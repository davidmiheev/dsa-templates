//! Breadth-First Search (BFS) on a graph.
//!
//! BFS visits neighbours of the current node before going deeper. It is
//! typically used to compute shortest paths on unweighted graphs, to
//! find connected components, or to perform level-order traversal.
//!
//! This module exposes BFS as two helpers operating on an adjacency list:
//!
//! * [`bfs_distances`] — returns the shortest distance (in edges) from a
//!   source vertex to every reachable vertex. Unreachable vertices are
//!   absent from the result.
//! * [`bfs_order`] — returns the vertices in the order they are
//!   dequeued (i.e., the BFS visitation order).
//!
//! The graph is passed as `&[Vec<usize>]`, indexed by vertex. Vertex ids
//! are `0..n` where `n = graph.len()`.

use std::collections::{HashMap, VecDeque};

/// Compute the shortest distance (in number of edges) from `source` to
/// every vertex reachable from it.
///
/// Returns a map from vertex id to distance. Vertices that are not
/// reachable from `source` are absent. `source` itself maps to `0`.
///
/// # Panics
///
/// Panics if `source >= graph.len()`.
pub fn bfs_distances(graph: &[Vec<usize>], source: usize) -> HashMap<usize, usize> {
    assert!(source < graph.len(), "source vertex out of range");
    let mut dist: HashMap<usize, usize> = HashMap::new();
    let mut queue: VecDeque<usize> = VecDeque::new();
    dist.insert(source, 0);
    queue.push_back(source);
    while let Some(v) = queue.pop_front() {
        let d = dist[&v];
        for &u in &graph[v] {
            if !dist.contains_key(&u) {
                dist.insert(u, d + 1);
                queue.push_back(u);
            }
        }
    }
    dist
}

/// Return the visitation order of BFS starting from `source`.
///
/// The first element is always `source`. Each subsequent element is the
/// next vertex dequeued from the BFS queue.
///
/// # Panics
///
/// Panics if `source >= graph.len()`.
pub fn bfs_order(graph: &[Vec<usize>], source: usize) -> Vec<usize> {
    assert!(source < graph.len(), "source vertex out of range");
    let mut visited = vec![false; graph.len()];
    let mut order = Vec::with_capacity(graph.len());
    let mut queue: VecDeque<usize> = VecDeque::new();
    visited[source] = true;
    queue.push_back(source);
    while let Some(v) = queue.pop_front() {
        order.push(v);
        for &u in &graph[v] {
            if !visited[u] {
                visited[u] = true;
                queue.push_back(u);
            }
        }
    }
    order
}
