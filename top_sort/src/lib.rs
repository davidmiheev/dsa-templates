//! Topological sort (Kahn's algorithm) on a directed acyclic graph.
//!
//! Given a DAG, a topological order is a linear ordering of the
//! vertices such that every edge `u -> v` has `u` before `v` in the
//! order. Kahn's algorithm produces such an order by repeatedly
//! removing vertices with in-degree zero.
//!
//! This module exposes [`topological_sort`], which returns an order and
//! asserts that the graph is a DAG, and [`is_acyclic`] to test that
//! precondition without panicking.

use std::collections::VecDeque;

/// Compute a topological order of a directed graph.
///
/// `graph[u]` is the list of outgoing neighbours of `u`.
///
/// # Panics
/// Panics if the graph contains a cycle (use [`is_acyclic`] to check first).
pub fn topological_sort(graph: &[Vec<usize>]) -> Vec<usize> {
    let order = kahn(graph);
    assert_eq!(order.len(), graph.len(), "the graph contains a cycle");
    order
}

/// Whether the directed graph `graph` has no cycle.
pub fn is_acyclic(graph: &[Vec<usize>]) -> bool {
    kahn(graph).len() == graph.len()
}

/// Kahn's algorithm: the vertices it manages to remove, in order. Fewer than `n` of them
/// means a cycle.
fn kahn(graph: &[Vec<usize>]) -> Vec<usize> {
    let n = graph.len();
    let mut indeg: Vec<usize> = vec![0; n];
    for u in 0..n {
        for &v in &graph[u] {
            indeg[v] += 1;
        }
    }
    let mut queue: VecDeque<usize> = VecDeque::new();
    for v in 0..n {
        if indeg[v] == 0 {
            queue.push_back(v);
        }
    }

    let mut order = Vec::with_capacity(n);
    while let Some(v) = queue.pop_front() {
        order.push(v);
        for &u in &graph[v] {
            indeg[u] -= 1;
            if indeg[u] == 0 {
                queue.push_back(u);
            }
        }
    }
    order
}
