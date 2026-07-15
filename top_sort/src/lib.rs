//! Topological sort (Kahn's algorithm) on a directed acyclic graph.
//!
//! Given a DAG, a topological order is a linear ordering of the
//! vertices such that every edge `u -> v` has `u` before `v` in the
//! order. Kahn's algorithm produces such an order by repeatedly
//! removing vertices with in-degree zero.
//!
//! This module exposes [`topological_sort`], which returns `Some(order)`
//! if the graph is a DAG and `None` if the graph contains a cycle (in
//! which case no topological order exists).

use std::collections::VecDeque;

/// Compute a topological order of a directed graph.
///
/// `graph[u]` is the list of outgoing neighbours of `u`. Returns
/// `Some(Vec<usize>)` with a valid topological order if the graph is a
/// DAG, or `None` if the graph contains a cycle (cycle detection:
/// any vertex not visited by the algorithm at the end indicates one).
pub fn topological_sort(graph: &[Vec<usize>]) -> Option<Vec<usize>> {
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
    if order.len() == n {
        Some(order)
    } else {
        None
    }
}
