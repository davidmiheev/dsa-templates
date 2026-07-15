//! Depth-First Search (DFS) on a graph.
//!
//! DFS explores as deep as possible along each branch before
//! backtracking. It is the natural building block for topological
//! order, cycle detection on directed/undirected graphs, strongly
//! connected components, and tree traversals.
//!
//! This module exposes DFS as two helpers operating on an adjacency
//! list:
//!
//! * [`dfs_order`] — returns vertices in the order they are first
//!   visited (iterative, stack-based; safe for deep graphs).
//! * [`has_cycle_undirected`] — detects whether an undirected graph
//!   contains a cycle.

use std::collections::HashSet;

/// Return the DFS visitation order starting from `source`.
///
/// Iterative implementation using an explicit stack, so it is safe for
/// graphs deeper than the default Rust stack limit.
///
/// # Panics
///
/// Panics if `source >= graph.len()`.
pub fn dfs_order(graph: &[Vec<usize>], source: usize) -> Vec<usize> {
    assert!(source < graph.len(), "source vertex out of range");
    let mut visited = vec![false; graph.len()];
    let mut order = Vec::with_capacity(graph.len());
    let mut stack: Vec<usize> = Vec::new();
    stack.push(source);
    while let Some(v) = stack.pop() {
        if visited[v] {
            continue;
        }
        visited[v] = true;
        order.push(v);
        // push neighbours in reverse so the natural left-to-right order
        // is preserved in the resulting traversal
        for &u in graph[v].iter().rev() {
            if !visited[u] {
                stack.push(u);
            }
        }
    }
    order
}

/// Detect whether an undirected graph contains a cycle.
///
/// Pass the adjacency list and edges must already be reciprocal (i.e.,
/// for any edge `u -- v` you must include both `v` in `graph[u]` and
/// `u` in `graph[v]`).
pub fn has_cycle_undirected(graph: &[Vec<usize>]) -> bool {
    let mut visited: HashSet<usize> = HashSet::new();
    for start in 0..graph.len() {
        if visited.contains(&start) {
            continue;
        }
        // iterative DFS on the connected component of `start`
        let mut stack: Vec<(usize, usize)> = Vec::new(); // (node, parent)
        stack.push((start, usize::MAX));
        while let Some((v, parent)) = stack.pop() {
            if visited.contains(&v) {
                // reached an already-visited node via a different path:
                // this is a back edge, hence a cycle
                return true;
            }
            visited.insert(v);
            for &u in &graph[v] {
                if u == parent {
                    continue;
                }
                stack.push((u, v));
            }
        }
    }
    false
}
