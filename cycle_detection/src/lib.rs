//! Find a simple cycle in a directed or an undirected multigraph, `O(n + m)`, iteratively.
//!
//! Both functions return the cycle as a list of **edge indices** in walking order, and the empty
//! list when the graph is acyclic (a cycle has at least one edge, so the empty list is
//! unambiguous). Use [`cycle_vertices_directed`] / [`cycle_vertices_undirected`] to turn the
//! edges into the vertex sequence.
//!
//! # Example
//!
//! ```
//! use cycle_detection::{find_cycle_directed, find_cycle_undirected};
//!
//! let edges = [(0, 1), (1, 2), (2, 0), (2, 3)];
//! assert_eq!(find_cycle_directed(4, &edges).len(), 3);
//! assert!(find_cycle_directed(4, &edges[..2]).is_empty());
//! // Undirected: the triangle again; parallel edges also form a cycle of length 2.
//! assert_eq!(find_cycle_undirected(4, &edges).len(), 3);
//! assert_eq!(find_cycle_undirected(2, &[(0, 1), (1, 0)]).len(), 2);
//! assert_eq!(find_cycle_undirected(1, &[(0, 0)]), vec![0]); // a self-loop
//! ```

use csr_graph::Csr;

const NONE: usize = usize::MAX;

/// A directed cycle (edge indices in order `e_0, e_1, ...` with `head(e_i) = tail(e_{i+1})`),
/// or an empty list if the digraph is acyclic. Self-loops count as cycles of length 1.
pub fn find_cycle_directed(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let g = Csr::directed(n, edges);
    // 0 = unseen, 1 = on the DFS stack, 2 = finished.
    let mut state = vec![0u8; n];
    let mut next_arc = g.start[..n].to_vec();
    let mut enter_edge = vec![NONE; n]; // edge used to enter the vertex
    let mut stack: Vec<usize> = Vec::new();
    for root in 0..n {
        if state[root] != 0 {
            continue;
        }
        state[root] = 1;
        stack.push(root);
        while let Some(&v) = stack.last() {
            if next_arc[v] < g.start[v + 1] {
                let a = next_arc[v];
                next_arc[v] += 1;
                let w = g.to[a];
                match state[w] {
                    0 => {
                        state[w] = 1;
                        enter_edge[w] = g.id[a];
                        stack.push(w);
                    }
                    1 => {
                        // Back edge to a vertex on the stack: walk the stack down to `w`.
                        let mut cycle = vec![g.id[a]];
                        let mut x = v;
                        while x != w {
                            cycle.push(enter_edge[x]);
                            x = edges[enter_edge[x]].0;
                        }
                        cycle.reverse();
                        return cycle;
                    }
                    _ => {}
                }
            } else {
                state[v] = 2;
                stack.pop();
            }
        }
    }
    Vec::new()
}

/// The vertex sequence `v_0, v_1, ...` of a directed cycle given by edge indices (`v_i` is the
/// tail of `cycle[i]`).
pub fn cycle_vertices_directed(edges: &[(usize, usize)], cycle: &[usize]) -> Vec<usize> {
    cycle.iter().map(|&e| edges[e].0).collect()
}

/// A simple cycle of an undirected multigraph (distinct vertices and distinct edges; two
/// parallel edges or one self-loop qualify), as edge indices in walking order, or an empty list
/// if the graph is a forest.
pub fn find_cycle_undirected(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let g = Csr::undirected(n, edges);
    let mut seen = vec![false; n];
    let mut on_stack = vec![false; n];
    let mut next_arc = g.start[..n].to_vec();
    let mut enter_edge = vec![NONE; n];
    let mut stack: Vec<usize> = Vec::new();
    for root in 0..n {
        if seen[root] {
            continue;
        }
        seen[root] = true;
        on_stack[root] = true;
        stack.push(root);
        while let Some(&v) = stack.last() {
            if next_arc[v] < g.start[v + 1] {
                let a = next_arc[v];
                next_arc[v] += 1;
                let (w, e) = (g.to[a], g.id[a]);
                if e == enter_edge[v] {
                    continue;
                }
                if !seen[w] {
                    seen[w] = true;
                    on_stack[w] = true;
                    enter_edge[w] = e;
                    stack.push(w);
                } else if on_stack[w] {
                    // Back edge to an ancestor (or a self-loop): close the cycle.
                    let mut cycle = vec![e];
                    let mut x = v;
                    while x != w {
                        cycle.push(enter_edge[x]);
                        let (a, b) = edges[enter_edge[x]];
                        x = if a == x { b } else { a };
                    }
                    cycle.reverse();
                    return cycle;
                }
            } else {
                on_stack[v] = false;
                stack.pop();
            }
        }
    }
    Vec::new()
}

/// The vertex sequence of an undirected cycle given by edge indices, starting at the vertex
/// where the walk begins (`v_i` and `v_{i+1}` are the ends of `cycle[i]`, `v_L = v_0`).
///
/// # Panics
/// Panics if `cycle` is not a closed walk.
pub fn cycle_vertices_undirected(edges: &[(usize, usize)], cycle: &[usize]) -> Vec<usize> {
    assert!(!cycle.is_empty(), "an empty cycle has no vertices");
    // Orient: the first vertex is the end of cycle[0] not shared with cycle[1] (for length 1
    // and 2 handle directly).
    let (a, b) = edges[cycle[0]];
    let start = if cycle.len() == 1 {
        a
    } else {
        let (c, d) = edges[cycle[1]];
        if b == c || b == d {
            a
        } else {
            b
        }
    };
    let mut vertices = vec![start];
    let mut v = start;
    for &e in cycle {
        let (x, y) = edges[e];
        assert!(v == x || v == y, "cycle edges are not consecutive");
        v = if v == x { y } else { x };
        vertices.push(v);
    }
    assert_eq!(v, start, "the walk is not closed");
    vertices.pop();
    vertices
}
