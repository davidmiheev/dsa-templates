//! Eulerian trails (a walk using every edge exactly once) by an iterative Hierholzer search,
//! `O(n + m)`, for directed and undirected multigraphs (loops and parallel edges allowed).
//!
//! A [`Trail`] lists the visited vertices `v_0 .. v_m` and the used edge indices `e_0 .. e_{m-1}`
//! (`e_i` joins `v_i` and `v_{i+1}`). When **no** Eulerian trail exists, `vertices` is empty
//! (a real trail always has at least one vertex, so this is unambiguous). A graph without edges
//! has the trivial trail consisting of vertex `0`.
//!
//! # Example
//!
//! ```
//! use eulerian_trail::{eulerian_trail_directed, eulerian_trail_undirected};
//!
//! // 0 -> 1 -> 2 -> 0 and 0 -> 3.
//! let t = eulerian_trail_directed(4, &[(0, 1), (1, 2), (2, 0), (0, 3)]);
//! assert_eq!(t.vertices, vec![0, 1, 2, 0, 3]);
//! assert_eq!(t.edges, vec![0, 1, 2, 3]);
//! // A star with three leaves has four odd vertices: no trail.
//! assert!(eulerian_trail_undirected(4, &[(0, 1), (0, 2), (0, 3)]).vertices.is_empty());
//! ```

use csr_graph::Csr;

/// An Eulerian trail. See the [module documentation](self).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Trail {
    /// Visited vertices, `edges.len() + 1` of them; empty if there is no Eulerian trail.
    pub vertices: Vec<usize>,
    /// Edge indices in walking order.
    pub edges: Vec<usize>,
}

impl Trail {
    fn none() -> Self {
        Trail { vertices: Vec::new(), edges: Vec::new() }
    }

    /// Whether an Eulerian trail was found.
    pub fn exists(&self) -> bool {
        !self.vertices.is_empty()
    }
}

/// An Eulerian trail of the directed multigraph on `0..n`, or the empty [`Trail`] if none
/// exists. If a circuit exists the trail starts and ends at the same vertex.
///
/// # Panics
/// Panics if `n == 0` or an endpoint is `>= n`.
pub fn eulerian_trail_directed(n: usize, edges: &[(usize, usize)]) -> Trail {
    assert!(n > 0, "the graph needs at least one vertex");
    let g = Csr::directed(n, edges);
    let mut balance = vec![0i64; n]; // out - in
    for &(u, v) in edges {
        balance[u] += 1;
        balance[v] -= 1;
    }
    let (mut start, mut ends) = (None, 0);
    for v in 0..n {
        match balance[v] {
            0 => {}
            1 if start.is_none() => start = Some(v),
            -1 if ends == 0 => ends = 1,
            _ => return Trail::none(),
        }
    }
    if start.is_some() != (ends == 1) {
        return Trail::none();
    }
    let start = start.or_else(|| (0..n).find(|&v| g.degree(v) > 0)).unwrap_or(0);
    // Hierholzer: walk until stuck, emit on backtrack, then reverse.
    let mut next_arc = g.start[..n].to_vec();
    let mut stack: Vec<(usize, usize)> = vec![(start, usize::MAX)];
    let (mut vs, mut es) = (Vec::with_capacity(edges.len() + 1), Vec::with_capacity(edges.len()));
    while let Some(&(v, e)) = stack.last() {
        if next_arc[v] < g.start[v + 1] {
            let a = next_arc[v];
            next_arc[v] += 1;
            stack.push((g.to[a], g.id[a]));
        } else {
            stack.pop();
            vs.push(v);
            if e != usize::MAX {
                es.push(e);
            }
        }
    }
    if es.len() != edges.len() {
        return Trail::none();
    }
    vs.reverse();
    es.reverse();
    Trail { vertices: vs, edges: es }
}

/// An Eulerian trail of the undirected multigraph on `0..n`, or the empty [`Trail`] if none
/// exists. If a circuit exists the trail starts and ends at the same vertex.
///
/// # Panics
/// Panics if `n == 0` or an endpoint is `>= n`.
pub fn eulerian_trail_undirected(n: usize, edges: &[(usize, usize)]) -> Trail {
    assert!(n > 0, "the graph needs at least one vertex");
    let g = Csr::undirected(n, edges);
    let odd: Vec<usize> = (0..n).filter(|&v| g.degree(v) % 2 == 1).collect();
    if odd.len() > 2 {
        return Trail::none();
    }
    let start = odd.first().copied().or_else(|| (0..n).find(|&v| g.degree(v) > 0)).unwrap_or(0);
    let mut next_arc = g.start[..n].to_vec();
    let mut used = vec![false; edges.len()];
    let mut stack: Vec<(usize, usize)> = vec![(start, usize::MAX)];
    let (mut vs, mut es) = (Vec::with_capacity(edges.len() + 1), Vec::with_capacity(edges.len()));
    while let Some(&(v, e)) = stack.last() {
        while next_arc[v] < g.start[v + 1] && used[g.id[next_arc[v]]] {
            next_arc[v] += 1;
        }
        if next_arc[v] < g.start[v + 1] {
            let a = next_arc[v];
            next_arc[v] += 1;
            used[g.id[a]] = true;
            stack.push((g.to[a], g.id[a]));
        } else {
            stack.pop();
            vs.push(v);
            if e != usize::MAX {
                es.push(e);
            }
        }
    }
    if es.len() != edges.len() {
        return Trail::none();
    }
    vs.reverse();
    es.reverse();
    Trail { vertices: vs, edges: es }
}
