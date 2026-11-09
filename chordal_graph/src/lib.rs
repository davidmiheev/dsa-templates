//! Chordal graphs: a graph is chordal iff it has a *perfect elimination ordering* (PEO), i.e.
//! an order in which every vertex together with its later neighbours forms a clique.
//!
//! [`recognize`] runs maximum cardinality search (Tarjan–Yannakakis, bucket lists, `O(n + m)`),
//! takes the reverse of the visiting order as the candidate PEO and checks it in `O(n + m)`
//! (each vertex's later neighbours must be adjacent to its *follower*, the first of them in the
//! order). If the check fails the result is a certificate instead: an **induced cycle of length
//! at least 4**, found as `v` + a shortest path between two non-adjacent later neighbours of
//! `v` that avoids the rest of `N(v)`.
//!
//! # Example
//!
//! ```
//! use chordal_graph::{recognize, Chordality};
//!
//! // A triangle with a pendant vertex is chordal.
//! assert!(matches!(recognize(4, &[(0, 1), (1, 2), (2, 0), (2, 3)]), Chordality::Chordal { .. }));
//! // The 4-cycle is not: the certificate is the cycle itself.
//! match recognize(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]) {
//!     Chordality::Chordal { .. } => unreachable!(),
//!     Chordality::Hole { cycle } => assert_eq!(cycle.len(), 4),
//! }
//! ```

use csr_graph::Csr;
use std::collections::VecDeque;

const NONE: usize = usize::MAX;

/// Result of [`recognize`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Chordality {
    /// Chordal, with a perfect elimination ordering (every vertex exactly once).
    Chordal { peo: Vec<usize> },
    /// Not chordal, with an induced cycle (length >= 4) as the vertex sequence around it.
    Hole { cycle: Vec<usize> },
}

/// Maximum cardinality search order: the vertices in the order they are visited (each time a
/// vertex with the most already-visited neighbours is taken). Its reverse is a PEO iff the graph
/// is chordal.
pub fn maximum_cardinality_search(g: &Csr) -> Vec<usize> {
    let n = g.len();
    // Doubly linked bucket lists: `head[c]` holds the unvisited vertices with `c` visited neighbours.
    let mut head = vec![NONE; n + 1];
    let (mut next, mut prev) = (vec![NONE; n], vec![NONE; n]);
    let mut count = vec![0usize; n];
    let mut done = vec![false; n];
    for v in (0..n).rev() {
        next[v] = head[0];
        if head[0] != NONE {
            prev[head[0]] = v;
        }
        head[0] = v;
    }
    let mut order = Vec::with_capacity(n);
    let mut best = 0usize;
    for _ in 0..n {
        while head[best] == NONE {
            best -= 1;
        }
        let v = head[best];
        // unlink v
        head[best] = next[v];
        if next[v] != NONE {
            prev[next[v]] = NONE;
        }
        done[v] = true;
        order.push(v);
        for &w in g.neighbors(v) {
            if done[w] {
                continue;
            }
            let c = count[w];
            // unlink w from bucket c
            if prev[w] != NONE {
                next[prev[w]] = next[w];
            } else {
                head[c] = next[w];
            }
            if next[w] != NONE {
                prev[next[w]] = prev[w];
            }
            count[w] = c + 1;
            next[w] = head[c + 1];
            prev[w] = NONE;
            if head[c + 1] != NONE {
                prev[head[c + 1]] = w;
            }
            head[c + 1] = w;
        }
        best += 1; // a neighbour may have moved one bucket up
        if best > n {
            best = n;
        }
    }
    order
}

/// Decide whether the simple undirected graph `edges` on `0..n` is chordal. Parallel edges and
/// self-loops are ignored.
pub fn recognize(n: usize, edges: &[(usize, usize)]) -> Chordality {
    let g = Csr::undirected(n, edges);
    let visit = maximum_cardinality_search(&g);
    // peo[i] = visit[n-1-i]; position in the PEO:
    let mut pos = vec![0usize; n];
    for (i, &v) in visit.iter().rev().enumerate() {
        pos[v] = i;
    }
    let peo: Vec<usize> = visit.iter().rev().copied().collect();
    // For v, follower = later neighbour with the smallest position.
    let mut mark = vec![NONE; n];
    for &v in &peo {
        let mut follower = NONE;
        for &w in g.neighbors(v) {
            if w != v && pos[w] > pos[v] && (follower == NONE || pos[w] < pos[follower]) {
                follower = w;
            }
        }
        if follower == NONE {
            continue;
        }
        for &w in g.neighbors(follower) {
            mark[w] = follower;
        }
        for &w in g.neighbors(v) {
            if w != v && w != follower && pos[w] > pos[v] && mark[w] != follower {
                return Chordality::Hole { cycle: find_hole(&g, v, follower, w) };
            }
        }
    }
    Chordality::Chordal { peo }
}

/// `v` has the non-adjacent neighbours `f` and `w`: the shortest `f`-`w` path avoiding the rest
/// of `N(v)` closes with `v` into an induced cycle.
fn find_hole(g: &Csr, v: usize, f: usize, w: usize) -> Vec<usize> {
    let n = g.len();
    let mut banned = vec![false; n];
    for &x in g.neighbors(v) {
        banned[x] = true;
    }
    banned[f] = false;
    banned[w] = false;
    banned[v] = true;
    let mut from = vec![NONE; n];
    from[f] = f;
    let mut queue = VecDeque::from([f]);
    while let Some(x) = queue.pop_front() {
        if x == w {
            break;
        }
        for &y in g.neighbors(x) {
            if !banned[y] && from[y] == NONE {
                from[y] = x;
                queue.push_back(y);
            }
        }
    }
    assert!(from[w] != NONE, "no path between the violating neighbours: the order was not an MCS order");
    let mut cycle = vec![v];
    let mut x = w;
    while x != f {
        cycle.push(x);
        x = from[x];
    }
    cycle.push(f);
    cycle
}
