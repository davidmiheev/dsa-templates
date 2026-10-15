//! Connected components of the **complement** of a graph without ever building it.
//!
//! BFS over the complement keeps the not-yet-visited vertices in a [`fast_set::FastSet`]. When a
//! vertex `v` is popped, its original neighbours are stamped and the unvisited set is scanned:
//! every unvisited vertex that is not stamped is a complement neighbour and is removed, every
//! stamped one is skipped (and charged to the edge that stamped it). Total `O((n + m) log_64 n)`.
//!
//! # Example
//!
//! ```
//! use complement_graph::complement_components;
//!
//! // Path 0-1-2-3: the complement has edges 0-2, 0-3, 1-3, so it is connected.
//! assert_eq!(complement_components(4, &[(0, 1), (1, 2), (2, 3)]).len(), 1);
//! // Complete graph K3: the complement has no edges.
//! assert_eq!(complement_components(3, &[(0, 1), (1, 2), (0, 2)]).len(), 3);
//! ```

use csr_graph::Csr;
use fast_set::{FastSet, NONE};

/// Components of the complement of the undirected graph `edges` on `0..n`, each as a list of
/// vertices (in BFS order, starting from its smallest vertex); components are ordered by their
/// smallest vertex. Parallel edges and self-loops are ignored.
pub fn complement_components(n: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    complement_components_csr(&Csr::undirected(n, edges))
}

/// Same as [`complement_components`] for a prebuilt undirected [`Csr`].
pub fn complement_components_csr(g: &Csr) -> Vec<Vec<usize>> {
    let n = g.len();
    let mut unvisited = FastSet::new(n.max(1));
    for v in 0..n {
        unvisited.insert(v);
    }
    let mut stamp = vec![NONE; n];
    let mut comps: Vec<Vec<usize>> = Vec::new();
    loop {
        let s = unvisited.next(0);
        if s == NONE {
            return comps;
        }
        unvisited.remove(s);
        let mut comp = vec![s];
        let mut head = 0;
        while head < comp.len() {
            let v = comp[head];
            head += 1;
            for &w in g.neighbors(v) {
                stamp[w] = v;
            }
            let mut u = unvisited.next(0);
            while u != NONE {
                if stamp[u] != v {
                    unvisited.remove(u);
                    comp.push(u);
                }
                u = unvisited.next(u + 1);
            }
        }
        comps.push(comp);
    }
}
