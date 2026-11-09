//! Strongly connected components of a directed graph, by an iterative Tarjan search (no
//! recursion, so million-vertex paths are fine), in `O(n + m)`.
//!
//! Components are numbered in **topological order of the condensation**: for every edge
//! `u -> v`, `comp[u] <= comp[v]`.
//!
//! # Example
//!
//! ```
//! use scc::Scc;
//!
//! // 0 <-> 1 -> 2 <-> 3
//! let s = Scc::from_edges(4, &[(0, 1), (1, 0), (1, 2), (2, 3), (3, 2)]);
//! assert_eq!(s.count, 2);
//! assert_eq!(s.comp, vec![0, 0, 1, 1]);
//! assert_eq!(s.groups(), vec![vec![0, 1], vec![2, 3]]);
//! ```

use csr_graph::Csr;

/// The result of the decomposition. See the [module documentation](self).
#[derive(Debug, Clone)]
pub struct Scc {
    /// Number of components.
    pub count: usize,
    /// Component id of every vertex, in topological order.
    pub comp: Vec<usize>,
}

impl Scc {
    /// Decompose the directed graph given by `edges` on `0..n`.
    pub fn from_edges(n: usize, edges: &[(usize, usize)]) -> Self {
        Self::new(&Csr::directed(n, edges))
    }

    /// Decompose a CSR graph (arcs are directed as stored).
    pub fn new(g: &Csr) -> Self {
        let n = g.len();
        const UNSEEN: usize = usize::MAX;
        let mut order = vec![UNSEEN; n]; // discovery index
        let mut low = vec![0usize; n];
        let mut comp = vec![UNSEEN; n];
        let mut next_arc = g.start[..n].to_vec();
        let mut stack: Vec<usize> = Vec::new(); // Tarjan stack
        let mut call: Vec<usize> = Vec::new(); // DFS stack
        let (mut time, mut count) = (0usize, 0usize);
        for root in 0..n {
            if order[root] != UNSEEN {
                continue;
            }
            call.push(root);
            order[root] = time;
            low[root] = time;
            time += 1;
            stack.push(root);
            while let Some(&v) = call.last() {
                if next_arc[v] < g.start[v + 1] {
                    let w = g.to[next_arc[v]];
                    next_arc[v] += 1;
                    if order[w] == UNSEEN {
                        order[w] = time;
                        low[w] = time;
                        time += 1;
                        stack.push(w);
                        call.push(w);
                    } else if comp[w] == UNSEEN {
                        low[v] = low[v].min(order[w]);
                    }
                } else {
                    call.pop();
                    if let Some(&p) = call.last() {
                        low[p] = low[p].min(low[v]);
                    }
                    if low[v] == order[v] {
                        loop {
                            let w = stack.pop().expect("v is on the stack");
                            comp[w] = count;
                            if w == v {
                                break;
                            }
                        }
                        count += 1;
                    }
                }
            }
        }
        // Tarjan numbers components in reverse topological order; flip.
        for c in comp.iter_mut() {
            *c = count - 1 - *c;
        }
        Scc { count, comp }
    }

    /// The vertices of every component, components in topological order, vertices ascending.
    pub fn groups(&self) -> Vec<Vec<usize>> {
        let mut groups = vec![Vec::new(); self.count];
        for (v, &c) in self.comp.iter().enumerate() {
            groups[c].push(v);
        }
        groups
    }
}
