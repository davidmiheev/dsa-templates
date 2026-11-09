//! Lowlink values of an undirected multigraph, and what they give in `O(n + m)`:
//!
//! * [`bridges`](Lowlink::bridges) and [`articulation_points`](Lowlink::articulation_points),
//! * [`two_edge_connected_components`](Lowlink::two_edge_connected_components),
//! * [`biconnected_components`](Lowlink::biconnected_components) (vertex sets of the blocks).
//!
//! The graph may have parallel edges and self-loops (the DFS skips only the exact edge it came
//! from) and may be disconnected. The search is iterative, so deep graphs are fine.
//!
//! # Example
//!
//! ```
//! use csr_graph::Csr;
//! use lowlink::Lowlink;
//!
//! // A triangle 0-1-2 with a pendant vertex 3 hanging off 2.
//! let edges = [(0, 1), (1, 2), (2, 0), (2, 3)];
//! let ll = Lowlink::new(&Csr::undirected(4, &edges));
//! assert_eq!(ll.bridges(), vec![3]);
//! assert_eq!(ll.articulation_points(), vec![2]);
//! assert_eq!(ll.biconnected_components(), vec![vec![0, 1, 2], vec![2, 3]]);
//! ```

use csr_graph::Csr;

/// Lowlink data of a DFS forest. See the [module documentation](self).
#[derive(Debug, Clone)]
pub struct Lowlink {
    /// DFS discovery index of every vertex.
    pub order: Vec<usize>,
    /// `low[v]`: the smallest discovery index reachable from the subtree of `v` using tree edges
    /// down and at most one back edge.
    pub low: Vec<usize>,
    /// DFS parent ([`Lowlink::NONE`] for roots).
    pub parent: Vec<usize>,
    /// Index of the edge to the parent ([`Lowlink::NONE`] for roots).
    pub parent_edge: Vec<usize>,
    /// Vertices in discovery order (parents before their descendants).
    pub preorder: Vec<usize>,
}

impl Lowlink {
    /// The parent of a DFS root.
    pub const NONE: usize = usize::MAX;

    /// Run the search on the undirected graph `g` (use [`Csr::undirected`]).
    pub fn new(g: &Csr) -> Self {
        let n = g.len();
        let mut order = vec![Self::NONE; n];
        let mut low = vec![0usize; n];
        let mut parent = vec![Self::NONE; n];
        let mut parent_edge = vec![Self::NONE; n];
        let mut preorder = Vec::with_capacity(n);
        let mut next_arc = g.start[..n].to_vec();
        let mut call: Vec<usize> = Vec::new();
        for root in 0..n {
            if order[root] != Self::NONE {
                continue;
            }
            order[root] = preorder.len();
            low[root] = order[root];
            preorder.push(root);
            call.push(root);
            while let Some(&v) = call.last() {
                if next_arc[v] < g.start[v + 1] {
                    let a = next_arc[v];
                    next_arc[v] += 1;
                    let (w, e) = (g.to[a], g.id[a]);
                    if e == parent_edge[v] {
                        continue;
                    }
                    if order[w] == Self::NONE {
                        order[w] = preorder.len();
                        low[w] = order[w];
                        parent[w] = v;
                        parent_edge[w] = e;
                        preorder.push(w);
                        call.push(w);
                    } else {
                        low[v] = low[v].min(order[w]);
                    }
                } else {
                    call.pop();
                    if parent[v] != Self::NONE {
                        let p = parent[v];
                        low[p] = low[p].min(low[v]);
                    }
                }
            }
        }
        Lowlink { order, low, parent, parent_edge, preorder }
    }

    /// The edges whose removal disconnects their component, in discovery order of the lower
    /// endpoint.
    pub fn bridges(&self) -> Vec<usize> {
        self.preorder.iter().filter(|&&v| self.parent[v] != Self::NONE && self.low[v] > self.order[self.parent[v]]).map(|&v| self.parent_edge[v]).collect()
    }

    /// The vertices whose removal disconnects their component, ascending.
    pub fn articulation_points(&self) -> Vec<usize> {
        let n = self.order.len();
        let mut children = vec![0usize; n];
        let mut is_art = vec![false; n];
        for &v in &self.preorder {
            let p = self.parent[v];
            if p != Self::NONE {
                children[p] += 1;
                if self.parent[p] != Self::NONE && self.low[v] >= self.order[p] {
                    is_art[p] = true;
                }
            }
        }
        for v in 0..n {
            if self.parent[v] == Self::NONE && children[v] >= 2 {
                is_art[v] = true;
            }
        }
        (0..n).filter(|&v| is_art[v]).collect()
    }

    /// The 2-edge-connected components: their number and the component id of every vertex
    /// (ids follow the discovery order of the component roots).
    pub fn two_edge_connected_components(&self) -> (usize, Vec<usize>) {
        let n = self.order.len();
        let mut comp = vec![0usize; n];
        let mut count = 0;
        for &v in &self.preorder {
            let p = self.parent[v];
            if p == Self::NONE || self.low[v] > self.order[p] {
                comp[v] = count;
                count += 1;
            } else {
                comp[v] = comp[p];
            }
        }
        (count, comp)
    }

    /// The biconnected components (blocks) as vertex sets: every maximal 2-vertex-connected
    /// piece with at least one edge, plus one singleton per vertex without a non-loop edge.
    /// Each set is ascending; an articulation point appears in every block around it.
    pub fn biconnected_components(&self) -> Vec<Vec<usize>> {
        let n = self.order.len();
        // `block[v]` = the block of the tree edge into `v`; a block is headed by the parent of
        // the vertex that opened it.
        let mut block = vec![Self::NONE; n];
        let mut head: Vec<usize> = Vec::new();
        for &v in &self.preorder {
            let p = self.parent[v];
            if p == Self::NONE {
                continue;
            }
            if self.low[v] >= self.order[p] {
                block[v] = head.len();
                head.push(p);
            } else {
                block[v] = block[p];
            }
        }
        let mut sets: Vec<Vec<usize>> = head.iter().map(|&h| vec![h]).collect();
        for v in 0..n {
            if block[v] != Self::NONE {
                sets[block[v]].push(v);
            }
        }
        let mut has_edge = vec![false; n];
        for v in 0..n {
            if self.parent[v] != Self::NONE {
                has_edge[v] = true;
                has_edge[self.parent[v]] = true;
            }
        }
        for s in sets.iter_mut() {
            s.sort_unstable();
        }
        for (v, &covered) in has_edge.iter().enumerate() {
            if !covered {
                sets.push(vec![v]);
            }
        }
        sets
    }
}
