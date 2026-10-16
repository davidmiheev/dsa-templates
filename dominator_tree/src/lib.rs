//! Dominator tree of a directed graph rooted at `root`, by the Lengauer–Tarjan algorithm
//! (semidominators + path-compressed `eval`), in `O((n + m) log n)`. Everything is iterative, so
//! deep graphs are fine.
//!
//! Vertex `d` *dominates* `v` if every path `root -> v` passes through `d`. The immediate
//! dominator `idom[v]` is the closest strict dominator; `idom[root] == root`, and unreachable
//! vertices have `idom == NONE`.
//!
//! # Example
//!
//! ```
//! use dominator_tree::{DominatorTree, NONE};
//!
//! // 0 -> 1 -> 3, 0 -> 2 -> 3, 3 -> 4; vertex 5 is unreachable.
//! let d = DominatorTree::from_edges(6, &[(0, 1), (0, 2), (1, 3), (2, 3), (3, 4)], 0);
//! assert_eq!(d.idom, vec![0, 0, 0, 0, 3, NONE]);
//! assert!(d.dominates(0, 4) && d.dominates(3, 4) && !d.dominates(1, 3));
//! ```

use csr_graph::Csr;

/// Marks "unreachable" in [`DominatorTree::idom`].
pub const NONE: usize = usize::MAX;

/// See the [module documentation](self).
#[derive(Debug, Clone)]
pub struct DominatorTree {
    pub root: usize,
    /// Immediate dominator of every vertex (`idom[root] == root`, `NONE` if unreachable).
    pub idom: Vec<usize>,
}

impl DominatorTree {
    /// Build from an edge list.
    pub fn from_edges(n: usize, edges: &[(usize, usize)], root: usize) -> Self {
        Self::new(&Csr::directed(n, edges), root)
    }

    /// Build from a directed CSR graph.
    pub fn new(g: &Csr, root: usize) -> Self {
        let n = g.len();
        assert!(root < n, "root {root} out of range for {n} vertices");
        let rg = g.reversed();
        // Iterative DFS numbering.
        let mut num = vec![NONE; n];
        let mut order = Vec::with_capacity(n);
        let mut parent: Vec<usize> = Vec::with_capacity(n); // by dfs number
        let mut next_arc = g.start[..n].to_vec();
        let mut stack = vec![root];
        num[root] = 0;
        order.push(root);
        parent.push(0);
        while let Some(&v) = stack.last() {
            if next_arc[v] < g.start[v + 1] {
                let w = g.to[next_arc[v]];
                next_arc[v] += 1;
                if num[w] == NONE {
                    num[w] = order.len();
                    parent.push(num[v]);
                    order.push(w);
                    stack.push(w);
                }
            } else {
                stack.pop();
            }
        }
        let cnt = order.len();
        // Lengauer–Tarjan, everything indexed by dfs number.
        let mut semi: Vec<usize> = (0..cnt).collect();
        let mut label: Vec<usize> = (0..cnt).collect();
        let mut ancestor = vec![NONE; cnt];
        let mut idom = vec![0usize; cnt];
        let mut head = vec![NONE; cnt]; // bucket lists
        let mut link = vec![NONE; cnt];
        let mut chain: Vec<usize> = Vec::new();
        for i in (1..cnt).rev() {
            let v = order[i];
            for &u in rg.neighbors(v) {
                if num[u] == NONE {
                    continue;
                }
                let j = num[u];
                let e = eval(j, &mut ancestor, &mut label, &semi, &mut chain);
                if semi[e] < semi[i] {
                    semi[i] = semi[e];
                }
            }
            link[i] = head[semi[i]];
            head[semi[i]] = i;
            let p = parent[i];
            ancestor[i] = p;
            let mut w = head[p];
            head[p] = NONE;
            while w != NONE {
                let nxt = link[w];
                let u = eval(w, &mut ancestor, &mut label, &semi, &mut chain);
                idom[w] = if semi[u] < semi[w] { u } else { p };
                w = nxt;
            }
        }
        for i in 1..cnt {
            if idom[i] != semi[i] {
                idom[i] = idom[idom[i]];
            }
        }
        let mut res = vec![NONE; n];
        res[root] = root;
        for i in 1..cnt {
            res[order[i]] = order[idom[i]];
        }
        Self { root, idom: res }
    }

    /// Whether `v` is reachable from the root.
    pub fn reachable(&self, v: usize) -> bool {
        self.idom[v] != NONE
    }

    /// Whether `d` dominates `v` (every vertex dominates itself). Both must be reachable. `O(depth)`.
    pub fn dominates(&self, d: usize, mut v: usize) -> bool {
        assert!(self.reachable(d) && self.reachable(v), "dominates() needs reachable vertices");
        loop {
            if v == d {
                return true;
            }
            if v == self.root {
                return false;
            }
            v = self.idom[v];
        }
    }
}

/// Path-compressed `eval` of the Lengauer–Tarjan forest (iterative).
fn eval(v: usize, ancestor: &mut [usize], label: &mut [usize], semi: &[usize], chain: &mut Vec<usize>) -> usize {
    if ancestor[v] == NONE {
        return v;
    }
    chain.clear();
    let mut x = v;
    while ancestor[ancestor[x]] != NONE {
        chain.push(x);
        x = ancestor[x];
    }
    while let Some(x) = chain.pop() {
        let a = ancestor[x];
        if semi[label[a]] < semi[label[x]] {
            label[x] = label[a];
        }
        ancestor[x] = ancestor[a];
    }
    label[v]
}
