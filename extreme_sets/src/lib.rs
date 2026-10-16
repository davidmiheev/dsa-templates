//! The laminar tree of **extreme sets** of a weighted undirected graph.
//!
//! Repeat `n - 1` times: on the current contracted graph, scan the vertices in the order "take the
//! unvisited vertex with the smallest remaining weighted degree" (its degree to the unvisited
//! vertices, updated as vertices are visited), then merge the last two vertices of the scan into
//! one new vertex. The merges form a binary tree with `2n - 1` nodes (leaves `0..n`, the merge
//! made in round `i` is node `n + i`, the root is `2n - 2`). Its nodes are the sets that matter for
//! cut problems: for every nonnegative vertex weighting `a`, the global minimum of
//! `cut(X) + a(X)` over all nonempty vertex sets `X` is attained at a node of the tree.
//! That turns "global minimum cut of the graph plus a star with changing leaf weights" into
//! "path-to-root add, global minimum" on the tree.
//!
//! Construction is `O(n (n + m) log n)`.
//!
//! # Example
//!
//! ```
//! use extreme_sets::extreme_set_tree;
//!
//! // A triangle 0-1-2 with weights 3, 3, 1 on a path 0-1-2... plus a heavy edge 0-1.
//! let t = extreme_set_tree(3, &[(0, 1, 5i64), (1, 2, 1), (0, 2, 1)]);
//! assert_eq!(t.nodes(), 5);
//! assert_eq!(t.cut[t.root()], 0);
//! // The leaves' cuts are the weighted degrees.
//! assert_eq!(&t.cut[..3], &[6, 6, 2]);
//! ```

use std::collections::BinaryHeap;
use std::cmp::Reverse;
use std::ops::{Add, Sub};

/// "No parent" (the root).
pub const NONE: usize = usize::MAX;

/// The tree of extreme sets. See the [module documentation](self).
#[derive(Debug, Clone)]
pub struct ExtremeSetTree<W> {
    /// Number of graph vertices (leaves).
    pub n: usize,
    /// Parent of every node ([`NONE`] for the root).
    pub parent: Vec<usize>,
    /// The two children of every internal node (empty for leaves).
    pub children: Vec<Vec<usize>>,
    /// Weight of the edges leaving the vertex set of every node.
    pub cut: Vec<W>,
}

impl<W> ExtremeSetTree<W> {
    /// Total number of nodes, `2n - 1`.
    pub fn nodes(&self) -> usize {
        2 * self.n - 1
    }

    /// The root (the whole vertex set).
    pub fn root(&self) -> usize {
        2 * self.n - 2
    }
}

/// Build the extreme-set tree of the undirected graph `edges = (u, v, weight)` on `0..n`
/// (`n >= 1`; parallel edges are summed, self-loops and zero weights ignored). Weights must be
/// non-negative.
pub fn extreme_set_tree<W>(n: usize, edges: &[(usize, usize, W)]) -> ExtremeSetTree<W>
where
    W: Copy + Ord + Default + Add<Output = W> + Sub<Output = W>,
{
    assert!(n >= 1, "need at least one vertex");
    let zero = W::default();
    let total = 2 * n - 1;
    let mut g: Vec<Vec<(usize, W)>> = vec![Vec::new(); total];
    for &(u, v, w) in edges {
        assert!(u < n && v < n, "edge ({u}, {v}) has an endpoint outside 0..{n}");
        if u != v && w != zero {
            g[u].push((v, w));
            g[v].push((u, w));
        }
    }
    let original = g.clone();
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); total];
    let mut used = vec![true; total];
    for i in 0..n - 1 {
        let live = n + i;
        let mut cost = vec![zero; live];
        let mut visited = vec![false; live];
        let mut heap: BinaryHeap<(Reverse<W>, usize)> = BinaryHeap::new();
        for j in 0..live {
            if used[j] {
                cost[j] = g[j].iter().fold(zero, |s, &(_, w)| s + w);
                heap.push((Reverse(cost[j]), j));
            }
        }
        let mut order = Vec::with_capacity(n - i);
        while order.len() < n - i {
            let (_, v) = heap.pop().expect("a live vertex remains");
            if visited[v] {
                continue;
            }
            visited[v] = true;
            order.push(v);
            for k in 0..g[v].len() {
                let (u, w) = g[v][k];
                if !visited[u] {
                    cost[u] = cost[u] - w;
                    heap.push((Reverse(cost[u]), u));
                }
            }
        }
        let (v0, v1) = (order[n - i - 1], order[n - i - 2]);
        g[v0].clear();
        g[v1].clear();
        used[v0] = false;
        used[v1] = false;
        let merged = n + i;
        let mut incoming = Vec::new();
        for j in 0..live {
            if !used[j] {
                continue;
            }
            for p in g[j].iter_mut() {
                if p.0 == v0 || p.0 == v1 {
                    p.0 = merged;
                    incoming.push((j, p.1));
                }
            }
        }
        g[merged] = incoming;
        children[merged] = vec![v0, v1];
    }
    let mut parent = vec![NONE; total];
    for v in n..total {
        for &c in &children[v] {
            parent[c] = v;
        }
    }
    // cut(node) = sum of w over edges with exactly one end in the node's leaf set: add w to every
    // node on the tree path between the two ends, excluding their lowest common ancestor.
    let mut depth = vec![0usize; total];
    for v in (0..total - 1).rev() {
        depth[v] = depth[parent[v]] + 1;
    }
    let mut cut = vec![zero; total];
    for u in 0..n {
        for &(v, w) in &original[u] {
            if u < v {
                let (mut a, mut b) = (u, v);
                while a != b {
                    if depth[a] < depth[b] {
                        std::mem::swap(&mut a, &mut b);
                    }
                    cut[a] = cut[a] + w;
                    a = parent[a];
                }
            }
        }
    }
    ExtremeSetTree { n, parent, children, cut }
}
