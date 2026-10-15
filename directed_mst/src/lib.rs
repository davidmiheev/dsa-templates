//! Minimum-cost arborescence (directed MST) rooted at `root`, by Chu–Liu/Edmonds in
//! `O(m log m)`: every vertex keeps a mergeable (skew) heap of its incoming edges with a lazy
//! "subtract" tag, cycles are contracted with a [`union_find::RollbackUnionFind`], and the
//! contractions are rolled back at the end to recover which edge enters each vertex. The heaps
//! are merged without recursion, so adversarial inputs cannot overflow the stack.
//!
//! # Example
//!
//! ```
//! use directed_mst::{directed_mst, NONE};
//!
//! // 0 -> 1 (5), 0 -> 2 (1), 2 -> 1 (1), 1 -> 2 (0)
//! let t = directed_mst(3, 0, &[(0, 1, 5i64), (0, 2, 1), (2, 1, 1), (1, 2, 0)]);
//! assert_eq!(t.cost, 2);
//! assert_eq!(t.parent, vec![0, 2, 0]);
//! assert_eq!(t.edge, vec![NONE, 2, 1]);
//! ```

use std::ops::{Add, Sub};
use union_find::RollbackUnionFind;

/// "No incoming edge" (the root) in [`Arborescence::edge`].
pub const NONE: usize = usize::MAX;
const NIL: usize = usize::MAX;

/// The result of [`directed_mst`].
#[derive(Debug, Clone)]
pub struct Arborescence<W> {
    /// Total weight.
    pub cost: W,
    /// Parent of every vertex; `parent[root] == root`.
    pub parent: Vec<usize>,
    /// Index (into the input edge list) of the edge entering every vertex; `NONE` for the root.
    pub edge: Vec<usize>,
}

struct Node<W> {
    w: W,
    a: usize, // tail
    b: usize, // head
    id: usize,
    delta: W,
    l: usize,
    r: usize,
}

struct Heaps<W> {
    nodes: Vec<Node<W>>,
}

impl<W: Copy + Ord + Add<Output = W> + Sub<Output = W>> Heaps<W> {
    fn prop(&mut self, x: usize) {
        let d = self.nodes[x].delta;
        let (l, r) = (self.nodes[x].l, self.nodes[x].r);
        self.nodes[x].w = self.nodes[x].w + d;
        if l != NIL {
            self.nodes[l].delta = self.nodes[l].delta + d;
        }
        if r != NIL {
            self.nodes[r].delta = self.nodes[r].delta + d;
        }
        self.nodes[x].delta = d - d;
    }

    /// Iterative skew-heap merge (min at the root).
    fn merge(&mut self, mut a: usize, mut b: usize) -> usize {
        let (mut root, mut hook) = (NIL, NIL);
        loop {
            let (x, last) = if a == NIL {
                (b, true)
            } else if b == NIL {
                (a, true)
            } else {
                self.prop(a);
                self.prop(b);
                if self.nodes[a].w > self.nodes[b].w {
                    std::mem::swap(&mut a, &mut b);
                }
                let (l, r) = (self.nodes[a].l, self.nodes[a].r);
                self.nodes[a].l = NIL; // filled by the next round
                self.nodes[a].r = l;
                let x = a;
                a = r;
                (x, false)
            };
            if hook == NIL {
                root = x;
            } else {
                self.nodes[hook].l = x;
            }
            if last {
                return root;
            }
            hook = x;
        }
    }

    /// Remove the root, returning the new root.
    fn pop(&mut self, x: usize) -> usize {
        self.prop(x);
        let (l, r) = (self.nodes[x].l, self.nodes[x].r);
        self.merge(l, r)
    }
}

/// Minimum arborescence of `edges = (tail, head, weight)` on `0..n` rooted at `root`.
/// Weights may be any ordered additive type with subtraction (integers; use a wide one).
///
/// # Panics
/// Panics if some vertex is unreachable from `root`, or an endpoint is out of range.
pub fn directed_mst<W>(n: usize, root: usize, edges: &[(usize, usize, W)]) -> Arborescence<W>
where
    W: Copy + Ord + Default + Add<Output = W> + Sub<Output = W>,
{
    assert!(root < n, "root {root} out of range for {n} vertices");
    let zero = W::default();
    let mut h = Heaps { nodes: Vec::with_capacity(edges.len()) };
    let mut heap = vec![NIL; n];
    for (id, &(a, b, w)) in edges.iter().enumerate() {
        assert!(a < n && b < n, "edge {id} has an endpoint outside 0..{n}");
        h.nodes.push(Node { w, a, b, id, delta: zero, l: NIL, r: NIL });
        let top = heap[b];
        heap[b] = h.merge(top, id);
    }
    let mut uf = RollbackUnionFind::new(n);
    let mut cost = zero;
    let mut seen = vec![NIL; n];
    let mut path = vec![0usize; n];
    // Chosen edges along the current walk: (tail, head, id).
    let mut chosen: Vec<(usize, usize, usize)> = vec![(0, 0, NIL); n];
    let mut inc = vec![(NIL, NIL, NIL); n];
    // (contracted representative, union-find snapshot before it, chosen edges of the cycle)
    let mut cycles: Vec<(usize, usize, Vec<(usize, usize, usize)>)> = Vec::new();
    seen[root] = root;
    for s in 0..n {
        let (mut u, mut qi) = (s, 0usize);
        while seen[u] == NIL {
            assert!(heap[u] != NIL, "vertex {u} is not reachable from the root");
            let top = heap[u];
            h.prop(top);
            let (w, a, b, id) = (h.nodes[top].w, h.nodes[top].a, h.nodes[top].b, h.nodes[top].id);
            h.nodes[top].delta = zero - w;
            heap[u] = h.pop(top);
            chosen[qi] = (a, b, id);
            path[qi] = u;
            qi += 1;
            seen[u] = s;
            cost = cost + w;
            u = uf.find(a);
            if seen[u] == s {
                // Found a cycle: contract it.
                let mut cyc = NIL;
                let end = qi;
                let snap = uf.snapshot();
                loop {
                    qi -= 1;
                    let w = path[qi];
                    cyc = h.merge(cyc, heap[w]);
                    if !uf.union(u, w) {
                        break;
                    }
                }
                u = uf.find(u);
                heap[u] = cyc;
                seen[u] = NIL;
                cycles.push((u, snap, chosen[qi..end].to_vec()));
            }
        }
        for &(a, b, id) in &chosen[..qi] {
            inc[uf.find(b)] = (a, b, id);
        }
    }
    // Expand the contractions in reverse order.
    while let Some((u, snap, comp)) = cycles.pop() {
        uf.rollback(snap);
        let entering = inc[u];
        for &(a, b, id) in &comp {
            inc[uf.find(b)] = (a, b, id);
        }
        inc[uf.find(entering.1)] = entering;
    }
    let mut parent = vec![root; n];
    let mut edge = vec![NONE; n];
    for v in 0..n {
        if v != root {
            parent[v] = inc[v].0;
            edge[v] = inc[v].2;
        }
    }
    Arborescence { cost, parent, edge }
}
