//! The `k` shortest **walks** (vertices and edges may repeat; equal lengths count separately)
//! from `s` to `t` in a directed graph with non-negative weights: Eppstein's algorithm,
//! `O(m log m + k log k)` after one Dijkstra run.
//!
//! Let `d(v)` be the distance from `v` to `t` and fix a shortest-path tree towards `t`. Every walk
//! is the tree walk from `s` with a sequence of *sidetracks* (non-tree edges `e = (u, v)` of
//! extra cost `delta(e) = w(e) + d(v) - d(u) >= 0`) spliced in. All sidetracks available after
//! arriving at a vertex `u` are the non-tree out-edges of every vertex on the tree path from `u`
//! to `t`; they are kept in one **persistent leftist heap** per vertex (`H(u)` = the out-edges of
//! `u` merged with `H(parent(u))`, sharing structure). Walks are then generated best-first: a
//! walk is a heap node, and its successors are the node's two heap children (swap the last
//! sidetrack for a slightly worse one) and the root of `H(head)` (append another sidetrack).
//!
//! # Example
//!
//! ```
//! use k_shortest_walks::k_shortest_walks;
//!
//! // 0 -> 1 (1), 1 -> 2 (1), 0 -> 2 (3), and a loop 1 -> 1 (1).
//! let lengths = k_shortest_walks(3, &[(0, 1, 1u64), (1, 2, 1), (0, 2, 3), (1, 1, 1)], 0, 2, 5);
//! assert_eq!(lengths, vec![2, 3, 3, 4, 5]);
//! ```

use csr_graph::Csr;
use dijkstra::dijkstra;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ops::{Add, Sub};

const NIL: u32 = u32::MAX;

#[derive(Clone, Copy)]
struct Node<W> {
    key: W,       // sidetrack cost
    head: u32,    // where the sidetrack edge arrives
    left: u32,
    right: u32,
    rank: u32,    // leftist "null path length"
}

struct Heaps<W> {
    nodes: Vec<Node<W>>,
}

impl<W: Copy + Ord> Heaps<W> {
    fn rank(&self, x: u32) -> u32 {
        if x == NIL {
            0
        } else {
            self.nodes[x as usize].rank
        }
    }

    /// Persistent merge: never modifies existing nodes. Iterative (the right spine is short, but
    /// the code stays recursion-free anyway).
    fn merge(&mut self, a: u32, b: u32) -> u32 {
        if a == NIL {
            return b;
        }
        if b == NIL {
            return a;
        }
        // Walk down the right spines, copying the nodes on the way, then fix ranks upwards.
        let mut spine: Vec<u32> = Vec::new();
        let (mut x, mut y) = (a, b);
        let result_tail;
        loop {
            if x == NIL {
                result_tail = y;
                break;
            }
            if y == NIL {
                result_tail = x;
                break;
            }
            if self.nodes[y as usize].key < self.nodes[x as usize].key {
                std::mem::swap(&mut x, &mut y);
            }
            let copy = self.nodes.len() as u32;
            let n = self.nodes[x as usize];
            self.nodes.push(n);
            spine.push(copy);
            x = n.right;
        }
        // Re-link: spine[i].right = spine[i+1] (or tail), then fix leftist property bottom-up.
        let mut below = result_tail;
        for &c in spine.iter().rev() {
            self.nodes[c as usize].right = below;
            let (l, r) = (self.nodes[c as usize].left, below);
            if self.rank(l) < self.rank(r) {
                self.nodes[c as usize].left = r;
                self.nodes[c as usize].right = l;
            }
            let rr = self.nodes[c as usize].right;
            self.nodes[c as usize].rank = self.rank(rr) + 1;
            below = c;
        }
        below
    }

    fn singleton(&mut self, key: W, head: u32) -> u32 {
        self.nodes.push(Node { key, head, left: NIL, right: NIL, rank: 1 });
        (self.nodes.len() - 1) as u32
    }
}

/// Lengths of the (at most) `k` shortest `s -> t` walks in increasing order; fewer than `k`
/// entries mean the graph has no more walks. Weights must be non-negative.
///
/// # Panics
/// Panics if `s`, `t` or an endpoint is out of range.
pub fn k_shortest_walks<W>(n: usize, edges: &[(usize, usize, W)], s: usize, t: usize, k: usize) -> Vec<W>
where
    W: Copy + Ord + Default + Add<Output = W> + Sub<Output = W>,
{
    assert!(s < n && t < n, "terminals out of range");
    let mut result = Vec::with_capacity(k.min(1 << 20));
    if k == 0 {
        return result;
    }
    let rev_pairs: Vec<(usize, usize)> = edges.iter().map(|&(u, v, _)| (v, u)).collect();
    let weight: Vec<W> = edges.iter().map(|&(_, _, w)| w).collect();
    let rg = Csr::directed(n, &rev_pairs);
    let sp = dijkstra(&rg, &weight, t); // distances *to* t; parent_edge = tree edge leaving v towards t
    if !sp.reachable(s) {
        return result;
    }
    let d = &sp.dist;
    // Tree children lists, then a breadth-first order from `t` so that parents come first.
    let mut tree_children: Vec<(usize, usize)> = Vec::new();
    for v in 0..n {
        if v != t && sp.reachable(v) {
            tree_children.push((edges[sp.parent_edge[v]].1, v));
        }
    }
    let tree = Csr::directed(n, &tree_children);
    let mut order = vec![t];
    let mut head = 0;
    while head < order.len() {
        let p = order[head];
        head += 1;
        order.extend_from_slice(tree.neighbors(p));
    }
    // Out-edge lists (original direction).
    let g = Csr::directed(n, &edges.iter().map(|&(u, v, _)| (u, v)).collect::<Vec<_>>());
    let mut heaps = Heaps { nodes: Vec::new() };
    let mut root = vec![NIL; n]; // persistent heap H(v): sidetracks along the tree path from v to t
    for &v in &order {
        // the non-tree out-edges of v, as a heap...
        let mut local = NIL;
        for a in g.arcs(v) {
            let (w, e) = (g.to[a], g.id[a]);
            if !sp.reachable(w) || (v != t && e == sp.parent_edge[v]) {
                continue;
            }
            let delta = weight[e] + d[w] - d[v];
            let single = heaps.singleton(delta, w as u32);
            local = heaps.merge(local, single);
        }
        // ... merged persistently with the heap of the tree parent
        root[v] = if v == t { local } else { heaps.merge(local, root[edges[sp.parent_edge[v]].1]) };
    }
    // Best-first enumeration.
    result.push(d[s]);
    if result.len() == k {
        return result;
    }
    let mut pq: BinaryHeap<Reverse<(W, u32)>> = BinaryHeap::new();
    if root[s] != NIL {
        pq.push(Reverse((d[s] + heaps.nodes[root[s] as usize].key, root[s])));
    }
    while let Some(Reverse((cost, x))) = pq.pop() {
        result.push(cost);
        if result.len() == k {
            break;
        }
        let node = heaps.nodes[x as usize];
        for child in [node.left, node.right] {
            if child != NIL {
                pq.push(Reverse((cost - node.key + heaps.nodes[child as usize].key, child)));
            }
        }
        let next = root[node.head as usize];
        if next != NIL {
            pq.push(Reverse((cost + heaps.nodes[next as usize].key, next)));
        }
    }
    result
}
