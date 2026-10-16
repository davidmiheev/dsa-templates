//! Minimum Steiner tree by the Dreyfus–Wagner dynamic programme, `O(3^k n + 2^k m log m)` for
//! `k` terminals, returning the chosen edges.
//!
//! `dp[T][v]` is the cheapest tree that spans the terminal set `T` and the vertex `v`. It is
//! obtained from `v` by either splitting `T` into two parts that both reach `v`, or by walking
//! one edge from a neighbour (the latter closed under shortest paths with a Dijkstra pass per
//! subset). Every cell remembers how it was derived so the tree can be read back.
//!
//! # Example
//!
//! ```
//! use steiner_tree::steiner_tree;
//!
//! // A star with centre 3: connecting the leaves 0, 1, 2 costs 3 through the centre.
//! let edges = [(0, 3, 1u64), (1, 3, 1), (2, 3, 1), (0, 1, 5), (1, 2, 5)];
//! let (cost, chosen) = steiner_tree(4, &edges, &[0, 1, 2]);
//! assert_eq!(cost, 3);
//! assert_eq!(chosen, vec![0, 1, 2]);
//! ```

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ops::Add;

#[derive(Clone, Copy)]
enum How {
    Terminal,
    Split(u32),
    Edge(usize), // arrived through this edge index
}

/// Minimum-weight subgraph of the undirected graph `edges = (u, v, weight)` on `0..n` that
/// connects all `terminals`. Returns `(total weight, sorted edge indices)`. Weights must be
/// non-negative.
///
/// # Panics
/// Panics if `terminals` is empty or some terminal is out of range, or if the terminals are not
/// all connected to each other.
pub fn steiner_tree<W>(n: usize, edges: &[(usize, usize, W)], terminals: &[usize]) -> (W, Vec<usize>)
where
    W: Copy + Ord + Default + Add<Output = W>,
{
    let k = terminals.len();
    assert!(k >= 1 && k <= 20, "need between 1 and 20 terminals, got {k}");
    let mut adj: Vec<Vec<(usize, W, usize)>> = vec![Vec::new(); n];
    for (i, &(u, v, w)) in edges.iter().enumerate() {
        assert!(u < n && v < n, "edge {i} has an endpoint outside 0..{n}");
        adj[u].push((v, w, i));
        adj[v].push((u, w, i));
    }
    let full = (1usize << k) - 1;
    let mut dp: Vec<Vec<Option<W>>> = vec![vec![None; n]; full + 1];
    let mut how: Vec<Vec<How>> = vec![vec![How::Terminal; n]; full + 1];
    for (i, &t) in terminals.iter().enumerate() {
        assert!(t < n, "terminal {t} out of range");
        dp[1 << i][t] = Some(W::default());
    }
    for mask in 1..=full {
        // merge two parts at the same vertex
        let mut sub = (mask - 1) & mask;
        while sub > 0 {
            let other = mask ^ sub;
            if sub < other {
                break; // each unordered split once
            }
            for v in 0..n {
                if let (Some(a), Some(b)) = (dp[sub][v], dp[other][v]) {
                    let c = a + b;
                    if dp[mask][v].map_or(true, |cur| c < cur) {
                        dp[mask][v] = Some(c);
                        how[mask][v] = How::Split(sub as u32);
                    }
                }
            }
            sub = (sub - 1) & mask;
        }
        // shortest-path closure
        let mut heap = BinaryHeap::new();
        for v in 0..n {
            if let Some(d) = dp[mask][v] {
                heap.push(Reverse((d, v)));
            }
        }
        while let Some(Reverse((d, v))) = heap.pop() {
            if dp[mask][v] != Some(d) {
                continue;
            }
            for &(w, c, id) in &adj[v] {
                let nd = d + c;
                if dp[mask][w].map_or(true, |cur| nd < cur) {
                    dp[mask][w] = Some(nd);
                    how[mask][w] = How::Edge(id);
                    heap.push(Reverse((nd, w)));
                }
            }
        }
    }
    let root = terminals[0];
    let total = dp[full][root].expect("the terminals are not connected");
    // read the tree back
    let mut chosen = Vec::new();
    let mut stack = vec![(full, root)];
    while let Some((mask, v)) = stack.pop() {
        match how[mask][v] {
            How::Terminal => {}
            How::Split(sub) => {
                stack.push((sub as usize, v));
                stack.push((mask ^ sub as usize, v));
            }
            How::Edge(id) => {
                chosen.push(id);
                let (a, b, _) = edges[id];
                stack.push((mask, if a == v { b } else { a }));
            }
        }
    }
    chosen.sort_unstable();
    chosen.dedup();
    (total, chosen)
}
