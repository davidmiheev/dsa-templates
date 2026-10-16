//! Maximum independent set and maximum clique on graphs with at most 64 vertices, as `u64`
//! adjacency bitmasks and branching:
//!
//! * a vertex of degree `0` or `1` is always taken (degree 1: taking it is never worse than
//!   taking its neighbour);
//! * otherwise branch on a vertex `v` of maximum degree: either take `v` (delete `N[v]`) or
//!   discard it. As long as the maximum degree is at least 3 this is the classic
//!   `T(n) = T(n-1) + T(n-4)` recursion, `O(1.39^n)`, which is instantaneous for `n <= 40`.
//!   When every degree is 2 the graph is a union of cycles and branching on any vertex
//!   immediately creates degree-1 vertices.
//!
//! # Example
//!
//! ```
//! use independent_set::{maximum_clique, maximum_independent_set};
//!
//! // The 5-cycle: independence number 2, clique number 2.
//! let c5 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)];
//! assert_eq!(maximum_independent_set(5, &c5).len(), 2);
//! assert_eq!(maximum_clique(5, &c5).len(), 2);
//! ```

fn masks(n: usize, edges: &[(usize, usize)]) -> Vec<u64> {
    assert!(n <= 64, "bitmask search is limited to 64 vertices, got {n}");
    let mut adj = vec![0u64; n];
    for &(u, v) in edges {
        assert!(u < n && v < n, "edge ({u}, {v}) has an endpoint outside 0..{n}");
        if u != v {
            adj[u] |= 1 << v;
            adj[v] |= 1 << u;
        }
    }
    adj
}

/// Best independent set inside `alive`, as a bitmask.
fn best(adj: &[u64], alive: u64) -> u64 {
    if alive == 0 {
        return 0;
    }
    // find a vertex of degree <= 1 or the one of maximum degree
    let (mut pick, mut pick_deg) = (0usize, 0u32);
    let mut rest = alive;
    while rest != 0 {
        let v = rest.trailing_zeros() as usize;
        rest &= rest - 1;
        let d = (adj[v] & alive).count_ones();
        if d <= 1 {
            return (1 << v) | best(adj, alive & !(1 << v) & !adj[v]);
        }
        if d > pick_deg {
            pick = v;
            pick_deg = d;
        }
    }
    let take = (1u64 << pick) | best(adj, alive & !(1 << pick) & !adj[pick]);
    let skip = best(adj, alive & !(1 << pick));
    if take.count_ones() >= skip.count_ones() {
        take
    } else {
        skip
    }
}

fn bits(mut mask: u64) -> Vec<usize> {
    let mut v = Vec::with_capacity(mask.count_ones() as usize);
    while mask != 0 {
        v.push(mask.trailing_zeros() as usize);
        mask &= mask - 1;
    }
    v
}

/// A maximum independent set of the graph `edges` on `0..n` (`n <= 64`), in increasing order.
pub fn maximum_independent_set(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let adj = masks(n, edges);
    let alive = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
    bits(best(&adj, alive))
}

/// A maximum clique (an independent set of the complement), in increasing order.
pub fn maximum_clique(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let adj = masks(n, edges);
    let alive = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
    let co: Vec<u64> = (0..n).map(|v| !adj[v] & alive & !(1 << v)).collect();
    bits(best(&co, alive))
}
