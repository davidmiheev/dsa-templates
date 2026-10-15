//! Matching algorithms.
//!
//! * [`hopcroft_karp`]: maximum matching in a bipartite graph, `O(m sqrt(n))`.
//! * [`blossom`]: maximum matching in a general graph (Edmonds' blossom algorithm), `O(n^3)`.
//! * [`hungarian`]: minimum-cost perfect matching of a square cost matrix, `O(n^3)`.
//!
//! Unmatched vertices carry [`NONE`] in the result arrays.
//!
//! # Example
//!
//! ```
//! use matching::{blossom, hopcroft_karp, hungarian, NONE};
//!
//! // Left {0, 1}, right {0, 1}: edges 0-0, 0-1, 1-0.
//! let m = hopcroft_karp(2, 2, &[(0, 0), (0, 1), (1, 0)]);
//! assert_eq!(m.size, 2);
//! assert_eq!(m.left_to_right, vec![1, 0]);
//!
//! // A triangle plus a pendant vertex: 4 vertices, perfect matching exists.
//! let g = blossom(4, &[(0, 1), (1, 2), (2, 0), (2, 3)]);
//! assert_eq!(g.size, 2);
//! assert_ne!(g.mate[3], NONE);
//!
//! let (cost, assignment) = hungarian(&[vec![4, 1], vec![2, 0]]);
//! assert_eq!((cost, assignment), (3, vec![1, 0]));
//! ```

use csr_graph::Csr;
use std::collections::VecDeque;

/// "Unmatched".
pub const NONE: usize = usize::MAX;

/// A bipartite matching. See [`hopcroft_karp`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BipartiteMatching {
    /// Number of matched pairs.
    pub size: usize,
    /// Partner of every left vertex ([`NONE`] if unmatched).
    pub left_to_right: Vec<usize>,
    /// Partner of every right vertex ([`NONE`] if unmatched).
    pub right_to_left: Vec<usize>,
}

/// Maximum matching of the bipartite graph with `l` left and `r` right vertices and edges
/// `(left, right)` (parallel edges are fine).
///
/// # Panics
/// Panics if an endpoint is out of range.
pub fn hopcroft_karp(l: usize, r: usize, edges: &[(usize, usize)]) -> BipartiteMatching {
    for &(a, b) in edges {
        assert!(a < l && b < r, "edge ({a}, {b}) is outside {l} x {r}");
    }
    // Left vertices are the tails of a directed CSR graph; heads index the right side (so the
    // graph is built over `max(l, r)` vertices and only its first `l` are used as tails).
    let g = Csr::directed(l.max(r).max(1), edges);
    let mut left_to_right = vec![NONE; l];
    let mut right_to_left = vec![NONE; r];
    let mut dist = vec![0u32; l];
    let mut next_arc = vec![0usize; l];
    let mut queue: Vec<usize> = Vec::with_capacity(l);
    let mut size = 0;
    loop {
        // BFS layers from the free left vertices.
        queue.clear();
        for u in 0..l {
            if left_to_right[u] == NONE {
                dist[u] = 0;
                queue.push(u);
            } else {
                dist[u] = u32::MAX;
            }
        }
        let mut head = 0;
        let mut found = false;
        while head < queue.len() {
            let u = queue[head];
            head += 1;
            for a in g.arcs(u) {
                let w = right_to_left[g.to[a]];
                if w == NONE {
                    found = true;
                } else if dist[w] == u32::MAX {
                    dist[w] = dist[u] + 1;
                    queue.push(w);
                }
            }
        }
        if !found {
            break;
        }
        // Iterative DFS for vertex-disjoint shortest augmenting paths.
        next_arc[..l].copy_from_slice(&g.start[..l]);
        for root in 0..l {
            if left_to_right[root] != NONE {
                continue;
            }
            let mut path: Vec<usize> = vec![root];
            'search: while let Some(&u) = path.last() {
                while next_arc[u] < g.start[u + 1] {
                    let a = next_arc[u];
                    next_arc[u] += 1;
                    let v = g.to[a];
                    let w = right_to_left[v];
                    if w == NONE {
                        // Augment along `path`: each left vertex takes the right vertex it pointed at.
                        let mut right = v;
                        while let Some(x) = path.pop() {
                            let previous = left_to_right[x];
                            left_to_right[x] = right;
                            right_to_left[right] = x;
                            right = previous;
                        }
                        size += 1;
                        break 'search;
                    } else if dist[w] == dist[u] + 1 {
                        path.push(w);
                        continue 'search;
                    }
                }
                dist[u] = u32::MAX; // dead end
                path.pop();
            }
        }
    }
    BipartiteMatching { size, left_to_right, right_to_left }
}

/// A matching in a general graph. See [`blossom`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Matching {
    /// Number of matched pairs.
    pub size: usize,
    /// Partner of every vertex ([`NONE`] if unmatched).
    pub mate: Vec<usize>,
}

/// Maximum matching of the simple undirected graph on `0..n` (Edmonds' blossom algorithm,
/// `O(n^3)`; self-loops are ignored, parallel edges are fine).
///
/// # Panics
/// Panics if an endpoint is `>= n`.
pub fn blossom(n: usize, edges: &[(usize, usize)]) -> Matching {
    let g = Csr::undirected(n, edges);
    let mut mate = vec![NONE; n];
    // Greedy start.
    for u in 0..n {
        if mate[u] == NONE {
            if let Some(&v) = g.neighbors(u).iter().find(|&&v| v != u && mate[v] == NONE) {
                mate[u] = v;
                mate[v] = u;
            }
        }
    }
    let mut parent = vec![NONE; n];
    let mut base: Vec<usize> = (0..n).collect();
    let mut used = vec![false; n];
    let mut in_blossom = vec![false; n];
    let mut queue: VecDeque<usize> = VecDeque::new();
    for root in 0..n {
        if mate[root] != NONE {
            continue;
        }
        // Search for an augmenting path from `root`.
        used.fill(false);
        parent.fill(NONE);
        for (i, b) in base.iter_mut().enumerate() {
            *b = i;
        }
        used[root] = true;
        queue.clear();
        queue.push_back(root);
        let mut end = NONE;
        'bfs: while let Some(v) = queue.pop_front() {
            for &to in g.neighbors(v) {
                if base[v] == base[to] || mate[v] == to || to == v {
                    continue;
                }
                if to == root || (mate[to] != NONE && parent[mate[to]] != NONE) {
                    // Found a blossom: contract it.
                    let cur = lowest_common_base(v, to, &mate, &parent, &base, n);
                    in_blossom.fill(false);
                    mark_path(v, cur, to, &mate, &mut parent, &base, &mut in_blossom);
                    mark_path(to, cur, v, &mate, &mut parent, &base, &mut in_blossom);
                    for i in 0..n {
                        if in_blossom[base[i]] {
                            base[i] = cur;
                            if !used[i] {
                                used[i] = true;
                                queue.push_back(i);
                            }
                        }
                    }
                } else if parent[to] == NONE {
                    parent[to] = v;
                    if mate[to] == NONE {
                        end = to;
                        break 'bfs;
                    }
                    let next = mate[to];
                    used[next] = true;
                    queue.push_back(next);
                }
            }
        }
        // Augment.
        let mut v = end;
        while v != NONE {
            let pv = parent[v];
            let ppv = mate[pv];
            mate[v] = pv;
            mate[pv] = v;
            v = ppv;
        }
    }
    let size = mate.iter().filter(|&&m| m != NONE).count() / 2;
    Matching { size, mate }
}

fn lowest_common_base(mut a: usize, mut b: usize, mate: &[usize], parent: &[usize], base: &[usize], n: usize) -> usize {
    let mut seen = vec![false; n];
    loop {
        a = base[a];
        seen[a] = true;
        if mate[a] == NONE {
            break;
        }
        a = parent[mate[a]];
    }
    loop {
        b = base[b];
        if seen[b] {
            return b;
        }
        b = parent[mate[b]];
    }
}

fn mark_path(mut v: usize, b: usize, mut child: usize, mate: &[usize], parent: &mut [usize], base: &[usize], in_blossom: &mut [bool]) {
    while base[v] != b {
        in_blossom[base[v]] = true;
        in_blossom[base[mate[v]]] = true;
        parent[v] = child;
        child = mate[v];
        v = parent[mate[v]];
    }
}

/// Minimum-cost perfect matching of the square `n x n` cost matrix (Hungarian algorithm with
/// potentials, `O(n^3)`): the total cost and `assignment[row] = column`. Costs may be negative.
///
/// # Panics
/// Panics if the matrix is empty or not square.
pub fn hungarian(cost: &[Vec<i64>]) -> (i64, Vec<usize>) {
    let n = cost.len();
    assert!(n > 0 && cost.iter().all(|row| row.len() == n), "the cost matrix must be square and non-empty");
    const INF: i64 = i64::MAX / 4;
    // 1-based arrays; column 0 is a sentinel.
    let (mut u, mut v) = (vec![0i64; n + 1], vec![0i64; n + 1]);
    let (mut p, mut way) = (vec![0usize; n + 1], vec![0usize; n + 1]);
    for i in 1..=n {
        p[0] = i;
        let mut j0 = 0;
        let mut minv = vec![INF; n + 1];
        let mut used = vec![false; n + 1];
        loop {
            used[j0] = true;
            let (i0, mut delta, mut j1) = (p[j0], INF, 0);
            for j in 1..=n {
                if !used[j] {
                    let cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
                    if cur < minv[j] {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if minv[j] < delta {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            for j in 0..=n {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
            if p[j0] == 0 {
                break;
            }
        }
        loop {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break;
            }
        }
    }
    let mut assignment = vec![0usize; n];
    for j in 1..=n {
        assignment[p[j] - 1] = j - 1;
    }
    let total = (0..n).map(|i| cost[i][assignment[i]]).sum();
    (total, assignment)
}
