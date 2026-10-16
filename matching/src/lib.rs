//! Matching algorithms.
//!
//! * [`hopcroft_karp`]: maximum matching in a bipartite graph, `O(m sqrt(n))`.
//! * [`blossom`]: maximum matching in a general graph (Edmonds' blossom algorithm), `O(n^3)`.
//! * [`hungarian`]: minimum-cost perfect matching of a square cost matrix, `O(n^3)`.
//! * [`max_weight_matching`]: maximum-weight matching of a general graph (weighted blossom), `O(n^3)`.
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

/// A maximum-weight matching in a general graph. See [`max_weight_matching`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WeightedMatching {
    /// Total weight of the matched edges.
    pub weight: i64,
    /// Number of matched pairs.
    pub size: usize,
    /// Partner of every vertex ([`NONE`] if unmatched).
    pub mate: Vec<usize>,
}

/// Maximum-**weight** matching (not necessarily perfect) of the simple undirected graph
/// `edges = (u, v, w)` on `0..n` with positive weights, by Edmonds' weighted blossom algorithm
/// with dual variables, `O(n^3)`. Parallel edges keep the heaviest one; self-loops are ignored.
///
/// # Panics
/// Panics on an out-of-range endpoint or a non-positive weight.
pub fn max_weight_matching(n: usize, edges: &[(usize, usize, i64)]) -> WeightedMatching {
    let mut b = WeightedBlossom::new(n);
    for &(u, v, w) in edges {
        assert!(u < n && v < n, "edge ({u}, {v}) has an endpoint outside 0..{n}");
        assert!(w > 0, "weights must be positive, got {w}");
        if u != v {
            let (a, c) = (u + 1, v + 1);
            if w > b.g[a * b.dim + c].w {
                b.g[a * b.dim + c].w = w;
                b.g[c * b.dim + a].w = w;
            }
        }
    }
    b.solve()
}

#[derive(Clone, Copy, Default)]
struct WEdge {
    u: usize,
    v: usize,
    w: i64,
}

/// State of the weighted blossom search. Vertices are `1..=n`, blossoms `n+1..=n_x`; index 0 is
/// the null vertex. (A port of the classical O(n^3) formulation with explicit `flower` lists.)
struct WeightedBlossom {
    n: usize,
    nx: usize,
    dim: usize,
    g: Vec<WEdge>,
    lab: Vec<i64>,
    mate: Vec<usize>,
    slack: Vec<usize>,
    st: Vec<usize>,
    pa: Vec<usize>,
    flower_from: Vec<usize>,
    s: Vec<i8>,
    vis: Vec<u32>,
    flower: Vec<Vec<usize>>,
    queue: VecDeque<usize>,
    stamp: u32,
}

impl WeightedBlossom {
    fn new(n: usize) -> Self {
        let dim = 2 * n + 2;
        let mut g = vec![WEdge::default(); dim * dim];
        for u in 1..=n {
            for v in 1..=n {
                g[u * dim + v] = WEdge { u, v, w: 0 };
            }
        }
        Self {
            n,
            nx: n,
            dim,
            g,
            lab: vec![0; dim],
            mate: vec![0; dim],
            slack: vec![0; dim],
            st: vec![0; dim],
            pa: vec![0; dim],
            flower_from: vec![0; dim * (n + 1)],
            s: vec![0; dim],
            vis: vec![0; dim],
            flower: vec![Vec::new(); dim],
            queue: VecDeque::new(),
            stamp: 0,
        }
    }

    fn e(&self, u: usize, v: usize) -> WEdge {
        self.g[u * self.dim + v]
    }

    fn dist(&self, e: WEdge) -> i64 {
        self.lab[e.u] + self.lab[e.v] - self.e(e.u, e.v).w * 2
    }

    fn ff(&self, b: usize, x: usize) -> usize {
        self.flower_from[b * (self.n + 1) + x]
    }

    fn update_slack(&mut self, u: usize, x: usize) {
        if self.slack[x] == 0 || self.dist(self.e(u, x)) < self.dist(self.e(self.slack[x], x)) {
            self.slack[x] = u;
        }
    }

    fn set_slack(&mut self, x: usize) {
        self.slack[x] = 0;
        for u in 1..=self.n {
            if self.e(u, x).w > 0 && self.st[u] != x && self.s[self.st[u]] == 0 {
                self.update_slack(u, x);
            }
        }
    }

    fn q_push(&mut self, x: usize) {
        if x <= self.n {
            self.queue.push_back(x);
        } else {
            for i in 0..self.flower[x].len() {
                let y = self.flower[x][i];
                self.q_push(y);
            }
        }
    }

    fn set_st(&mut self, x: usize, b: usize) {
        self.st[x] = b;
        if x > self.n {
            for i in 0..self.flower[x].len() {
                let y = self.flower[x][i];
                self.set_st(y, b);
            }
        }
    }

    fn get_pr(&mut self, b: usize, xr: usize) -> usize {
        let pr = self.flower[b].iter().position(|&y| y == xr).unwrap();
        if pr % 2 == 1 {
            self.flower[b][1..].reverse();
            self.flower[b].len() - pr
        } else {
            pr
        }
    }

    fn set_match(&mut self, u: usize, v: usize) {
        self.mate[u] = self.e(u, v).v;
        if u > self.n {
            let e = self.e(u, v);
            let xr = self.ff(u, e.u);
            let pr = self.get_pr(u, xr);
            for i in 0..pr {
                let (a, b) = (self.flower[u][i], self.flower[u][i ^ 1]);
                self.set_match(a, b);
            }
            self.set_match(xr, v);
            self.flower[u].rotate_left(pr);
        }
    }

    fn augment(&mut self, mut u: usize, mut v: usize) {
        loop {
            let xnv = self.st[self.mate[u]];
            self.set_match(u, v);
            if xnv == 0 {
                return;
            }
            let t = self.st[self.pa[xnv]];
            self.set_match(xnv, t);
            u = self.st[self.pa[xnv]];
            v = xnv;
        }
    }

    fn get_lca(&mut self, mut u: usize, mut v: usize) -> usize {
        self.stamp += 1;
        while u != 0 || v != 0 {
            if u != 0 {
                if self.vis[u] == self.stamp {
                    return u;
                }
                self.vis[u] = self.stamp;
                u = self.st[self.mate[u]];
                if u != 0 {
                    u = self.st[self.pa[u]];
                }
            }
            std::mem::swap(&mut u, &mut v);
        }
        0
    }

    fn add_blossom(&mut self, u: usize, lca: usize, v: usize) {
        let mut b = self.n + 1;
        while b <= self.nx && self.st[b] != 0 {
            b += 1;
        }
        if b > self.nx {
            self.nx += 1;
        }
        self.lab[b] = 0;
        self.s[b] = 0;
        self.mate[b] = self.mate[lca];
        self.flower[b].clear();
        self.flower[b].push(lca);
        let mut x = u;
        while x != lca {
            self.flower[b].push(x);
            let y = self.st[self.mate[x]];
            self.flower[b].push(y);
            self.q_push(y);
            x = self.st[self.pa[y]];
        }
        self.flower[b][1..].reverse();
        let mut x = v;
        while x != lca {
            self.flower[b].push(x);
            let y = self.st[self.mate[x]];
            self.flower[b].push(y);
            self.q_push(y);
            x = self.st[self.pa[y]];
        }
        self.set_st(b, b);
        for x in 1..=self.nx {
            self.g[b * self.dim + x].w = 0;
            self.g[x * self.dim + b].w = 0;
        }
        for x in 1..=self.n {
            self.flower_from[b * (self.n + 1) + x] = 0;
        }
        for i in 0..self.flower[b].len() {
            let xs = self.flower[b][i];
            for x in 1..=self.nx {
                if self.e(b, x).w == 0 || self.dist(self.e(xs, x)) < self.dist(self.e(b, x)) {
                    self.g[b * self.dim + x] = self.e(xs, x);
                    self.g[x * self.dim + b] = self.e(x, xs);
                }
            }
            for x in 1..=self.n {
                if self.ff(xs, x) != 0 {
                    self.flower_from[b * (self.n + 1) + x] = xs;
                }
            }
        }
        self.set_slack(b);
    }

    fn expand_blossom(&mut self, b: usize) {
        for i in 0..self.flower[b].len() {
            let y = self.flower[b][i];
            self.set_st(y, y);
        }
        let xr = self.ff(b, self.e(b, self.pa[b]).u);
        let pr = self.get_pr(b, xr);
        let mut i = 0;
        while i < pr {
            let xs = self.flower[b][i];
            let xns = self.flower[b][i + 1];
            self.pa[xs] = self.e(xns, xs).u;
            self.s[xs] = 1;
            self.s[xns] = 0;
            self.slack[xs] = 0;
            self.set_slack(xns);
            self.q_push(xns);
            i += 2;
        }
        self.s[xr] = 1;
        self.pa[xr] = self.pa[b];
        for i in pr + 1..self.flower[b].len() {
            let xs = self.flower[b][i];
            self.s[xs] = -1;
            self.set_slack(xs);
        }
        self.st[b] = 0;
    }

    fn on_found_edge(&mut self, e: WEdge) -> bool {
        let (u, v) = (self.st[e.u], self.st[e.v]);
        if self.s[v] == -1 {
            self.pa[v] = e.u;
            self.s[v] = 1;
            let nu = self.st[self.mate[v]];
            self.slack[v] = 0;
            self.slack[nu] = 0;
            self.s[nu] = 0;
            self.q_push(nu);
        } else if self.s[v] == 0 {
            let lca = self.get_lca(u, v);
            if lca == 0 {
                self.augment(u, v);
                self.augment(v, u);
                return true;
            }
            self.add_blossom(u, lca, v);
        }
        false
    }

    /// One phase: find an augmenting path of the equality graph, adjusting the duals as needed.
    fn matching(&mut self) -> bool {
        for x in 1..=self.nx {
            self.s[x] = -1;
            self.slack[x] = 0;
        }
        self.queue.clear();
        for x in 1..=self.nx {
            if self.st[x] == x && self.mate[x] == 0 {
                self.pa[x] = 0;
                self.s[x] = 0;
                self.q_push(x);
            }
        }
        if self.queue.is_empty() {
            return false;
        }
        loop {
            while let Some(u) = self.queue.pop_front() {
                if self.s[self.st[u]] == 1 {
                    continue;
                }
                for v in 1..=self.n {
                    if self.e(u, v).w > 0 && self.st[u] != self.st[v] {
                        if self.dist(self.e(u, v)) == 0 {
                            if self.on_found_edge(self.e(u, v)) {
                                return true;
                            }
                        } else {
                            let sv = self.st[v];
                            self.update_slack(u, sv);
                        }
                    }
                }
            }
            let mut d = i64::MAX / 4;
            for b in self.n + 1..=self.nx {
                if self.st[b] == b && self.s[b] == 1 {
                    d = d.min(self.lab[b] / 2);
                }
            }
            for x in 1..=self.nx {
                if self.st[x] == x && self.slack[x] != 0 {
                    let sl = self.dist(self.e(self.slack[x], x));
                    if self.s[x] == -1 {
                        d = d.min(sl);
                    } else if self.s[x] == 0 {
                        d = d.min(sl / 2);
                    }
                }
            }
            for u in 1..=self.n {
                let su = self.s[self.st[u]];
                if su == 0 {
                    if self.lab[u] <= d {
                        return false;
                    }
                    self.lab[u] -= d;
                } else if su == 1 {
                    self.lab[u] += d;
                }
            }
            for b in self.n + 1..=self.nx {
                if self.st[b] == b {
                    if self.s[self.st[b]] == 0 {
                        self.lab[b] += d * 2;
                    } else if self.s[self.st[b]] == 1 {
                        self.lab[b] -= d * 2;
                    }
                }
            }
            self.queue.clear();
            for x in 1..=self.nx {
                if self.st[x] == x && self.slack[x] != 0 && self.st[self.slack[x]] != x && self.dist(self.e(self.slack[x], x)) == 0 {
                    let e = self.e(self.slack[x], x);
                    if self.on_found_edge(e) {
                        return true;
                    }
                }
            }
            for b in self.n + 1..=self.nx {
                if self.st[b] == b && self.s[b] == 1 && self.lab[b] == 0 {
                    self.expand_blossom(b);
                }
            }
        }
    }

    fn solve(mut self) -> WeightedMatching {
        let n = self.n;
        self.nx = n;
        for u in 0..=n {
            self.st[u] = u;
            self.flower[u].clear();
        }
        let mut w_max = 0;
        for u in 1..=n {
            for v in 1..=n {
                self.flower_from[u * (n + 1) + v] = if u == v { u } else { 0 };
                w_max = w_max.max(self.e(u, v).w);
            }
        }
        for u in 1..=n {
            self.lab[u] = w_max;
        }
        while self.matching() {}
        let mut mate = vec![NONE; n];
        let (mut weight, mut size) = (0i64, 0usize);
        for u in 1..=n {
            if self.mate[u] != 0 {
                mate[u - 1] = self.mate[u] - 1;
                if self.mate[u] < u {
                    weight += self.e(u, self.mate[u]).w;
                    size += 1;
                }
            }
        }
        WeightedMatching { weight, size, mate }
    }
}
