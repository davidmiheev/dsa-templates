//! Small-subgraph enumeration on sparse graphs by **degree ordering**: rank vertices by
//! `(degree, index)` and only look "upwards" (towards higher rank) or "downwards", so that every
//! vertex has at most `O(sqrt m)` higher-ranked neighbours and the total work is `O(m sqrt m)`.
//!
//! * [`for_each_triangle`] — every triangle exactly once.
//! * [`for_each_clique`] — every non-empty clique exactly once (exponential in the clique
//!   number, as the output is).
//! * [`four_cycles_per_edge`] — for every edge of a multigraph, the number of 4-cycles through it
//!   (cycles of four distinct edges on four distinct vertices).
//!
//! # Example
//!
//! ```
//! use subgraph_counting::{for_each_triangle, four_cycles_per_edge};
//!
//! let mut triangles = Vec::new();
//! for_each_triangle(4, &[(0, 1), (1, 2), (2, 0), (2, 3), (3, 0)], |a, b, c| {
//!     let mut t = [a, b, c];
//!     t.sort();
//!     triangles.push(t);
//! });
//! triangles.sort();
//! assert_eq!(triangles, vec![[0, 1, 2], [0, 2, 3]]);
//!
//! // The square 0-1-2-3 has a single 4-cycle through each edge.
//! assert_eq!(four_cycles_per_edge(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]), vec![1, 1, 1, 1]);
//! ```

use csr_graph::Csr;

/// Simple undirected graph (parallel edges and loops removed) oriented from lower to higher
/// `(degree, index)` rank; `out(v)` lists the higher-ranked neighbours.
fn orient(g: &Csr) -> Vec<Vec<usize>> {
    let n = g.len();
    let rank_less = |a: usize, b: usize| (g.degree(a), a) < (g.degree(b), b);
    let mut out: Vec<Vec<usize>> = vec![Vec::new(); n];
    for v in 0..n {
        for &w in g.neighbors(v) {
            if w != v && rank_less(v, w) {
                out[v].push(w);
            }
        }
        out[v].sort_unstable();
        out[v].dedup();
    }
    out
}

/// Call `f(a, b, c)` once for every triangle of the simple graph `edges` on `0..n` (parallel
/// edges are merged, loops ignored). `O(m sqrt m)`.
pub fn for_each_triangle(n: usize, edges: &[(usize, usize)], mut f: impl FnMut(usize, usize, usize)) {
    let g = Csr::undirected(n, edges);
    let out = orient(&g);
    let mut mark = vec![usize::MAX; n];
    for u in 0..n {
        for &v in &out[u] {
            mark[v] = u;
        }
        for &v in &out[u] {
            for &w in &out[v] {
                if mark[w] == u {
                    f(u, v, w);
                }
            }
        }
    }
}

/// Call `f(clique)` once for every non-empty clique (including single vertices and edges) of
/// the simple graph `edges` on `0..n`. The slice lists the vertices in increasing rank order.
pub fn for_each_clique(n: usize, edges: &[(usize, usize)], mut f: impl FnMut(&[usize])) {
    let g = Csr::undirected(n, edges);
    let out = orient(&g);
    let mut cur = Vec::new();
    // `mark[w] == v` exactly when `w` is a higher-ranked neighbour of `v` (stamped just before use).
    fn rec(cand: &[usize], out: &[Vec<usize>], mark: &mut [usize], cur: &mut Vec<usize>, f: &mut dyn FnMut(&[usize])) {
        for &v in cand {
            cur.push(v);
            f(cur);
            for &w in &out[v] {
                mark[w] = v;
            }
            let next: Vec<usize> = cand.iter().copied().filter(|&w| mark[w] == v).collect();
            if !next.is_empty() {
                rec(&next, out, mark, cur, f);
            }
            cur.pop();
        }
    }
    let mut mark = vec![usize::MAX; n];
    for v in 0..n {
        cur.push(v);
        f(&cur);
        if !out[v].is_empty() {
            rec(&out[v], &out, &mut mark, &mut cur, &mut f);
        }
        cur.pop();
    }
}

/// For every edge of the undirected multigraph `edges` (no loops), the number of 4-cycles
/// containing it: unordered sets of four distinct edges forming a cycle through four distinct
/// vertices. Parallel edges multiply the counts as expected. `O(m sqrt m)`.
///
/// # Panics
/// Panics on a self-loop or an out-of-range endpoint.
pub fn four_cycles_per_edge(n: usize, edges: &[(usize, usize)]) -> Vec<u64> {
    for &(u, v) in edges {
        assert!(u < n && v < n && u != v, "edge ({u}, {v}) is a loop or out of range");
    }
    // Merge parallel edges: pair id and multiplicity.
    let mut keyed: Vec<(usize, usize, usize)> = edges
        .iter()
        .enumerate()
        .map(|(i, &(u, v))| (u.min(v), u.max(v), i))
        .collect();
    keyed.sort_unstable();
    let mut pair_of_edge = vec![0usize; edges.len()];
    let mut pairs: Vec<(usize, usize, u64)> = Vec::new();
    for &(u, v, i) in &keyed {
        match pairs.last_mut() {
            Some(last) if last.0 == u && last.1 == v => last.2 += 1,
            _ => pairs.push((u, v, 1)),
        }
        pair_of_edge[i] = pairs.len() - 1;
    }
    let simple: Vec<(usize, usize)> = pairs.iter().map(|&(u, v, _)| (u, v)).collect();
    let g = Csr::undirected(n, &simple);
    let weight = |arc: usize| pairs[g.id[arc]].2;
    let rank_less = |a: usize, b: usize| (g.degree(a), a) < (g.degree(b), b);
    let mut value = vec![0u64; pairs.len()];
    let mut cnt = vec![0u64; n];
    let mut touched: Vec<usize> = Vec::new();
    for a in 0..n {
        // Wedges a - b - c with b, c ranked below a.
        for ab in g.arcs(a) {
            let b = g.to[ab];
            if !rank_less(b, a) {
                continue;
            }
            let wab = weight(ab);
            for bc in g.arcs(b) {
                let c = g.to[bc];
                if c == a || !rank_less(c, a) {
                    continue;
                }
                if cnt[c] == 0 {
                    touched.push(c);
                }
                cnt[c] += wab * weight(bc);
            }
        }
        for ab in g.arcs(a) {
            let b = g.to[ab];
            if !rank_less(b, a) {
                continue;
            }
            let wab = weight(ab);
            for bc in g.arcs(b) {
                let c = g.to[bc];
                if c == a || !rank_less(c, a) {
                    continue;
                }
                let wbc = weight(bc);
                let other = cnt[c] - wab * wbc; // wedges through a different middle vertex
                value[g.id[ab]] += wbc * other;
                value[g.id[bc]] += wab * other;
            }
        }
        for c in touched.drain(..) {
            cnt[c] = 0;
        }
    }
    pair_of_edge.iter().map(|&p| value[p]).collect()
}
