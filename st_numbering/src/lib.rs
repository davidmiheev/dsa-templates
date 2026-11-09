//! st-numbering: a numbering of the vertices `0..n` such that `s` gets the smallest number, `t`
//! the largest, and every other vertex has a neighbour with a smaller and one with a larger
//! number. Equivalently, orienting every edge from the smaller to the larger number puts every
//! vertex on some `s -> t` path. It exists iff the graph plus the edge `s-t` is biconnected.
//!
//! The algorithm is Tarjan's: a DFS from `s` whose first tree edge is the (virtual) edge `s-t`,
//! then the vertices are inserted into a linked list in preorder — before or after their parent,
//! depending on a sign that records on which side the lowpoint of the vertex was placed.
//! `O(n + m)`, no recursion.
//!
//! # Example
//!
//! ```
//! use st_numbering::st_numbering;
//!
//! // A 4-cycle 0-1-2-3 with s = 0 and t = 2.
//! let p = st_numbering(4, &[(0, 1), (1, 2), (2, 3), (3, 0)], 0, 2);
//! assert_eq!((p[0], p[2]), (0, 3));
//! // A path is not biconnected: no numbering (empty result).
//! assert!(st_numbering(3, &[(0, 1), (1, 2)], 0, 1).is_empty());
//! ```

use csr_graph::Csr;

const NONE: usize = usize::MAX;

/// A valid numbering `p` (a permutation of `0..n`, `p[s] = 0`, `p[t] = n - 1`) of the undirected
/// multigraph `edges`, or an **empty vector** if none exists. Requires `s != t`.
///
/// # Panics
/// Panics if `s == t` or a vertex is out of range.
pub fn st_numbering(n: usize, edges: &[(usize, usize)], s: usize, t: usize) -> Vec<usize> {
    assert!(s < n && t < n && s != t, "need two distinct terminals inside 0..{n}");
    // With two vertices the virtual edge is not enough: some real edge has to carry the path.
    if n == 2 && edges.is_empty() {
        return Vec::new();
    }
    // The virtual edge s-t goes first so that the DFS takes it from `s`.
    let virtual_id = edges.len();
    let mut all = Vec::with_capacity(edges.len() + 1);
    all.push((s, t));
    all.extend_from_slice(edges);
    let g = Csr::undirected(n, &all);
    let eid = |a: usize| if g.id[a] == 0 { virtual_id } else { g.id[a] - 1 };

    // Iterative DFS with preorder, parent, lowpoint vertex.
    let mut pre = vec![NONE; n];
    let mut parent = vec![NONE; n];
    let mut parent_edge = vec![NONE; n];
    let mut low: Vec<usize> = (0..n).collect();
    let mut order = Vec::with_capacity(n);
    let mut next_arc = g.start[..n].to_vec();
    let mut stack = vec![s];
    pre[s] = 0;
    order.push(s);
    while let Some(&v) = stack.last() {
        if next_arc[v] < g.start[v + 1] {
            let a = next_arc[v];
            next_arc[v] += 1;
            let w = g.to[a];
            if w == v || eid(a) == parent_edge[v] {
                continue;
            }
            if pre[w] == NONE {
                pre[w] = order.len();
                parent[w] = v;
                parent_edge[w] = eid(a);
                order.push(w);
                stack.push(w);
            } else if pre[w] < pre[v] && pre[w] < pre[low[v]] {
                low[v] = w;
            }
        } else {
            stack.pop();
            let p = parent[v];
            if p != NONE && pre[low[v]] < pre[low[p]] {
                low[p] = low[v];
            }
        }
    }
    // Connected, and biconnected once s-t is present: the root has one child (t) and no other
    // vertex has a child whose subtree cannot escape above it.
    if order.len() != n || parent[t] != s {
        return Vec::new();
    }
    for &v in &order {
        let p = parent[v];
        if p == NONE || v == t {
            continue;
        }
        if p == s || pre[low[v]] >= pre[p] {
            return Vec::new();
        }
    }
    // Build the list.
    let (mut prev, mut next) = (vec![NONE; n], vec![NONE; n]);
    next[s] = t;
    prev[t] = s;
    let mut minus = vec![false; n]; // sign of every vertex
    minus[s] = true;
    for &v in &order {
        if v == s || v == t {
            continue;
        }
        let p = parent[v];
        if minus[low[v]] {
            // insert v before p
            let q = prev[p];
            prev[v] = q;
            next[v] = p;
            prev[p] = v;
            if q != NONE {
                next[q] = v;
            }
            minus[p] = false;
        } else {
            // insert v after p
            let q = next[p];
            next[v] = q;
            prev[v] = p;
            next[p] = v;
            if q != NONE {
                prev[q] = v;
            }
            minus[p] = true;
        }
    }
    let mut num = vec![0usize; n];
    let (mut v, mut k) = (s, 0);
    while v != NONE {
        num[v] = k;
        k += 1;
        v = next[v];
    }
    debug_assert_eq!(k, n);
    num
}
