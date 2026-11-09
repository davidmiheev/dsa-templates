//! Counting structures with determinants, modulo a prime `M`.
//!
//! * [`count_spanning_trees`] — Kirchhoff's matrix-tree theorem for undirected multigraphs.
//! * [`count_arborescences`] — directed version: spanning trees in which every vertex is
//!   reachable from `root`, i.e. every non-root vertex picks one incoming edge.
//! * [`count_eulerian_circuits`] — BEST theorem: `t_w(G) * prod_v (deg(v) - 1)!`, counting
//!   circuits as cyclic edge sequences (equivalently, circuits starting with edge 0).
//!
//! All three are `O(n^3 + m)`.
//!
//! # Example
//!
//! ```
//! use modint::ModInt998244353 as Mint;
//! use spanning_count::{count_arborescences, count_eulerian_circuits, count_spanning_trees};
//!
//! // K4 has 4^2 = 16 spanning trees (Cayley).
//! let k4: Vec<(usize, usize)> = (0..4).flat_map(|u| (u + 1..4).map(move |v| (u, v))).collect();
//! assert_eq!(count_spanning_trees::<998_244_353>(4, &k4), Mint::new(16));
//! // A directed 3-cycle has exactly one arborescence per root and one Eulerian circuit.
//! let cyc = [(0, 1), (1, 2), (2, 0)];
//! assert_eq!(count_arborescences::<998_244_353>(3, 0, &cyc), Mint::new(1));
//! assert_eq!(count_eulerian_circuits::<998_244_353>(3, &cyc), Mint::new(1));
//! ```

use modint::StaticModInt;
use modmatrix::determinant;

type Mint<const M: u32> = StaticModInt<M>;

/// Number of spanning trees of the undirected multigraph `edges` on `0..n` (loops ignored,
/// parallel edges counted separately), modulo `M`. A single vertex has exactly one.
pub fn count_spanning_trees<const M: u32>(n: usize, edges: &[(usize, usize)]) -> Mint<M> {
    if n == 0 {
        return Mint::one();
    }
    let mut lap = vec![vec![Mint::<M>::zero(); n]; n];
    for &(u, v) in edges {
        assert!(u < n && v < n, "edge ({u}, {v}) has an endpoint outside 0..{n}");
        if u != v {
            lap[u][u] += Mint::one();
            lap[v][v] += Mint::one();
            lap[u][v] -= Mint::one();
            lap[v][u] -= Mint::one();
        }
    }
    lap.pop();
    for row in &mut lap {
        row.pop();
    }
    determinant(lap)
}

/// Number of spanning arborescences rooted at `root` (all vertices reachable from it through
/// directed tree edges) of the directed multigraph `edges`, modulo `M`.
pub fn count_arborescences<const M: u32>(n: usize, root: usize, edges: &[(usize, usize)]) -> Mint<M> {
    assert!(root < n, "root {root} out of range for {n} vertices");
    let mut lap = vec![vec![Mint::<M>::zero(); n]; n];
    for &(u, v) in edges {
        assert!(u < n && v < n, "edge ({u}, {v}) has an endpoint outside 0..{n}");
        if u != v {
            lap[v][v] += Mint::one(); // in-degree
            lap[u][v] -= Mint::one();
        }
    }
    lap.remove(root);
    for row in &mut lap {
        row.remove(root);
    }
    determinant(lap)
}

/// Number of Eulerian circuits of the directed multigraph `edges` on `0..n`, modulo `M`, with
/// circuits that differ by a cyclic shift identified. Returns 0 if there is no circuit; a graph
/// without edges has none to count either (the caller guarantees at least one edge).
pub fn count_eulerian_circuits<const M: u32>(n: usize, edges: &[(usize, usize)]) -> Mint<M> {
    assert!(!edges.is_empty(), "the graph needs at least one edge");
    let (mut indeg, mut outdeg) = (vec![0usize; n], vec![0usize; n]);
    let mut uf = union_find::plain(n);
    for &(u, v) in edges {
        assert!(u < n && v < n, "edge ({u}, {v}) has an endpoint outside 0..{n}");
        outdeg[u] += 1;
        indeg[v] += 1;
        uf.merge(u, v);
    }
    let start = edges[0].0;
    for v in 0..n {
        if indeg[v] != outdeg[v] || (indeg[v] > 0 && !uf.are_connected(v, start)) {
            return Mint::zero();
        }
    }
    // Restrict to vertices with edges: isolated ones would make the Laplacian singular.
    let used: Vec<usize> = (0..n).filter(|&v| outdeg[v] > 0).collect();
    let mut id = vec![usize::MAX; n];
    for (i, &v) in used.iter().enumerate() {
        id[v] = i;
    }
    let sub: Vec<(usize, usize)> = edges.iter().map(|&(u, v)| (id[u], id[v])).collect();
    let trees = count_arborescences::<M>(used.len(), id[start], &sub);
    let mut fact = vec![Mint::<M>::one(); edges.len() + 1];
    for i in 1..fact.len() {
        fact[i] = fact[i - 1] * Mint::from(i);
    }
    used.iter().fold(trees, |acc, &v| acc * fact[outdeg[v] - 1])
}
