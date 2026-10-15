//! Minimum diameter spanning tree of a connected weighted undirected graph (Hassin–Tamir):
//! the optimum is a shortest-path tree rooted at the **absolute centre** — the point of the graph
//! (a vertex or an interior point of an edge) minimising the largest distance to a vertex — and
//! its diameter is twice that radius.
//!
//! With all-pairs distances (`n` Dijkstra runs), an edge `{u, v}` of weight `w` is handled by the
//! classic sweep: sort the vertices by `d(u, .)` descending, keep the *skyline* of vertices with
//! strictly increasing `d(v, .)`, and every consecutive skyline pair `(j, j+1)` yields the best
//! point of the edge between the descending leg of `j` and the ascending leg of `j+1`, of doubled
//! radius `d(v, j) + d(u, j+1) + w`. The whole thing is `O(n m log)` for the distances plus
//! `O(n m)` for the sweeps.
//!
//! Everything is computed with *doubled* quantities, so no division is needed and integer
//! weights stay exact.
//!
//! # Example
//!
//! ```
//! use graph_center::minimum_diameter_spanning_tree;
//!
//! // A path 0 - 1 - 2 with an extra long edge 0 - 2: the best tree is the path.
//! let (diameter, tree) = minimum_diameter_spanning_tree(3, &[(0, 1, 2u64), (1, 2, 2), (0, 2, 9)]);
//! assert_eq!(diameter, 4);
//! assert_eq!(tree, vec![0, 1]);
//! ```

use csr_graph::Csr;
use dijkstra::dijkstra;
use std::ops::{Add, Sub};

const NONE: usize = usize::MAX;

/// A minimum-diameter spanning tree of the connected multigraph `edges = (u, v, weight)` on
/// `0..n` with non-negative weights. Returns `(diameter, sorted edge indices)`; the tree has
/// `n - 1` edges.
///
/// # Panics
/// Panics if the graph is not connected or `n == 0`.
pub fn minimum_diameter_spanning_tree<W>(n: usize, edges: &[(usize, usize, W)]) -> (W, Vec<usize>)
where
    W: Copy + Ord + Default + Add<Output = W> + Sub<Output = W>,
{
    assert!(n >= 1, "need at least one vertex");
    let zero = W::default();
    let m = edges.len();
    let pairs: Vec<(usize, usize)> = edges.iter().map(|&(u, v, _)| (u, v)).collect();
    let weight: Vec<W> = edges.iter().map(|&(_, _, w)| w).collect();
    let g = Csr::undirected(n, &pairs);
    let dist: Vec<Vec<W>> = (0..n)
        .map(|s| {
            let sp = dijkstra(&g, &weight, s);
            assert!((0..n).all(|v| sp.reachable(v)), "the graph is not connected");
            sp.dist
        })
        .collect();
    let order: Vec<Vec<u32>> = (0..n)
        .map(|u| {
            let mut o: Vec<u32> = (0..n as u32).collect();
            o.sort_unstable_by(|&a, &b| dist[u][b as usize].cmp(&dist[u][a as usize]));
            o
        })
        .collect();
    // Best centre so far: doubled radius = diameter; `centre_edge == NONE` means a vertex.
    let mut best = {
        let ecc = |u: usize| dist[u][order[u][0] as usize];
        let (u, e) = (0..n).map(|u| (u, ecc(u))).min_by_key(|&(_, e)| e).unwrap();
        (e + e, NONE, zero, u) // (diameter, edge, doubled offset from the edge's first end, vertex)
    };
    for (i, &(u, v, w)) in edges.iter().enumerate() {
        if u == v {
            continue;
        }
        let (da, db) = (&dist[u], &dist[v]);
        let mut skyline: Vec<usize> = Vec::new();
        let mut max_b: Option<W> = None;
        for &k in &order[u] {
            let k = k as usize;
            if max_b < Some(db[k]) {
                // the previous skyline vertex with the same distance from u is dominated
                if let Some(&last) = skyline.last() {
                    if da[last] == da[k] {
                        skyline.pop();
                    }
                }
                skyline.push(k);
                max_b = Some(db[k]);
            }
        }
        for pair in skyline.windows(2) {
            let (j, l) = (pair[0], pair[1]);
            // the descending leg of j meets the ascending leg of l
            let reach = db[j] + w; // d(v, j) + w
            if da[l] > reach {
                continue; // the crossing would lie before u
            }
            let twice_x = reach - da[l];
            if twice_x > w + w {
                continue; // ... or beyond v
            }
            let diameter = db[j] + da[l] + w;
            if diameter < best.0 {
                best = (diameter, i, twice_x, u);
            }
        }
    }
    // Shortest-path tree from the centre.
    let (_, centre_edge, twice_x, centre_vertex) = best;
    let mut tree: Vec<usize> = Vec::with_capacity(n - 1);
    if centre_edge == NONE {
        let sp = dijkstra(&g, &weight, centre_vertex);
        for v in 0..n {
            if v != centre_vertex {
                tree.push(sp.parent_edge[v]);
            }
        }
    } else {
        // Subdivide the centre edge with a new vertex `n` (all weights doubled).
        let (u, v, w) = edges[centre_edge];
        let mut pairs2 = pairs.clone();
        let mut weight2: Vec<W> = weight.iter().map(|&x| x + x).collect();
        pairs2[centre_edge] = (n, u);
        weight2[centre_edge] = twice_x;
        pairs2.push((n, v));
        weight2.push(w + w - twice_x);
        let g2 = Csr::undirected(n + 1, &pairs2);
        let sp = dijkstra(&g2, &weight2, n);
        let (mut used_u, mut used_v) = (false, false);
        for x in 0..n {
            let e = sp.parent_edge[x];
            if e == centre_edge {
                used_u = true;
            } else if e == m {
                used_v = true;
            } else {
                tree.push(e);
            }
        }
        if used_u && used_v {
            tree.push(centre_edge);
        }
    }
    tree.sort_unstable();
    tree.dedup();
    assert_eq!(tree.len(), n - 1, "the shortest-path tree must span all vertices");
    // Actual diameter of the tree (two sweeps).
    let tree_pairs: Vec<(usize, usize)> = tree.iter().map(|&e| pairs[e]).collect();
    let tree_weight: Vec<W> = tree.iter().map(|&e| weight[e]).collect();
    let tg = Csr::undirected(n, &tree_pairs);
    let far = |s: usize| {
        let sp = dijkstra(&tg, &tree_weight, s);
        (0..n).map(|v| (sp.dist[v], v)).max().unwrap()
    };
    let (_, a) = far(0);
    let (diameter, _) = far(a);
    (diameter, tree)
}
