//! 3-edge-connected components of an undirected multigraph, in `O((n + m) log m)`.
//!
//! Two vertices are 3-edge-connected iff no bridge and no pair of edges separates them. The
//! method is a randomised cycle-space hash on a DFS tree (fixed seed, so runs are reproducible;
//! a wrong answer needs a collision of 64-bit values, probability `~ m^2 / 2^64`):
//!
//! * every non-tree edge gets a random 64-bit value; a tree edge gets the XOR of the values of
//!   the non-tree edges covering it (`0` exactly for bridges);
//! * two edges form a cut pair iff their values are equal. In a DFS tree the tree edges of such a
//!   class are nested on one root path `e_1 ⊃ … ⊃ e_t`, and removing the whole class leaves the
//!   pieces "above `e_1`", "between `e_i` and `e_{i+1}`" and "below `e_t`" — the first and last
//!   one are glued together by the covering non-tree edges unless a non-tree edge belongs to
//!   the class itself;
//! * a vertex's piece in every class (and "side" of every bridge) is folded into one 64-bit
//!   fingerprint by a top-down sum along the tree, and equal fingerprints mean one component.
//!
//! # Example
//!
//! ```
//! use three_edge_cc::three_edge_connected_components;
//!
//! // Two triangles sharing the edge 1-2 and a doubled edge: {0,1,2,3} is a K4 minus an edge.
//! let (count, comp) = three_edge_connected_components(4, &[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (1, 2)]);
//! assert_eq!(count, 3);
//! assert_eq!(comp[1], comp[2]);
//! assert_ne!(comp[0], comp[1]);
//! ```

use csr_graph::Csr;

const NONE: usize = usize::MAX;

fn splitmix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

/// Components of the multigraph `edges` on `0..n` (self-loops ignored): `(count, comp)` with
/// `comp[v]` in `0..count`, numbered by smallest member.
pub fn three_edge_connected_components(n: usize, edges: &[(usize, usize)]) -> (usize, Vec<usize>) {
    let g = Csr::undirected(n, edges);
    // Iterative DFS: preorder, parent, tree-edge flags.
    let mut order = Vec::with_capacity(n);
    let mut parent = vec![NONE; n];
    let mut depth = vec![0usize; n];
    let mut seen = vec![false; n];
    let mut is_tree = vec![false; edges.len()];
    let mut next_arc = g.start[..n].to_vec();
    let mut stack = Vec::new();
    for root in 0..n {
        if seen[root] {
            continue;
        }
        seen[root] = true;
        order.push(root);
        stack.push(root);
        while let Some(&v) = stack.last() {
            if next_arc[v] < g.start[v + 1] {
                let a = next_arc[v];
                next_arc[v] += 1;
                let w = g.to[a];
                if !seen[w] {
                    seen[w] = true;
                    parent[w] = v;
                    depth[w] = depth[v] + 1;
                    is_tree[g.id[a]] = true;
                    order.push(w);
                    stack.push(w);
                }
            } else {
                stack.pop();
            }
        }
    }
    // sig[v] = XOR of the values of the non-tree edges covering the edge (parent[v], v).
    let mut sig = vec![0u64; n];
    let mut keyed: Vec<(u64, usize, usize)> = Vec::new(); // (signature, depth, child vertex | NONE for a non-tree edge)
    for (e, &(u, v)) in edges.iter().enumerate() {
        if u != v && !is_tree[e] {
            let r = splitmix(e as u64 ^ 0xA5A5_5A5A).max(1);
            sig[u] ^= r;
            sig[v] ^= r;
            keyed.push((r, 0, NONE));
        }
    }
    for &v in order.iter().rev() {
        if parent[v] != NONE {
            sig[parent[v]] ^= sig[v];
            keyed.push((sig[v], depth[v], v));
        }
    }
    // The bridges are the tree edges with sig 0; every other signature is a class of equal edges.
    let mut val = vec![0u64; n];
    for &(s, _, v) in &keyed {
        if s == 0 && v != NONE {
            val[v] = splitmix(v as u64 ^ 0x1234_5678_9ABC);
        }
    }
    keyed.retain(|&(s, _, _)| s != 0);
    keyed.sort_unstable();
    let mut i = 0;
    while i < keyed.len() {
        let mut j = i;
        while j < keyed.len() && keyed[j].0 == keyed[i].0 {
            j += 1;
        }
        if j - i >= 2 {
            let r = splitmix(keyed[i].0);
            let has_non_tree = keyed[i..j].iter().any(|&(_, _, v)| v == NONE);
            let tree: Vec<&(u64, usize, usize)> = keyed[i..j].iter().filter(|x| x.2 != NONE).collect();
            let deepest = tree.iter().max_by_key(|x| x.1).map(|x| x.2);
            for x in &tree {
                val[x.2] = r;
            }
            if !has_non_tree {
                // Glue the part below the deepest class edge to the part above the top one.
                let d = deepest.expect("a class without a non-tree edge has tree edges");
                val[d] = r.wrapping_mul(1u64.wrapping_sub(tree.len() as u64));
            }
        }
        i = j;
    }
    // Fingerprints, top-down.
    let mut fp = vec![0u64; n];
    for &v in &order {
        fp[v] = if parent[v] == NONE { splitmix(v as u64 ^ 0xFEED_FACE) } else { fp[parent[v]].wrapping_add(val[v]) };
    }
    let mut by_fp: Vec<usize> = (0..n).collect();
    by_fp.sort_unstable_by_key(|&v| (fp[v], v));
    let mut raw = vec![0usize; n];
    let mut label = 0usize;
    for k in 0..n {
        if k > 0 && fp[by_fp[k]] != fp[by_fp[k - 1]] {
            label += 1;
        }
        raw[by_fp[k]] = label;
    }
    // Renumber by smallest member.
    let mut rename = vec![NONE; n];
    let mut comp = vec![0usize; n];
    let mut count = 0;
    for v in 0..n {
        if rename[raw[v]] == NONE {
            rename[raw[v]] = count;
            count += 1;
        }
        comp[v] = rename[raw[v]];
    }
    (count, comp)
}
