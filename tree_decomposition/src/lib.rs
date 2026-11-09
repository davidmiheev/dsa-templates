//! Tree decomposition of width at most 2 in `O(n + m)` expected time.
//!
//! A graph has treewidth `<= 2` iff repeatedly deleting a vertex of degree `<= 2` — joining the
//! two neighbours of a degree-2 vertex by an edge — removes everything (no `K4` minor). The
//! deleted vertex `v` together with its neighbours at that moment is a bag of size `<= 3`, and the
//! bag hangs below the bag of the neighbour that is deleted first afterwards (the other
//! neighbour is adjacent to it by construction, so it lies in that bag too). Components are
//! chained together at their roots.
//!
//! # Example
//!
//! ```
//! use tree_decomposition::tree_decomposition_width2;
//!
//! // A 4-cycle has treewidth 2.
//! let d = tree_decomposition_width2(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]);
//! assert!(!d.bags.is_empty() && d.bags.iter().all(|b| b.len() <= 3));
//! assert_eq!(d.edges.len(), d.bags.len() - 1);
//! // K4 has treewidth 3: no decomposition (empty result).
//! let k4: Vec<(usize, usize)> = (0..4).flat_map(|u| (u + 1..4).map(move |v| (u, v))).collect();
//! assert!(tree_decomposition_width2(4, &k4).bags.is_empty());
//! ```

use hash_map::{pair_key, FastMap};

/// A tree decomposition: bag `i` is `bags[i]`, `edges` connect bags (a tree).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Decomposition {
    pub bags: Vec<Vec<usize>>,
    pub edges: Vec<(usize, usize)>,
}

/// A tree decomposition with bags of size at most 3 of the simple graph `edges` on `0..n`
/// (`n >= 1`), or a decomposition with **no bags** when the treewidth exceeds 2. Parallel edges
/// and loops are ignored.
pub fn tree_decomposition_width2(n: usize, edges: &[(usize, usize)]) -> Decomposition {
    assert!(n >= 1, "need at least one vertex");
    let mut present: FastMap<()> = FastMap::with_capacity(edges.len() + n);
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut deg = vec![0usize; n];
    for &(u, v) in edges {
        assert!(u < n && v < n, "edge ({u}, {v}) has an endpoint outside 0..{n}");
        if u != v && !present.contains(pair_key(u, v)) {
            present.insert(pair_key(u, v), ());
            adj[u].push(v);
            adj[v].push(u);
            deg[u] += 1;
            deg[v] += 1;
        }
    }
    let mut alive = vec![true; n];
    let mut queued = vec![false; n];
    let mut queue: Vec<usize> = (0..n).filter(|&v| deg[v] <= 2).collect();
    for &v in &queue {
        queued[v] = true;
    }
    let mut order = Vec::with_capacity(n); // elimination order
    let mut bag_nbrs: Vec<Vec<usize>> = vec![Vec::new(); n]; // neighbours at elimination time
    while let Some(v) = queue.pop() {
        let nbrs: Vec<usize> = adj[v]
            .iter()
            .copied()
            .filter(|&w| alive[w] && present.contains(pair_key(v, w)))
            .collect();
        debug_assert!(nbrs.len() <= 2);
        alive[v] = false;
        for &w in &nbrs {
            present.remove(pair_key(v, w));
            deg[w] -= 1;
        }
        if let [a, b] = nbrs[..] {
            if !present.contains(pair_key(a, b)) {
                present.insert(pair_key(a, b), ());
                adj[a].push(b);
                adj[b].push(a);
                deg[a] += 1;
                deg[b] += 1;
            }
        }
        for &w in &nbrs {
            if deg[w] <= 2 && !queued[w] {
                queued[w] = true;
                queue.push(w);
            }
        }
        bag_nbrs[v] = nbrs;
        order.push(v);
    }
    if order.len() < n {
        return Decomposition { bags: Vec::new(), edges: Vec::new() };
    }
    let mut pos = vec![0usize; n];
    for (i, &v) in order.iter().enumerate() {
        pos[v] = i;
    }
    // bag i belongs to order[i]; parent = bag of the earliest-eliminated neighbour
    let mut bags = Vec::with_capacity(n);
    let mut tree_edges = Vec::with_capacity(n);
    let mut roots = Vec::new();
    for (i, &v) in order.iter().enumerate() {
        let mut bag = vec![v];
        bag.extend_from_slice(&bag_nbrs[v]);
        bags.push(bag);
        match bag_nbrs[v].iter().map(|&w| pos[w]).min() {
            Some(p) => tree_edges.push((i, p)),
            None => roots.push(i),
        }
    }
    for pair in roots.windows(2) {
        tree_edges.push((pair[0], pair[1]));
    }
    Decomposition { bags, edges: tree_edges }
}
