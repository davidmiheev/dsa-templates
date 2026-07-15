//! Heavy-Light Decomposition (HLD).
//!
//! Heavy-Light Decomposition splits a tree into "heavy" paths (chains of
//! vertices where each parent has its largest subtree as its first
//! child) and "light" edges (connecting chains). Once decomposed, any
//! path from `u` to `v` can be split into O(log n) segments, each lying
//! within a single chain. By storing per-chain data in a Fenwick tree
//! or segment tree, you can answer path queries / updates in
//! O(log² n) overall.
//!
//! This crate ships the decomposition itself: it computes for each
//! vertex its parent, depth, subtree size, the head of the chain it
//! belongs to, and its position in the linearised chain order. You can
//! then plug in any range-query data structure of your choice on top
//! of that linearisation.
//!
//! # Example
//!
//! ```
//! use heavy_light_decomposition::Hld;
//!
//! // Tree:
//! //      0
//! //      |
//! //      1
//! //     / \
//! //    2   3
//! //    |   |
//! //    4   5
//! let adj = vec![
//!     vec![1],
//!     vec![0, 2, 3],
//!     vec![1, 4],
//!     vec![1, 5],
//!     vec![2],
//!     vec![3],
//! ];
//! let hld = Hld::new(&adj, 0);
//! assert_eq!(hld.parent[1], Some(0));
//! assert_eq!(hld.depth[4], 3);
//! ```

/// Heavy-Light Decomposition over a rooted tree.
#[derive(Debug, Clone)]
pub struct Hld {
    /// Parent of each vertex in the rooted tree. `None` for the root.
    pub parent: Vec<Option<usize>>,
    /// Depth of each vertex (root has depth 0).
    pub depth: Vec<usize>,
    /// Subtree size of each vertex (including itself).
    pub size: Vec<usize>,
    /// Head of the heavy chain containing each vertex.
    pub head: Vec<usize>,
    /// Position of each vertex in the linearised order of chains.
    pub pos: Vec<usize>,
}

impl Hld {
    /// Build the decomposition for the tree `adj` rooted at `root`.
    ///
    /// `adj[v]` is the list of neighbours of `v`. The tree is assumed
    /// to be connected and undirected; edges to a vertex's parent are
    /// filtered automatically.
    ///
    /// # Panics
    /// Panics if `root >= adj.len()`.
    pub fn new(adj: &[Vec<usize>], root: usize) -> Self {
        let n = adj.len();
        assert!(root < n, "root vertex out of range");

        // Phase 1: parent + depth + subtree size via iterative DFS.
        let mut parent: Vec<Option<usize>> = vec![None; n];
        let mut depth = vec![0usize; n];
        let mut size = vec![1usize; n];
        // We store (vertex, parent, state) where state = 0 means "enter"
        // and state = 1 means "exit".
        let mut stack: Vec<(usize, Option<usize>, u8)> = vec![(root, None, 0)];
        let mut order: Vec<usize> = Vec::with_capacity(n);
        while let Some((v, p, state)) = stack.pop() {
            if state == 0 {
                parent[v] = p;
                depth[v] = match p {
                    Some(pv) => depth[pv] + 1,
                    None => 0,
                };
                stack.push((v, p, 1));
                for &u in &adj[v] {
                    if Some(u) != p {
                        stack.push((u, Some(v), 0));
                    }
                }
            } else {
                order.push(v);
            }
        }
        for &v in &order {
            if v != root {
                if let Some(p) = parent[v] {
                    size[p] += size[v];
                }
            }
        }

        // Phase 2: chain heads + positions via DFS from root, picking
        // the heavy child as the next link.
        let mut head = vec![root; n];
        let mut pos = vec![0usize; n];
        let mut cur_pos = 0usize;
        dfs_decompose(adj, root, usize::MAX, root, &size, &mut head, &mut pos, &mut cur_pos);

        Hld {
            parent,
            depth,
            size,
            head,
            pos,
        }
    }

    /// Decompose the path between `u` and `v` into half-open segments
    /// `[lo, hi)` in the linearised chain order. Pass these segments to
    /// a Fenwick / segment tree built on `pos` to answer path queries.
    pub fn path_segments(&self, u: usize, v: usize) -> Vec<(usize, usize)> {
        let mut segs = Vec::new();
        let (mut u, mut v) = (u, v);
        while self.head[u] != self.head[v] {
            if self.depth[self.head[u]] > self.depth[self.head[v]] {
                let h = self.head[u];
                segs.push((self.pos[h], self.pos[u] + 1));
                u = self.parent[h].unwrap();
            } else {
                let h = self.head[v];
                segs.push((self.pos[h], self.pos[v] + 1));
                v = self.parent[h].unwrap();
            }
        }
        let (lo, hi) = if self.pos[u] <= self.pos[v] {
            (self.pos[u], self.pos[v] + 1)
        } else {
            (self.pos[v], self.pos[u] + 1)
        };
        segs.push((lo, hi));
        segs
    }
}

fn dfs_decompose(
    adj: &[Vec<usize>],
    v: usize,
    p: usize,
    head: usize,
    size: &[usize],
    head_out: &mut [usize],
    pos_out: &mut [usize],
    cur_pos: &mut usize,
) {
    head_out[v] = head;
    pos_out[v] = *cur_pos;
    *cur_pos += 1;

    // Find the heavy child (largest subtree).
    let mut heavy: Option<usize> = None;
    let mut heavy_size = 0usize;
    for &u in &adj[v] {
        if u != p && size[u] > heavy_size {
            heavy_size = size[u];
            heavy = Some(u);
        }
    }

    // Recurse into the heavy child continuing the chain.
    if let Some(h) = heavy {
        dfs_decompose(adj, h, v, head, size, head_out, pos_out, cur_pos);
    }

    // Recurse into light children starting new chains.
    for &u in &adj[v] {
        if u != p && Some(u) != heavy {
            dfs_decompose(adj, u, v, u, size, head_out, pos_out, cur_pos);
        }
    }
}
