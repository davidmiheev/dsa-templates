//! Heavy-Light Decomposition (HLD) of a rooted tree.
//!
//! The tree is cut into *heavy chains* (every vertex continues the chain into its child with the
//! largest subtree) and the vertices are numbered in a DFS preorder that visits the heavy child
//! first. Consequently
//!
//! * every chain occupies consecutive positions,
//! * the subtree of `v` is the position range `pos[v] .. pos[v] + size[v]`,
//! * any path is a union of `O(log n)` chain pieces, i.e. `O(log n)` position ranges.
//!
//! Put a Fenwick / segment tree over the positions and path or subtree queries become range
//! queries. On top of that the crate answers [`lca`](Hld::lca), [`dist`](Hld::dist),
//! [`jump`](Hld::jump) (k-th ancestor) and [`path_vertex`](Hld::path_vertex) (k-th vertex on a
//! path) in `O(log n)`.
//!
//! Everything is iterative, so million-vertex paths are fine. The root's parent is
//! [`Hld::NONE`]. Edges can be given as adjacency lists, an edge list, or a parent array.
//!
//! # Example
//!
//! ```
//! use heavy_light_decomposition::Hld;
//!
//! // 0 - 1 - {2 - 4, 3 - 5}
//! let hld = Hld::from_edges(6, &[(0, 1), (1, 2), (1, 3), (2, 4), (3, 5)], 0);
//! assert_eq!(hld.parent[1], 0);
//! assert_eq!(hld.parent[0], Hld::NONE);
//! assert_eq!(hld.lca(4, 5), 1);
//! assert_eq!(hld.dist(4, 5), 4);
//! assert_eq!(hld.path_vertex(4, 5, 2), 1); // 4, 2, 1, 3, 5
//! // The path 4 -> 5 as position ranges (one per chain it touches).
//! assert_eq!(hld.path_segments(4, 5).len(), 2);
//! ```

/// Heavy-light decomposition over a tree rooted at `root`.
#[derive(Debug, Clone)]
pub struct Hld {
    /// Parent of each vertex; [`Hld::NONE`] for the root.
    pub parent: Vec<usize>,
    /// Depth of each vertex (the root has depth 0).
    pub depth: Vec<usize>,
    /// Subtree size of each vertex (including itself).
    pub size: Vec<usize>,
    /// Head (top vertex) of the heavy chain containing each vertex.
    pub head: Vec<usize>,
    /// Position of each vertex in the DFS preorder (heavy child first).
    pub pos: Vec<usize>,
    /// Inverse of `pos`: the vertex at each position.
    pub order: Vec<usize>,
    /// By position: `(position of the chain head, position of the head's parent)`, packed so
    /// that one hop of a chain walk touches a single cache line. The root chain's parent is
    /// `u32::MAX`.
    chain: Vec<(u32, u32)>,
}

impl Hld {
    /// The parent of the root.
    pub const NONE: usize = usize::MAX;

    /// Build from adjacency lists (`adj[v]` = neighbours of `v`).
    ///
    /// # Panics
    /// Panics if `root >= adj.len()`, or the graph is not a tree.
    pub fn new(adj: &[Vec<usize>], root: usize) -> Self {
        let n = adj.len();
        let mut offsets = vec![0usize; n + 1];
        for (v, a) in adj.iter().enumerate() {
            offsets[v + 1] = offsets[v] + a.len();
        }
        let targets: Vec<usize> = adj.iter().flatten().copied().collect();
        Self::from_csr(&offsets, &targets, root)
    }

    /// Build from an undirected edge list over vertices `0..n`.
    ///
    /// # Panics
    /// Panics if `root >= n`, an endpoint is out of range, or the edges do not form a tree.
    pub fn from_edges(n: usize, edges: &[(usize, usize)], root: usize) -> Self {
        assert_eq!(edges.len() + 1, n, "a tree on {n} vertices has {} edges", n.wrapping_sub(1));
        let mut offsets = vec![0usize; n + 2];
        for &(u, v) in edges {
            assert!(u < n && v < n, "edge endpoint out of range");
            offsets[u + 2] += 1;
            offsets[v + 2] += 1;
        }
        for i in 2..n + 2 {
            offsets[i] += offsets[i - 1];
        }
        let mut targets = vec![0usize; 2 * edges.len()];
        for &(u, v) in edges {
            targets[offsets[u + 1]] = v;
            offsets[u + 1] += 1;
            targets[offsets[v + 1]] = u;
            offsets[v + 1] += 1;
        }
        offsets.pop();
        Self::from_csr(&offsets, &targets, root)
    }

    /// Build from a parent array: `parent[v]` is the parent of `v`, and the entry of `root` is
    /// ignored.
    ///
    /// # Panics
    /// Panics if `root` is out of range or the array does not describe a tree rooted at `root`.
    pub fn from_parents(parent: &[usize], root: usize) -> Self {
        let n = parent.len();
        assert!(root < n, "root vertex out of range");
        let mut offsets = vec![0usize; n + 2];
        for (v, &p) in parent.iter().enumerate() {
            if v != root {
                assert!(p < n, "parent out of range");
                offsets[p + 2] += 1;
            }
        }
        for i in 2..n + 2 {
            offsets[i] += offsets[i - 1];
        }
        let mut children = vec![0usize; n.saturating_sub(1)];
        for (v, &p) in parent.iter().enumerate() {
            if v != root {
                children[offsets[p + 1]] = v;
                offsets[p + 1] += 1;
            }
        }
        offsets.pop();
        Self::from_rooted(parent, &offsets, &children, root)
    }

    fn from_csr(offsets: &[usize], targets: &[usize], root: usize) -> Self {
        let n = offsets.len() - 1;
        assert!(root < n, "root vertex out of range");
        // BFS to find parents; the tree is checked to be connected and acyclic.
        let mut parent = vec![Self::NONE; n];
        let mut order = Vec::with_capacity(n);
        let mut seen = vec![false; n];
        seen[root] = true;
        order.push(root);
        let mut i = 0;
        while i < order.len() {
            let v = order[i];
            i += 1;
            for &u in &targets[offsets[v]..offsets[v + 1]] {
                if u != parent[v] {
                    assert!(!seen[u], "the graph contains a cycle");
                    seen[u] = true;
                    parent[u] = v;
                    order.push(u);
                }
            }
        }
        assert_eq!(order.len(), n, "the graph is not connected");
        // Children lists: neighbours except the parent, in the same CSR slots.
        let mut child_offsets = vec![0usize; n + 1];
        for v in 0..n {
            child_offsets[v + 1] = child_offsets[v] + (offsets[v + 1] - offsets[v]) - (v != root) as usize;
        }
        let mut children = vec![0usize; n.saturating_sub(1)];
        for v in 0..n {
            let mut k = child_offsets[v];
            for &u in &targets[offsets[v]..offsets[v + 1]] {
                if u != parent[v] {
                    children[k] = u;
                    k += 1;
                }
            }
        }
        Self::from_rooted(&parent, &child_offsets, &children, root)
    }

    fn from_rooted(parent: &[usize], offsets: &[usize], children: &[usize], root: usize) -> Self {
        let n = parent.len();
        // BFS order, depth.
        let mut bfs = Vec::with_capacity(n);
        let mut depth = vec![0usize; n];
        bfs.push(root);
        let mut i = 0;
        while i < bfs.len() {
            let v = bfs[i];
            i += 1;
            for &u in &children[offsets[v]..offsets[v + 1]] {
                depth[u] = depth[v] + 1;
                bfs.push(u);
            }
        }
        assert_eq!(bfs.len(), n, "the parent array does not describe a tree rooted at {root}");
        let mut size = vec![1usize; n];
        for &v in bfs.iter().skip(1).rev() {
            size[parent[v]] += size[v];
        }
        let mut real_parent = parent.to_vec();
        real_parent[root] = Self::NONE;
        // Preorder positions with the heavy child first, assigned top-down in BFS order: the
        // heavy child directly follows its parent, the light children follow in turn, each
        // after the whole subtree of the previous one.
        assert!(n < u32::MAX as usize, "too many vertices");
        let mut head = vec![root; n];
        let mut pos = vec![0usize; n];
        let mut chain = vec![(0u32, u32::MAX); n];
        for &v in &bfs {
            let kids = &children[offsets[v]..offsets[v + 1]];
            let Some(heavy) = kids.iter().copied().max_by_key(|&u| size[u]) else { continue };
            pos[heavy] = pos[v] + 1;
            head[heavy] = head[v];
            chain[pos[heavy]] = chain[pos[v]];
            let mut next = pos[v] + 1 + size[heavy];
            for &u in kids {
                if u != heavy {
                    pos[u] = next;
                    head[u] = u;
                    chain[next] = (next as u32, pos[v] as u32);
                    next += size[u];
                }
            }
        }
        let mut order = vec![0usize; n];
        for (v, &p) in pos.iter().enumerate() {
            order[p] = v;
        }
        Hld { parent: real_parent, depth, size, head, pos, order, chain }
    }

    /// Number of vertices.
    pub fn len(&self) -> usize {
        self.parent.len()
    }

    pub fn is_empty(&self) -> bool {
        self.parent.is_empty()
    }

    /// The positions `pos[v] .. pos[v] + size[v]` covered by the subtree of `v`.
    pub fn subtree(&self, v: usize) -> (usize, usize) {
        (self.pos[v], self.pos[v] + self.size[v])
    }

    /// Whether `u` is an ancestor of `v` (every vertex is its own ancestor).
    pub fn is_ancestor(&self, u: usize, v: usize) -> bool {
        self.pos[u] <= self.pos[v] && self.pos[v] < self.pos[u] + self.size[u]
    }

    /// Lowest common ancestor of `u` and `v`.
    pub fn lca(&self, u: usize, v: usize) -> usize {
        // Work in position space: ancestors have smaller positions, and the chain whose head
        // lies later in the preorder cannot contain the LCA above its head.
        let (mut a, mut b) = (self.pos[u], self.pos[v]);
        loop {
            let (ha, hb) = (self.chain[a].0, self.chain[b].0);
            if ha == hb {
                return self.order[a.min(b)];
            }
            if ha > hb {
                a = self.chain[a].1 as usize;
            } else {
                b = self.chain[b].1 as usize;
            }
        }
    }

    /// Number of edges on the path between `u` and `v`.
    pub fn dist(&self, u: usize, v: usize) -> usize {
        self.depth[u] + self.depth[v] - 2 * self.depth[self.lca(u, v)]
    }

    /// The `k`-th ancestor of `v` (`k = 0` is `v` itself).
    ///
    /// # Panics
    /// Panics if `k > depth[v]`.
    pub fn jump(&self, mut v: usize, mut k: usize) -> usize {
        assert!(k <= self.depth[v], "vertex {v} has no ancestor {k} levels up");
        loop {
            let h = self.head[v];
            let up = self.depth[v] - self.depth[h];
            if k <= up {
                return self.order[self.pos[v] - k];
            }
            k -= up + 1;
            v = self.parent[h];
        }
    }

    /// The `k`-th vertex on the path from `u` to `v` (`k = 0` is `u`, `k = dist(u, v)` is `v`).
    ///
    /// # Panics
    /// Panics if `k > dist(u, v)`.
    pub fn path_vertex(&self, u: usize, v: usize, k: usize) -> usize {
        let l = self.lca(u, v);
        let up = self.depth[u] - self.depth[l];
        let down = self.depth[v] - self.depth[l];
        assert!(k <= up + down, "the path has only {} vertices", up + down + 1);
        if k <= up {
            self.jump(u, k)
        } else {
            self.jump(v, up + down - k)
        }
    }

    /// Calls `f(lo, hi)` for each of the `O(log n)` position ranges `[lo, hi)` that together
    /// make up the path between `u` and `v`, in no particular order and without allocating.
    /// Suitable for commutative aggregates.
    pub fn for_each_path_segment(&self, u: usize, v: usize, mut f: impl FnMut(usize, usize)) {
        let (mut a, mut b) = (self.pos[u], self.pos[v]);
        loop {
            let (ha, hb) = (self.chain[a].0 as usize, self.chain[b].0 as usize);
            if ha == hb {
                f(a.min(b), a.max(b) + 1);
                return;
            }
            if ha > hb {
                f(ha, a + 1);
                a = self.chain[a].1 as usize;
            } else {
                f(hb, b + 1);
                b = self.chain[b].1 as usize;
            }
        }
    }

    /// The vertices of the path between `u` and `v` as position ranges `[lo, hi)`, in no
    /// particular order (see [`for_each_path_segment`](Self::for_each_path_segment)).
    pub fn path_segments(&self, u: usize, v: usize) -> Vec<(usize, usize)> {
        let mut segs = Vec::new();
        self.for_each_path_segment(u, v, |lo, hi| segs.push((lo, hi)));
        segs
    }

    /// Calls `f(lo, hi, reversed)` for the position ranges of the path from `u` to `v` in
    /// *walking order*, without allocating; `reversed` is `true` when the range is walked from
    /// `hi - 1` down to `lo` (towards the root). Use this for non-commutative aggregates: fold
    /// the range forwards or backwards accordingly.
    pub fn for_each_path_segment_ordered(&self, u: usize, v: usize, mut f: impl FnMut(usize, usize, bool)) {
        let (mut a, mut b) = (self.pos[u], self.pos[v]);
        // Pieces on the `v` side are discovered bottom-up but walked top-down: buffer them (at
        // most one per light edge, so well below 64).
        let mut back = [(0usize, 0usize); 64];
        let mut count = 0;
        loop {
            let (ha, hb) = (self.chain[a].0 as usize, self.chain[b].0 as usize);
            if ha == hb {
                if a >= b {
                    f(b, a + 1, true);
                } else {
                    f(a, b + 1, false);
                }
                break;
            }
            if ha > hb {
                f(ha, a + 1, true);
                a = self.chain[a].1 as usize;
            } else {
                back[count] = (hb, b + 1);
                count += 1;
                b = self.chain[b].1 as usize;
            }
        }
        for &(lo, hi) in back[..count].iter().rev() {
            f(lo, hi, false);
        }
    }

    /// The path from `u` to `v` as position ranges in walking order, see
    /// [`for_each_path_segment_ordered`](Self::for_each_path_segment_ordered).
    pub fn path_ordered(&self, u: usize, v: usize) -> Vec<(usize, usize, bool)> {
        let mut segs = Vec::new();
        self.for_each_path_segment_ordered(u, v, |lo, hi, rev| segs.push((lo, hi, rev)));
        segs
    }
}
