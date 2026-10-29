//! Lowest common ancestors and level ancestors of a rooted tree.
//!
//! * [`RootedTree`]: the tree relabelled in BFS order (ancestors get smaller labels and the
//!   parents of consecutive labels never decrease, so every table below is built by sequential
//!   passes). Queries keep using the *original* vertex labels.
//! * [`Lca`]: `O(n)` build, `O(1)` query. A sparse table over blocks of 16 preorder slots plus
//!   in-block prefix / suffix minima; the preorder numbering (`tin`, `size`) is exposed too, so
//!   subtree queries and path-sum Fenwick trees need no separate DFS or HLD.
//! * [`LevelAncestor`]: `O(n log n)` build, `O(1)` k-th ancestor (jump pointers + ladders of the
//!   long-path decomposition), and the k-th vertex of a path.
//! * [`BinaryLifting`]: the textbook `O(log n)` structure, kept for comparison; it is 2-5x slower
//!   than [`Lca`] / [`LevelAncestor`] on deep trees because every step is a cache miss.
//!
//! # Example
//!
//! ```
//! use lca::{Lca, LevelAncestor, RootedTree};
//!
//! // 0 - 1 - {2 - 4, 3 - 5}
//! let t = RootedTree::from_edges(6, &[(0, 1), (1, 2), (1, 3), (2, 4), (3, 5)], 0);
//! let l = Lca::new(&t);
//! assert_eq!(l.lca(4, 5), 1);
//! assert_eq!(l.dist(4, 5), 4);
//! let la = LevelAncestor::new(&t);
//! assert_eq!(la.kth_ancestor(4, 2), 1);
//! assert_eq!(la.path_vertex(&l, 4, 5, 2), 1); // 4, 2, 1, 3, 5
//! ```

/// A rooted tree relabelled in BFS order.
#[derive(Debug, Clone)]
pub struct RootedTree {
    n: usize,
    new_of_old: Vec<u32>,
    old_of_new: Vec<u32>,
    parent: Vec<u32>, // by new label, parent[0] = 0
    depth: Vec<u32>,  // by new label
    depth_old: Vec<u32>,
}

impl RootedTree {
    /// The parent of the root.
    pub const NONE: usize = usize::MAX;

    /// Build from an undirected edge list over `0..n`.
    ///
    /// # Panics
    /// Panics if the edges do not form a tree on `0..n`, `root >= n` or `n >= 2^32`.
    pub fn from_edges(n: usize, edges: &[(usize, usize)], root: usize) -> Self {
        assert!(
            root < n && n < u32::MAX as usize,
            "root out of range or too many vertices"
        );
        assert_eq!(
            edges.len() + 1,
            n,
            "a tree on {n} vertices has {} edges",
            n.wrapping_sub(1)
        );
        let mut start = vec![0u32; n + 2];
        for &(u, v) in edges {
            assert!(u < n && v < n, "edge endpoint out of range");
            start[u + 2] += 1;
            start[v + 2] += 1;
        }
        for i in 2..n + 2 {
            start[i] += start[i - 1];
        }
        let mut to = vec![0u32; 2 * edges.len()];
        for &(u, v) in edges {
            to[start[u + 1] as usize] = v as u32;
            start[u + 1] += 1;
            to[start[v + 1] as usize] = u as u32;
            start[v + 1] += 1;
        }
        let mut new_of_old = vec![u32::MAX; n];
        let mut old_of_new = Vec::with_capacity(n);
        let mut parent = Vec::with_capacity(n);
        let mut depth = Vec::with_capacity(n);
        new_of_old[root] = 0;
        old_of_new.push(root as u32);
        parent.push(0u32);
        depth.push(0u32);
        let mut i = 0;
        while i < old_of_new.len() {
            let v = old_of_new[i] as usize;
            for &u in &to[start[v] as usize..start[v + 1] as usize] {
                if new_of_old[u as usize] == u32::MAX {
                    new_of_old[u as usize] = old_of_new.len() as u32;
                    old_of_new.push(u);
                    parent.push(i as u32);
                    depth.push(depth[i] + 1);
                }
            }
            i += 1;
        }
        assert_eq!(
            old_of_new.len(),
            n,
            "the edges do not form a tree (not connected)"
        );
        Self::finish(n, new_of_old, old_of_new, parent, depth)
    }

    /// Build from a parent array (`parent[root]` is ignored).
    ///
    /// # Panics
    /// Panics if the array does not describe a tree rooted at `root`.
    pub fn from_parents(parent: &[usize], root: usize) -> Self {
        let n = parent.len();
        assert!(
            root < n && n < u32::MAX as usize,
            "root out of range or too many vertices"
        );
        let edges: Vec<(usize, usize)> = (0..n)
            .filter(|&v| v != root)
            .map(|v| (v, parent[v]))
            .collect();
        Self::from_edges(n, &edges, root)
    }

    /// Build from a parent array in which `parent[v] < v` for `v > 0` (root 0): the labels are
    /// already a valid BFS-like order, so nothing is relabelled and no search is needed.
    ///
    /// # Panics
    /// Panics if some `parent[v] >= v`.
    pub fn from_increasing_parents(parent: &[usize]) -> Self {
        let n = parent.len();
        assert!(
            n >= 1 && n < u32::MAX as usize,
            "need 1 <= n < 2^32 vertices"
        );
        let mut p = vec![0u32; n];
        let mut depth = vec![0u32; n];
        for v in 1..n {
            assert!(
                parent[v] < v,
                "parent[{v}] = {} is not smaller than {v}",
                parent[v]
            );
            p[v] = parent[v] as u32;
            depth[v] = depth[parent[v]] + 1;
        }
        let id: Vec<u32> = (0..n as u32).collect();
        RootedTree {
            n,
            new_of_old: id.clone(),
            old_of_new: id,
            parent: p,
            depth: depth.clone(),
            depth_old: depth,
        }
    }

    fn finish(
        n: usize,
        new_of_old: Vec<u32>,
        old_of_new: Vec<u32>,
        parent: Vec<u32>,
        depth: Vec<u32>,
    ) -> Self {
        let mut depth_old = vec![0u32; n];
        for (new, &old) in old_of_new.iter().enumerate() {
            depth_old[old as usize] = depth[new];
        }
        RootedTree {
            n,
            new_of_old,
            old_of_new,
            parent,
            depth,
            depth_old,
        }
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Depth of `v` (the root has depth 0).
    pub fn depth(&self, v: usize) -> usize {
        self.depth_old[v] as usize
    }

    /// Parent of `v`, [`RootedTree::NONE`] for the root.
    pub fn parent(&self, v: usize) -> usize {
        let new = self.new_of_old[v] as usize;
        if new == 0 {
            Self::NONE
        } else {
            self.old_of_new[self.parent[new] as usize] as usize
        }
    }

    fn max_depth(&self) -> usize {
        self.depth.iter().copied().max().unwrap_or(0) as usize
    }

    fn levels(&self) -> usize {
        (usize::BITS - self.max_depth().leading_zeros()).max(1) as usize
    }

    /// Jump-pointer table `up[k * n + v] = 2^k-th ancestor` over BFS labels.
    fn jump_table(&self) -> (usize, Vec<u32>) {
        let (n, levels) = (self.n, self.levels());
        let mut up = vec![0u32; levels * n];
        up[..n].copy_from_slice(&self.parent);
        for k in 1..levels {
            let (lo, hi) = up.split_at_mut(k * n);
            let prev = &lo[(k - 1) * n..];
            for v in 0..n {
                hi[v] = prev[prev[v] as usize];
            }
        }
        (levels, up)
    }
}

const BLOCK_SHIFT: usize = 4;

/// `O(1)` lowest common ancestor with a preorder numbering.
#[derive(Debug, Clone)]
pub struct Lca {
    tin_of_old: Vec<u32>,
    size_old: Vec<u32>,
    depth_old: Vec<u32>,
    old_by_tin: Vec<u32>,
    key: Vec<u32>,
    prefix: Vec<u32>,
    suffix: Vec<u32>,
    table: Vec<u32>,
    blocks: usize,
}

impl Lca {
    pub fn new(t: &RootedTree) -> Self {
        let n = t.n;
        // subtree sizes and preorder numbers by sequential passes over the BFS labels: the
        // children of a vertex get consecutive blocks right after it
        let mut size = vec![1u32; n];
        for v in (1..n).rev() {
            size[t.parent[v] as usize] += size[v];
        }
        let mut tin = vec![0u32; n];
        let mut next_free = vec![0u32; n];
        next_free[0] = 1;
        for v in 1..n {
            let p = t.parent[v] as usize;
            tin[v] = next_free[p];
            next_free[p] += size[v];
            next_free[v] = tin[v] + 1;
        }
        let mut key = vec![0u32; n]; // key[tin[v]] = tin of the parent of v
        let mut old_by_tin = vec![0u32; n];
        let mut tin_of_old = vec![0u32; n];
        let mut size_old = vec![0u32; n];
        for v in 0..n {
            let old = t.old_of_new[v] as usize;
            old_by_tin[tin[v] as usize] = old as u32;
            tin_of_old[old] = tin[v];
            size_old[old] = size[v];
            if v > 0 {
                key[tin[v] as usize] = tin[t.parent[v] as usize];
            }
        }
        let b = 1usize << BLOCK_SHIFT;
        let blocks = n.div_ceil(b);
        let mut prefix = key.clone();
        let mut suffix = key.clone();
        for i in 0..n {
            if i & (b - 1) != 0 {
                prefix[i] = prefix[i].min(prefix[i - 1]);
            }
        }
        for i in (0..n.saturating_sub(1)).rev() {
            if i & (b - 1) != b - 1 {
                suffix[i] = suffix[i].min(suffix[i + 1]);
            }
        }
        let levels = (usize::BITS - blocks.leading_zeros()) as usize;
        let mut table = vec![0u32; levels * blocks];
        for j in 0..blocks {
            table[j] = suffix[j * b];
        }
        for k in 1..levels {
            let half = 1usize << (k - 1);
            let (lo, hi) = table.split_at_mut(k * blocks);
            let prev = &lo[(k - 1) * blocks..];
            for i in 0..(blocks + 1).saturating_sub(1 << k) {
                hi[i] = prev[i].min(prev[i + half]);
            }
        }
        Lca {
            tin_of_old,
            size_old,
            depth_old: t.depth_old.clone(),
            old_by_tin,
            key,
            prefix,
            suffix,
            table,
            blocks,
        }
    }

    /// The lowest common ancestor of `u` and `v`.
    #[inline]
    pub fn lca(&self, u: usize, v: usize) -> usize {
        if u == v {
            return u;
        }
        let (mut l, mut r) = (self.tin_of_old[u] as usize, self.tin_of_old[v] as usize);
        if l > r {
            std::mem::swap(&mut l, &mut r);
        }
        l += 1; // the parents of the vertices at preorder (tin[u], tin[v]]
        let (bl, br) = (l >> BLOCK_SHIFT, r >> BLOCK_SHIFT);
        let m = if bl == br {
            self.key[l..=r]
                .iter()
                .copied()
                .min()
                .expect("non-empty range")
        } else {
            let mut m = self.suffix[l].min(self.prefix[r]);
            if bl + 1 < br {
                let (a, b) = (bl + 1, br - 1);
                let k = (usize::BITS - 1 - (b - a + 1).leading_zeros()) as usize;
                m = m
                    .min(self.table[k * self.blocks + a])
                    .min(self.table[k * self.blocks + b + 1 - (1 << k)]);
            }
            m
        };
        self.old_by_tin[m as usize] as usize
    }

    /// Number of edges on the path between `u` and `v`.
    pub fn dist(&self, u: usize, v: usize) -> usize {
        let l = self.lca(u, v);
        (self.depth_old[u] + self.depth_old[v] - 2 * self.depth_old[l]) as usize
    }

    /// Depth of `v`.
    pub fn depth(&self, v: usize) -> usize {
        self.depth_old[v] as usize
    }

    /// Preorder number of `v` (the root has 0; a subtree is a contiguous range).
    pub fn tin(&self, v: usize) -> usize {
        self.tin_of_old[v] as usize
    }

    /// Number of vertices in the subtree of `v`.
    pub fn size(&self, v: usize) -> usize {
        self.size_old[v] as usize
    }

    /// The preorder range `[l, r)` of the subtree of `v`.
    pub fn subtree(&self, v: usize) -> (usize, usize) {
        let l = self.tin(v);
        (l, l + self.size(v))
    }

    /// Whether `u` is an ancestor of `v` (every vertex is its own ancestor).
    pub fn is_ancestor(&self, u: usize, v: usize) -> bool {
        let (l, r) = self.subtree(u);
        (l..r).contains(&self.tin(v))
    }

    /// The vertex with preorder number `i`.
    pub fn vertex_at(&self, i: usize) -> usize {
        self.old_by_tin[i] as usize
    }
}

/// `O(1)` k-th ancestor (jump pointers + ladders).
#[derive(Debug, Clone)]
pub struct LevelAncestor {
    n: usize,
    up: Vec<u32>,
    slot: Vec<u32>,   // by BFS label: index of the vertex inside `ladder`
    ladder: Vec<u32>, // original labels
    new_of_old: Vec<u32>,
    depth_old: Vec<u32>,
}

impl LevelAncestor {
    pub fn new(t: &RootedTree) -> Self {
        let n = t.n;
        let (_, up) = t.jump_table();
        // long-path decomposition: each vertex continues into its tallest child
        let mut height = vec![0u32; n];
        let mut long = vec![u32::MAX; n];
        for v in (1..n).rev() {
            let p = t.parent[v] as usize;
            if long[p] == u32::MAX || height[v] + 1 > height[p] {
                height[p] = height[v] + 1;
                long[p] = v as u32;
            }
        }
        // a path of length L is stored with up to L ancestors above its top: any vertex whose
        // 2^h-th descendant exists then finds its remaining k < 2^h ancestors inside its ladder
        let mut slot = vec![0u32; n];
        let mut ladder: Vec<u32> = Vec::with_capacity(2 * n);
        for top in 0..n {
            if top != 0 && long[t.parent[top] as usize] == top as u32 {
                continue;
            }
            let len = height[top] as usize + 1;
            let ext = len.min(t.depth[top] as usize);
            let base = ladder.len();
            ladder.resize(base + ext, 0);
            let mut a = top;
            for i in 0..ext {
                a = t.parent[a] as usize;
                ladder[base + ext - 1 - i] = t.old_of_new[a];
            }
            let mut x = top;
            loop {
                slot[x] = ladder.len() as u32;
                ladder.push(t.old_of_new[x]);
                if long[x] == u32::MAX {
                    break;
                }
                x = long[x] as usize;
            }
        }
        LevelAncestor {
            n,
            up,
            slot,
            ladder,
            new_of_old: t.new_of_old.clone(),
            depth_old: t.depth_old.clone(),
        }
    }

    /// The `k`-th ancestor of `v` (`k = 0` is `v`).
    ///
    /// # Panics
    /// Panics if `k` exceeds the depth of `v`.
    #[inline]
    pub fn kth_ancestor(&self, v: usize, k: usize) -> usize {
        assert!(
            k <= self.depth_old[v] as usize,
            "vertex {v} has no ancestor {k} levels up"
        );
        if k == 0 {
            return v;
        }
        let h = (usize::BITS - 1 - k.leading_zeros()) as usize;
        let w = self.up[h * self.n + self.new_of_old[v] as usize] as usize;
        self.ladder[self.slot[w] as usize - (k - (1 << h))] as usize
    }

    /// The `k`-th vertex on the path from `u` to `v` (`k = 0` is `u`, `k = dist(u, v)` is `v`).
    ///
    /// # Panics
    /// Panics if `k > dist(u, v)`.
    pub fn path_vertex(&self, lca: &Lca, u: usize, v: usize, k: usize) -> usize {
        let l = lca.lca(u, v);
        let up = lca.depth(u) - lca.depth(l);
        let down = lca.depth(v) - lca.depth(l);
        assert!(
            k <= up + down,
            "the path has only {} vertices",
            up + down + 1
        );
        if k <= up {
            self.kth_ancestor(u, k)
        } else {
            self.kth_ancestor(v, up + down - k)
        }
    }
}

/// Textbook binary lifting, `O(n log n)` memory and `O(log n)` per query.
#[derive(Debug, Clone)]
pub struct BinaryLifting {
    n: usize,
    levels: usize,
    up: Vec<u32>,
    depth: Vec<u32>, // by BFS label
    new_of_old: Vec<u32>,
    old_of_new: Vec<u32>,
}

impl BinaryLifting {
    pub fn new(t: &RootedTree) -> Self {
        let (levels, up) = t.jump_table();
        BinaryLifting {
            n: t.n,
            levels,
            up,
            depth: t.depth.clone(),
            new_of_old: t.new_of_old.clone(),
            old_of_new: t.old_of_new.clone(),
        }
    }

    fn kth(&self, mut v: usize, mut k: usize) -> usize {
        let mut level = 0;
        while k > 0 {
            if k & 1 == 1 {
                v = self.up[level * self.n + v] as usize;
            }
            k >>= 1;
            level += 1;
        }
        v
    }

    /// The `k`-th ancestor of `v`.
    ///
    /// # Panics
    /// Panics if `k` exceeds the depth of `v`.
    pub fn kth_ancestor(&self, v: usize, k: usize) -> usize {
        let new = self.new_of_old[v] as usize;
        assert!(
            k <= self.depth[new] as usize,
            "vertex {v} has no ancestor {k} levels up"
        );
        self.old_of_new[self.kth(new, k)] as usize
    }

    /// The lowest common ancestor of `u` and `v`.
    pub fn lca(&self, u: usize, v: usize) -> usize {
        let (mut u, mut v) = (self.new_of_old[u] as usize, self.new_of_old[v] as usize);
        if self.depth[u] < self.depth[v] {
            std::mem::swap(&mut u, &mut v);
        }
        u = self.kth(u, (self.depth[u] - self.depth[v]) as usize);
        if u != v {
            for k in (0..self.levels).rev() {
                let (a, b) = (
                    self.up[k * self.n + u] as usize,
                    self.up[k * self.n + v] as usize,
                );
                if a != b {
                    u = a;
                    v = b;
                }
            }
            u = self.up[u] as usize;
        }
        self.old_of_new[u] as usize
    }
}
