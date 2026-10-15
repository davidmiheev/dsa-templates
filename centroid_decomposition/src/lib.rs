//! Centroid decomposition of a tree, with distances to every centroid ancestor, and two
//! ready-made structures on top of it:
//!
//! * [`ContourSum`]: add to a vertex value, and sum the values of all vertices whose distance
//!   from a given vertex lies in `[l, r)`;
//! * [`ContourAdd`]: add `x` to all vertices whose distance from a given vertex lies in
//!   `[l, r)`, and read a single vertex value.
//!
//! Recursively splitting a tree at a centroid (a vertex whose removal leaves pieces of at most
//! half the size) gives the *centroid tree* of height `O(log n)`. Every vertex has `O(log n)`
//! centroid ancestors, and the path between two vertices passes through their lowest common
//! centroid ancestor, which is the basis of most "count / sum over paths" divide and conquer
//! on trees. Construction is iterative and runs in `O(n log n)`; `n * (height + 1)` distances
//! are stored (as `u32`).
//!
//! The centroid-tree root has parent [`CentroidTree::NONE`].
//!
//! # Example
//!
//! ```
//! use centroid_decomposition::{CentroidTree, ContourSum};
//!
//! // Path 0 - 1 - 2 - 3 - 4 with values 1, 10, 100, 1000, 10000.
//! let edges = [(0, 1), (1, 2), (2, 3), (3, 4)];
//! let tree = CentroidTree::from_edges(5, &edges);
//! assert_eq!(tree.parent[2], CentroidTree::NONE); // the middle vertex is the root
//! let mut c = ContourSum::<i64>::new(tree);
//! for (v, x) in [1, 10, 100, 1000, 10000].into_iter().enumerate() {
//!     c.add(v, x);
//! }
//! // Vertices at distance 1 or 2 from vertex 1: 0, 2 (distance 1) and 3 (distance 2).
//! assert_eq!(c.sum(1, 1, 3), 1 + 100 + 1000);
//! ```

use std::ops::{AddAssign, Sub};

/// A centroid decomposition. See the [module documentation](self).
#[derive(Debug, Clone)]
pub struct CentroidTree {
    /// Parent in the centroid tree ([`CentroidTree::NONE`] for the root).
    pub parent: Vec<usize>,
    /// Depth in the centroid tree (the root has depth 0).
    pub depth: Vec<usize>,
    /// For a centroid `c`: the largest distance from `c` to a vertex of its component.
    pub radius: Vec<usize>,
    /// The root of the centroid tree.
    pub root: usize,
    /// `dist[level][v]`: distance from `v` to its centroid ancestor at depth `level`
    /// (valid for `level <= depth[v]`).
    dist: Vec<Vec<u32>>,
}

impl CentroidTree {
    /// The parent of the root.
    pub const NONE: usize = usize::MAX;

    /// Decompose the tree given as adjacency lists.
    ///
    /// # Panics
    /// Panics if `adj` is empty or not a tree.
    pub fn new(adj: &[Vec<usize>]) -> Self {
        let n = adj.len();
        let mut offsets = vec![0usize; n + 1];
        for (v, a) in adj.iter().enumerate() {
            offsets[v + 1] = offsets[v] + a.len();
        }
        let targets: Vec<usize> = adj.iter().flatten().copied().collect();
        Self::from_csr(&offsets, &targets)
    }

    /// Decompose the tree on `0..n` with the given undirected edges.
    ///
    /// # Panics
    /// Panics if `n == 0` or the edges do not form a tree.
    pub fn from_edges(n: usize, edges: &[(usize, usize)]) -> Self {
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
        Self::from_csr(&offsets, &targets)
    }

    fn from_csr(offsets: &[usize], targets: &[usize]) -> Self {
        let n = offsets.len() - 1;
        assert!(n > 0, "the tree must have at least one vertex");
        let none = Self::NONE;
        let mut parent = vec![none; n];
        let mut depth = vec![0usize; n];
        let mut radius = vec![0usize; n];
        let mut dist: Vec<Vec<u32>> = Vec::new();
        let mut removed = vec![false; n];
        let mut size = vec![0usize; n];
        let mut bfs_parent = vec![none; n];
        let mut queue: Vec<usize> = Vec::with_capacity(n);
        let mut root = none;
        let mut decomposed = 0;
        // Components still to split: (entry vertex, centroid parent, level).
        let mut todo: Vec<(usize, usize, usize)> = vec![(0, none, 0)];
        while let Some((entry, centroid_parent, level)) = todo.pop() {
            // Sizes of the component, rooted at `entry`.
            queue.clear();
            queue.push(entry);
            bfs_parent[entry] = none;
            let mut i = 0;
            while i < queue.len() {
                let v = queue[i];
                i += 1;
                size[v] = 1;
                for &u in &targets[offsets[v]..offsets[v + 1]] {
                    if !removed[u] && u != bfs_parent[v] {
                        assert!(bfs_parent[u] == none && u != entry, "the graph contains a cycle");
                        bfs_parent[u] = v;
                        queue.push(u);
                    }
                }
            }
            let total = queue.len();
            for &v in queue.iter().skip(1).rev() {
                size[bfs_parent[v]] += size[v];
            }
            // Walk towards the heavy child until no child is larger than half.
            let mut c = entry;
            loop {
                let heavy = targets[offsets[c]..offsets[c + 1]].iter().copied().find(|&u| !removed[u] && bfs_parent[u] == c && size[u] * 2 > total);
                match heavy {
                    Some(u) => c = u,
                    None => break,
                }
            }
            for &v in &queue {
                bfs_parent[v] = none; // reset for later components
            }
            parent[c] = centroid_parent;
            depth[c] = level;
            if centroid_parent == none {
                root = c;
            }
            decomposed += 1;
            // Distances from `c` within the component.
            if dist.len() == level {
                dist.push(vec![0u32; n]);
            }
            let d = &mut dist[level];
            queue.clear();
            queue.push(c);
            d[c] = 0;
            bfs_parent[c] = none;
            let mut i = 0;
            let mut far = 0;
            while i < queue.len() {
                let v = queue[i];
                i += 1;
                far = far.max(d[v] as usize);
                for &u in &targets[offsets[v]..offsets[v + 1]] {
                    if !removed[u] && u != bfs_parent[v] {
                        bfs_parent[u] = v;
                        d[u] = d[v] + 1;
                        queue.push(u);
                    }
                }
            }
            for &v in &queue {
                bfs_parent[v] = none;
            }
            radius[c] = far;
            removed[c] = true;
            for &u in &targets[offsets[c]..offsets[c + 1]] {
                if !removed[u] {
                    todo.push((u, c, level + 1));
                }
            }
        }
        assert_eq!(decomposed, n, "the graph is not connected");
        CentroidTree { parent, depth, radius, root, dist }
    }

    /// Number of vertices.
    pub fn len(&self) -> usize {
        self.parent.len()
    }

    pub fn is_empty(&self) -> bool {
        self.parent.is_empty()
    }

    /// The distance from `v` to its centroid ancestor at depth `level` (`level == depth[v]`
    /// is `v` itself, at distance 0).
    ///
    /// # Panics
    /// Panics if `level > depth[v]`.
    pub fn dist_to_ancestor(&self, v: usize, level: usize) -> usize {
        assert!(level <= self.depth[v], "vertex {v} has no centroid ancestor at depth {level}");
        self.dist[level][v] as usize
    }

    /// The centroid ancestors of `v`, from `v` itself up to the root, as `(centroid, distance
    /// from v to it)`.
    pub fn ancestors(&self, v: usize) -> impl Iterator<Item = (usize, usize)> + '_ {
        let mut cur = v;
        let mut level = self.depth[v] as isize;
        std::iter::from_fn(move || {
            if level < 0 {
                return None;
            }
            let item = (cur, self.dist[level as usize][v] as usize);
            cur = self.parent[cur];
            level -= 1;
            Some(item)
        })
    }
}

/// Fenwick trees stored back to back: tree `t` has `len(t)` slots.
#[derive(Debug, Clone)]
struct FlatFenwick<T> {
    start: Vec<usize>,
    data: Vec<T>,
}

impl<T: Copy + Default + AddAssign + Sub<Output = T>> FlatFenwick<T> {
    fn new(lens: impl Iterator<Item = usize>) -> Self {
        let mut start = vec![0usize];
        for l in lens {
            start.push(start.last().unwrap() + l);
        }
        FlatFenwick { data: vec![T::default(); *start.last().unwrap()], start }
    }

    fn len(&self, t: usize) -> usize {
        self.start[t + 1] - self.start[t]
    }

    /// Add `x` at slot `i` of tree `t` (ignored if `i` is past the end).
    fn add(&mut self, t: usize, i: usize, x: T) {
        let (base, n) = (self.start[t], self.len(t));
        let mut i = i + 1;
        while i <= n {
            self.data[base + i - 1] += x;
            i += i & i.wrapping_neg();
        }
    }

    /// Sum of the first `k` slots of tree `t` (`k` is clamped to the length).
    fn prefix(&self, t: usize, k: usize) -> T {
        let base = self.start[t];
        let mut k = k.min(self.len(t));
        let mut s = T::default();
        while k > 0 {
            s += self.data[base + k - 1];
            k &= k - 1;
        }
        s
    }
}

/// Point add / sum of the values within a distance range, on a fixed tree.
///
/// ```
/// use centroid_decomposition::{CentroidTree, ContourSum};
/// let tree = CentroidTree::from_edges(3, &[(0, 1), (1, 2)]);
/// let mut c = ContourSum::<i64>::new(tree);
/// c.add(2, 5);
/// assert_eq!(c.sum(0, 2, 3), 5);
/// assert_eq!(c.sum(0, 0, 2), 0);
/// ```
#[derive(Debug, Clone)]
pub struct ContourSum<T> {
    tree: CentroidTree,
    /// Tree `c`: values of the component of `c` by distance from `c`.
    own: FlatFenwick<T>,
    /// Tree `c`: values of the component of `c` by distance from `parent[c]`.
    up: FlatFenwick<T>,
}

impl<T: Copy + Default + AddAssign + Sub<Output = T>> ContourSum<T> {
    /// All values start at zero.
    pub fn new(tree: CentroidTree) -> Self {
        let own = FlatFenwick::new((0..tree.len()).map(|c| tree.radius[c] + 1));
        let up = FlatFenwick::new((0..tree.len()).map(|c| if tree.parent[c] == CentroidTree::NONE { 0 } else { tree.radius[tree.parent[c]] + 1 }));
        ContourSum { tree, own, up }
    }

    /// The underlying decomposition.
    pub fn tree(&self) -> &CentroidTree {
        &self.tree
    }

    /// Add `x` to the value of vertex `v`.
    pub fn add(&mut self, v: usize, x: T) {
        let mut level = self.tree.depth[v];
        let mut c = v;
        loop {
            self.own.add(c, self.tree.dist[level][v] as usize, x);
            if level == 0 {
                break;
            }
            self.up.add(c, self.tree.dist[level - 1][v] as usize, x);
            c = self.tree.parent[c];
            level -= 1;
        }
    }

    /// The sum of the values of all vertices `u` with `l <= dist(v, u) < r`.
    pub fn sum(&self, v: usize, l: usize, r: usize) -> T {
        let mut total = T::default();
        if l >= r {
            return total;
        }
        let mut level = self.tree.depth[v];
        let mut c = v;
        let mut below = CentroidTree::NONE;
        loop {
            let d = self.tree.dist[level][v] as usize;
            // Distances from `c` that put `u` at distance [l, r) from `v` through `c`.
            if r > d {
                let (lo, hi) = (l.saturating_sub(d), r - d);
                total += self.own.prefix(c, hi) - self.own.prefix(c, lo);
                if below != CentroidTree::NONE {
                    total = total - (self.up.prefix(below, hi) - self.up.prefix(below, lo));
                }
            }
            if level == 0 {
                break;
            }
            below = c;
            c = self.tree.parent[c];
            level -= 1;
        }
        total
    }
}

/// Add to all vertices within a distance range, read one vertex, on a fixed tree.
///
/// ```
/// use centroid_decomposition::{CentroidTree, ContourAdd};
/// let tree = CentroidTree::from_edges(4, &[(0, 1), (1, 2), (1, 3)]);
/// let mut c = ContourAdd::<i64>::new(tree);
/// c.add(0, 1, 3, 7); // every vertex at distance 1 or 2 from vertex 0: 1, 2, 3
/// assert_eq!((c.get(0), c.get(1), c.get(3)), (0, 7, 7));
/// ```
#[derive(Debug, Clone)]
pub struct ContourAdd<T> {
    tree: CentroidTree,
    /// Difference arrays by distance from `c` / from `parent[c]`.
    own: FlatFenwick<T>,
    up: FlatFenwick<T>,
}

impl<T: Copy + Default + AddAssign + Sub<Output = T>> ContourAdd<T> {
    /// All values start at zero.
    pub fn new(tree: CentroidTree) -> Self {
        let own = FlatFenwick::new((0..tree.len()).map(|c| tree.radius[c] + 2));
        let up = FlatFenwick::new((0..tree.len()).map(|c| if tree.parent[c] == CentroidTree::NONE { 0 } else { tree.radius[tree.parent[c]] + 2 }));
        ContourAdd { tree, own, up }
    }

    /// The underlying decomposition.
    pub fn tree(&self) -> &CentroidTree {
        &self.tree
    }

    /// Add `x` to every vertex `u` with `l <= dist(v, u) < r`.
    pub fn add(&mut self, v: usize, l: usize, r: usize, x: T) {
        if l >= r {
            return;
        }
        let neg = T::default() - x;
        let mut level = self.tree.depth[v];
        let mut c = v;
        let mut below = CentroidTree::NONE;
        loop {
            let d = self.tree.dist[level][v] as usize;
            if r > d {
                let (lo, hi) = (l.saturating_sub(d), r - d);
                self.own.add(c, lo, x);
                self.own.add(c, hi, neg);
                if below != CentroidTree::NONE {
                    // Cancel the contribution of `u` inside the component we came from.
                    self.up.add(below, lo, neg);
                    self.up.add(below, hi, x);
                }
            }
            if level == 0 {
                break;
            }
            below = c;
            c = self.tree.parent[c];
            level -= 1;
        }
    }

    /// The current value of vertex `v`.
    pub fn get(&self, v: usize) -> T {
        let mut total = T::default();
        let mut level = self.tree.depth[v];
        let mut c = v;
        let mut below = CentroidTree::NONE;
        loop {
            let d = self.tree.dist[level][v] as usize;
            total += self.own.prefix(c, d + 1);
            if below != CentroidTree::NONE {
                total += self.up.prefix(below, d + 1);
            }
            if level == 0 {
                break;
            }
            below = c;
            c = self.tree.parent[c];
            level -= 1;
        }
        total
    }
}
