//! "Parallel" union of ranges: merge `a + i` with `b + i` for every `i < k` in `O(log n)`
//! amortized, using one union-find per power of two.
//!
//! Level `l` has one node per start position `i`, standing for the window `[i, i + 2^l)`.
//! Uniting windows `(a, b)` at level `l` and the two halves one level below are equivalent,
//! so a range of length `k` is two overlapping windows of `2^floor(log2 k)` (like a sparse
//! table), and a window union that finds its two nodes already connected stops immediately.
//! Each level is a plain DSU whose merges happen at most `n` times, so the total work for any
//! sequence of range unions is `O(n log n α(n))`.
//!
//! Only level 0 is the real vertex-level union-find; [`RangeParallelUnionFind::unite_range`]
//! reports each level-0 merge to the caller so it can maintain per-component aggregates.
//!
//! # Example
//!
//! ```
//! use union_find::RangeParallelUnionFind;
//!
//! let mut uf = RangeParallelUnionFind::new(6);
//! let mut merges = vec![];
//! // 0~3, 1~4, 2~5
//! uf.unite_range(0, 3, 3, |absorbed, root| merges.push((absorbed, root)));
//! assert_eq!(merges.len(), 3);
//! assert!(uf.same(1, 4) && !uf.same(0, 1));
//! ```

/// A stack of union-finds supporting `unite_range`.
pub struct RangeParallelUnionFind {
    n: usize,
    /// `levels[l][i] < 0`: root with size `-value` for the window `[i, i + 2^l)`.
    levels: Vec<Vec<i32>>,
}

fn find(p: &mut [i32], x: usize) -> usize {
    let mut root = x;
    while p[root] >= 0 {
        root = p[root] as usize;
    }
    let mut cur = x;
    while cur != root {
        let next = p[cur] as usize;
        p[cur] = root as i32;
        cur = next;
    }
    root
}

impl RangeParallelUnionFind {
    /// `n` singleton vertices.
    pub fn new(n: usize) -> Self {
        assert!(n < i32::MAX as usize, "too many elements: {n}");
        let top = if n == 0 { 0 } else { (usize::BITS - 1 - n.leading_zeros()) as usize };
        Self { n, levels: (0..=top).map(|_| vec![-1; n]).collect() }
    }

    /// Representative of vertex `x`.
    pub fn find(&mut self, x: usize) -> usize {
        find(&mut self.levels[0], x)
    }

    pub fn same(&mut self, x: usize, y: usize) -> bool {
        self.find(x) == self.find(y)
    }

    /// Size of the component of `x`.
    pub fn size(&mut self, x: usize) -> usize {
        let r = self.find(x);
        (-self.levels[0][r]) as usize
    }

    /// Add the edges `(a + i, b + i)` for `i = 0..k`.
    ///
    /// `on_merge(absorbed, root)` is called for each *vertex-level* merge, just after it
    /// happens, with the old root that was absorbed and the surviving root. The caller can fold
    /// its per-component data of `absorbed` into `root`'s.
    ///
    /// # Panics
    /// Panics unless `a + k <= n` and `b + k <= n`.
    pub fn unite_range<G: FnMut(usize, usize)>(&mut self, a: usize, b: usize, k: usize, mut on_merge: G) {
        assert!(a + k <= self.n && b + k <= self.n, "range out of bounds: a={a} b={b} k={k} n={}", self.n);
        if k == 0 || a == b {
            return;
        }
        let level = (usize::BITS - 1 - k.leading_zeros()) as usize;
        let width = 1 << level;
        self.unite_window(level, a, b, &mut on_merge);
        self.unite_window(level, a + k - width, b + k - width, &mut on_merge);
    }

    fn unite_window<G: FnMut(usize, usize)>(&mut self, level: usize, a: usize, b: usize, on_merge: &mut G) {
        // Iterative: each popped (level, a, b) merges at most one pair of nodes and pushes the
        // two child windows, so the stack holds at most 2 * level + 1 entries.
        let mut stack = vec![(level, a, b)];
        while let Some((l, a, b)) = stack.pop() {
            let p = &mut self.levels[l];
            let (mut ra, mut rb) = (find(p, a), find(p, b));
            if ra == rb {
                continue;
            }
            if p[ra] > p[rb] {
                std::mem::swap(&mut ra, &mut rb); // ra is the larger root
            }
            p[ra] += p[rb];
            p[rb] = ra as i32;
            if l == 0 {
                on_merge(rb, ra);
            } else {
                let half = 1 << (l - 1);
                stack.push((l - 1, a, b));
                stack.push((l - 1, a + half, b + half));
            }
        }
    }
}
