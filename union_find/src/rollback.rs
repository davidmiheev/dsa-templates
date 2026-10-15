//! Union-Find with undo (no path compression, so every union can be reverted).
//!
//! `find` takes `O(log n)` (union by size), `union` pushes one record on a history stack, and
//! [`RollbackUnionFind::rollback`] pops back to a [`RollbackUnionFind::snapshot`]. This is the
//! tool for offline "versioned" or divide-and-conquer-over-time problems: walk a tree of
//! versions depth first, `union` on the way down and `rollback` on the way up.
//!
//! # Example
//!
//! ```
//! use union_find::RollbackUnionFind;
//!
//! let mut uf = RollbackUnionFind::new(4);
//! let before = uf.snapshot();
//! assert!(uf.union(0, 1));
//! assert!(uf.same(0, 1));
//! uf.rollback(before);
//! assert!(!uf.same(0, 1));
//! ```

/// Union-find supporting `rollback`.
#[derive(Clone, Debug)]
pub struct RollbackUnionFind {
    /// `< 0`: root with size `-value`; otherwise the parent.
    parent_or_size: Vec<i32>,
    /// `(absorbed root, new root, absorbed root's old parent_or_size)` per successful union.
    history: Vec<(u32, u32, i32)>,
    components: usize,
}

impl RollbackUnionFind {
    /// `n` singleton vertices.
    pub fn new(n: usize) -> Self {
        assert!(n < i32::MAX as usize, "too many elements: {n}");
        Self { parent_or_size: vec![-1; n], history: Vec::new(), components: n }
    }

    /// Representative of `x` (`O(log n)`).
    pub fn find(&self, mut x: usize) -> usize {
        while self.parent_or_size[x] >= 0 {
            x = self.parent_or_size[x] as usize;
        }
        x
    }

    pub fn same(&self, x: usize, y: usize) -> bool {
        self.find(x) == self.find(y)
    }

    /// Size of the component of `x`.
    pub fn size(&self, x: usize) -> usize {
        (-self.parent_or_size[self.find(x)]) as usize
    }

    /// Number of components.
    pub fn components(&self) -> usize {
        self.components
    }

    /// Merge the components of `x` and `y`; returns whether they were different. A union that
    /// changes nothing leaves no history entry.
    pub fn union(&mut self, x: usize, y: usize) -> bool {
        let (mut a, mut b) = (self.find(x), self.find(y));
        if a == b {
            return false;
        }
        if self.parent_or_size[a] > self.parent_or_size[b] {
            std::mem::swap(&mut a, &mut b); // a is the larger root
        }
        self.history.push((b as u32, a as u32, self.parent_or_size[b]));
        self.parent_or_size[a] += self.parent_or_size[b];
        self.parent_or_size[b] = a as i32;
        self.components -= 1;
        true
    }

    /// A marker for the current state, to pass to [`Self::rollback`].
    pub fn snapshot(&self) -> usize {
        self.history.len()
    }

    /// Undo the most recent successful union.
    ///
    /// # Panics
    /// Panics if there is nothing to undo.
    pub fn undo(&mut self) {
        let (b, a, old_b) = self.history.pop().expect("nothing to undo");
        let (a, b) = (a as usize, b as usize);
        self.parent_or_size[a] -= old_b; // `a` stops counting the absorbed component
        self.parent_or_size[b] = old_b; // `b` is a root of its own component again
        self.components += 1;
    }

    /// Undo successful unions until the state of `snapshot` is restored.
    ///
    /// # Panics
    /// Panics if `snapshot` is newer than the current state.
    pub fn rollback(&mut self, snapshot: usize) {
        assert!(snapshot <= self.history.len(), "snapshot {snapshot} is from the future");
        while self.history.len() > snapshot {
            self.undo();
        }
    }
}
