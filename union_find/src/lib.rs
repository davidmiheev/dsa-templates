//! Generic Union-Find (Disjoint Set) data structure.
//!
//! Root and size share one array: following the trick used by KACTL's `UF` and AtCoder
//! Library's `Dsu`, a root stores `-size` and a non-root stores its parent's index, both as
//! `i32` in the same `Vec`. That is 4 bytes per node, with iterative path compression.
//!
//! On top of connectivity and component sizes, every component carries a `character` value of
//! any `Copy` type `C`, combined by a user-supplied closure `op`.
//! `union(x, y, w)` merges the components of `x` and `y`, applying `op` to combine their
//! characters along with the edge weight `w`. When `x` and `y` are already in the same
//! component only `character[root] = op(character[root], w)` is updated and the function
//! returns `false`; otherwise it returns `true`.
//!
//! When no character is needed use [`plain`]: its character type is `()`, so the character
//! array has zero size and all character work compiles away, giving a classic DSU with
//! [`UnionFind::merge`].
//!
//! Variants for harder problems live in the same crate:
//! [`PotentialUnionFind`] (relations in an arbitrary group),
//! [`RollbackUnionFind`] (undo, for offline version trees), and
//! [`RangeParallelUnionFind`] (`a + i ~ b + i` for a whole range in `O(log n)` amortized).
//!
//! # Example
//!
//! ```
//! let mut uf = union_find::plain(5);
//! assert!(uf.merge(1, 2));
//! assert!(!uf.merge(2, 1));
//! assert!(uf.are_connected(1, 2) && !uf.are_connected(1, 3));
//! assert_eq!(uf.get_size(1), 2);
//! ```

mod potential;
mod range_parallel;
mod rollback;

pub use potential::PotentialUnionFind;
pub use range_parallel::RangeParallelUnionFind;
pub use rollback::{RollbackSumUnionFind, RollbackUnionFind};

/// Union-Find with per-component character values of type `C`.
pub struct UnionFind<C, F> {
    /// `parent_or_size[x] < 0` means `x` is a root and `-parent_or_size[x]` is the size of its
    /// component; otherwise it holds the parent of `x`.
    parent_or_size: Vec<i32>,
    components: usize,
    character: Vec<C>,
    op: F,
}

/// A Union-Find without characters (`C = ()`): a plain disjoint-set structure over `0..=n`.
pub fn plain(n: usize) -> UnionFind<(), impl Fn((), ())> {
    UnionFind::new(n, |_, _| (), ())
}

impl<C: Copy, F: Fn(C, C) -> C> UnionFind<C, F> {
    /// Build a Union-Find over the indices `0..=n`. The closure `op` defines how character
    /// values combine during `union`. `init` seeds every component's character.
    ///
    /// # Panics
    /// Panics if `n >= i32::MAX`.
    pub fn new(n: usize, op: F, init: C) -> Self {
        assert!(n < i32::MAX as usize, "too many elements: {n}");
        Self { parent_or_size: vec![-1; n + 1], components: n, character: vec![init; n + 1], op }
    }

    /// Find the representative of `x` with path compression.
    pub fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent_or_size[root] >= 0 {
            root = self.parent_or_size[root] as usize;
        }
        let mut cur = x;
        while cur != root {
            let next = self.parent_or_size[cur] as usize;
            self.parent_or_size[cur] = root as i32;
            cur = next;
        }
        root
    }

    /// Union the components of `x` and `y`, combining characters via `op` and `w`.
    /// Returns `true` if a merge happened, `false` if `x` and `y` were already
    /// in the same component.
    pub fn union(&mut self, x: usize, y: usize, w: C) -> bool {
        let root_x = self.find(x);
        let root_y = self.find(y);

        if root_x == root_y {
            self.character[root_x] = (self.op)(self.character[root_x], w);
            return false;
        }

        // Character combine order is always (x, y) regardless of which
        // root absorbs the other, since `op` need not be commutative.
        let combined = (self.op)((self.op)(self.character[root_x], self.character[root_y]), w);

        // Union by size: the larger component's root absorbs the smaller one.
        let (new_root, absorbed) =
            if self.parent_or_size[root_x] <= self.parent_or_size[root_y] { (root_x, root_y) } else { (root_y, root_x) };

        self.parent_or_size[new_root] += self.parent_or_size[absorbed];
        self.parent_or_size[absorbed] = new_root as i32;
        self.character[new_root] = combined;

        self.components -= 1;
        true
    }

    /// Whether `x` and `y` are in the same component.
    pub fn are_connected(&mut self, x: usize, y: usize) -> bool {
        self.find(x) == self.find(y)
    }

    /// Number of distinct components remaining.
    pub fn get_components(&self) -> usize {
        self.components
    }

    /// Representative of `x`.
    pub fn get_root(&mut self, x: usize) -> usize {
        self.find(x)
    }

    /// Character value associated with the component containing `x`.
    pub fn get_character(&mut self, x: usize) -> C {
        let root_x = self.find(x);
        self.character[root_x]
    }

    /// Size of the component containing `x`.
    pub fn get_size(&mut self, x: usize) -> usize {
        let root_x = self.find(x);
        (-self.parent_or_size[root_x]) as usize
    }
}

impl<F: Fn((), ())> UnionFind<(), F> {
    /// Merge the components of `x` and `y` (no character involved). Returns `true` if they
    /// were different components.
    pub fn merge(&mut self, x: usize, y: usize) -> bool {
        self.union(x, y, ())
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn singleton_has_size_one() {
        let mut uf = UnionFind::new(5, |a: usize, b: usize| a + b, 0);
        assert_eq!(uf.get_size(3), 1);
        assert_eq!(uf.get_root(3), 3);
    }

    #[test]
    fn union_by_size_tracks_component_size() {
        let mut uf = UnionFind::new(6, |a, b| a + b, 0);
        uf.union(0, 1, 0);
        uf.union(1, 2, 0);
        assert_eq!(uf.get_size(0), 3);
        assert_eq!(uf.get_size(1), 3);
        assert_eq!(uf.get_size(2), 3);
        assert_eq!(uf.get_size(4), 1);

        uf.union(3, 4, 0);
        uf.union(0, 3, 0);
        assert_eq!(uf.get_size(4), 5);
        assert_eq!(uf.get_root(0), uf.get_root(4));
    }

    #[test]
    fn character_combine_order_is_not_swapped_by_which_root_absorbs() {
        // Subtraction is not commutative, so this pins down that `op` is
        // always applied as op(char[x], char[y]) regardless of which
        // component ends up as the new root.
        let mut uf = UnionFind::new(3, |a: usize, b: usize| a.wrapping_sub(b), 10);
        // Force component {1} to be the larger one so it becomes the root
        // when merged with {0}.
        uf.union(1, 2, 0);
        let before = uf.get_character(1);
        uf.union(0, 1, 0);
        assert_eq!(uf.get_character(0), (10usize).wrapping_sub(before));
    }
}
