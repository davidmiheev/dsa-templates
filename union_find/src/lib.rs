//! Generic Union-Find (Disjoint Set) data structure.
//!
//! Root and size share one array: following the trick used by KACTL's `UF`
//! and AtCoder Library's `Dsu`, a root stores `-size` and a non-root stores
//! its parent's index, both as `i64` in the same `Vec`. This halves the
//! per-node storage compared to keeping `root` and `size` as separate
//! arrays, while every node also carries a `character` value (an arbitrary
//! `usize` updated through the user-supplied `op` closure) that neither of
//! those references generalizes.
//!
//! `union(x, y, w)` merges the components of `x` and `y`, applying the closure
//! `op` to combine their characters along with the edge weight `w`. When `x` and
//! `y` are already in the same component, only `character[root] = op(character[root], w)`
//! is updated and the function returns `false`; otherwise it returns `true`.

/// Generic Union-Find with weighted union and per-component character.
pub struct UnionFind<F>
where
    F: Fn(usize, usize) -> usize,
{
    /// `parent_or_size[x] < 0` means `x` is a root and `-parent_or_size[x]`
    /// is the size of its component; otherwise it holds the parent of `x`.
    parent_or_size: Vec<i64>,
    components: usize,
    character: Vec<usize>,
    op: F,
}

impl<F> UnionFind<F>
where
    F: Fn(usize, usize) -> usize,
{
    /// Build a Union-Find over the indices `0..=n`. The closure `op` defines how
    /// character values combine during `union`. `init` seeds every component's
    /// character.
    pub fn new(n: usize, op: F, init: usize) -> Self {
        Self {
            parent_or_size: vec![-1; n + 1],
            components: n,
            character: vec![init; n + 1],
            op,
        }
    }

    /// Find the representative of `x` with path compression.
    pub fn find(&mut self, x: usize) -> usize {
        if self.parent_or_size[x] < 0 {
            return x;
        }
        let root = self.find(self.parent_or_size[x] as usize);
        self.parent_or_size[x] = root as i64;
        root
    }

    /// Union the components of `x` and `y`, combining characters via `op` and `w`.
    /// Returns `true` if a merge happened, `false` if `x` and `y` were already
    /// in the same component.
    pub fn union(&mut self, x: usize, y: usize, w: usize) -> bool {
        let root_x = self.find(x);
        let root_y = self.find(y);

        if root_x == root_y {
            self.character[root_x] = (self.op)(self.character[root_x], w);
            return false;
        }

        // Character combine order is always (x, y) regardless of which
        // root absorbs the other, since `op` need not be commutative.
        let combined =
            (self.op)((self.op)(self.character[root_x], self.character[root_y]), w);

        // Union by size: the larger component's root absorbs the smaller one.
        let (new_root, absorbed) =
            if self.parent_or_size[root_x] <= self.parent_or_size[root_y] {
                (root_x, root_y)
            } else {
                (root_y, root_x)
            };

        self.parent_or_size[new_root] += self.parent_or_size[absorbed];
        self.parent_or_size[absorbed] = new_root as i64;
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
    pub fn get_character(&mut self, x: usize) -> usize {
        let root_x = self.find(x);
        self.character[root_x]
    }

    /// Size of the component containing `x`.
    pub fn get_size(&mut self, x: usize) -> usize {
        let root_x = self.find(x);
        (-self.parent_or_size[root_x]) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn singleton_has_size_one() {
        let mut uf = UnionFind::new(5, |a, b| a + b, 0);
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
