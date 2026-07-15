//! Generic Union-Find (Disjoint Set) data structure.
//!
//! Each component carries:
//! - a `size` (used for union-by-size),
//! - a `character` value (an arbitrary `usize` updated through the user-supplied `op` closure).
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
    root: Vec<usize>,
    size: Vec<usize>,
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
            root: (0..=n).collect(),
            size: vec![1; n + 1],
            components: n,
            character: vec![init; n + 1],
            op,
        }
    }

    /// Find the representative of `x` with path compression.
    pub fn find(&mut self, x: usize) -> usize {
        if self.root[x] == x {
            return x;
        }
        self.root[x] = self.find(self.root[x]);
        self.root[x]
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

        if self.size[root_x] > self.size[root_y] {
            self.size[root_x] += self.size[root_y];
            self.character[root_x] =
                (self.op)((self.op)(self.character[root_x], self.character[root_y]), w);
            self.root[root_y] = root_x;
        } else {
            self.size[root_y] += self.size[root_x];
            self.character[root_y] =
                (self.op)((self.op)(self.character[root_x], self.character[root_y]), w);
            self.root[root_x] = root_y;
        }
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
        self.size[root_x]
    }
}