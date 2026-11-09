//! A generic, iterative (bottom-up) segment tree.
//!
//! `T` is the element type and `op` an associative combining operation with identity
//! `identity` (sum with `0`, min with `T::MAX`, max with `T::MIN`, function composition with
//! the identity function, ...). `op` does **not** need to be commutative: queries combine
//! elements strictly left to right.
//!
//! The tree is stored in a flat `Vec<T>` of size `2 * n` (leaves at `n..2n`), so it uses half
//! the memory of the usual recursive `4 * n` layout and needs no recursion. `update` and
//! `query` run in `O(log n)`; building takes `O(n)`.
//!
//! Ranges are half-open `[l, r)`, as everywhere else in this repository's range structures.
//!
//! # Example
//!
//! ```
//! use segment_tree::SegmentTree;
//!
//! let mut st = SegmentTree::new(&[1, 3, 5, 7, 9, 11], 0, |a, b| a + b);
//! assert_eq!(st.query(1, 4), 15); // 3 + 5 + 7
//! st.update(2, 6);
//! assert_eq!(st.query(1, 4), 16);
//! assert_eq!(st.query(2, 2), 0); // empty range gives the identity
//! ```

/// Generic segment tree over an associative operation `F` with identity element.
pub struct SegmentTree<T, F> {
    n: usize,
    tree: Vec<T>,
    identity: T,
    op: F,
}

impl<T, F> SegmentTree<T, F>
where
    T: Copy,
    F: Fn(T, T) -> T,
{
    /// Build a tree over `arr` in `O(n)`. `identity` must satisfy `op(identity, x) == x ==
    /// op(x, identity)`, and `op` must be associative. An empty `arr` is allowed.
    pub fn new(arr: &[T], identity: T, op: F) -> Self {
        let n = arr.len();
        let mut tree = vec![identity; 2 * n];
        tree[n..].copy_from_slice(arr);
        for i in (1..n).rev() {
            tree[i] = op(tree[2 * i], tree[2 * i + 1]);
        }
        SegmentTree { n, tree, identity, op }
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// The element at index `i`.
    ///
    /// # Panics
    /// Panics if `i >= len()`.
    pub fn get(&self, i: usize) -> T {
        assert!(i < self.n, "index {i} out of range for length {}", self.n);
        self.tree[self.n + i]
    }

    /// Point update: set element `i` to `val` and recompute its ancestors.
    ///
    /// # Panics
    /// Panics if `i >= len()`.
    pub fn update(&mut self, i: usize, val: T) {
        assert!(i < self.n, "index {i} out of range for length {}", self.n);
        let mut p = i + self.n;
        self.tree[p] = val;
        while p > 1 {
            p >>= 1;
            self.tree[p] = (self.op)(self.tree[2 * p], self.tree[2 * p + 1]);
        }
    }

    /// `op(a[l], a[l + 1], ..., a[r - 1])` over the half-open range `[l, r)`;
    /// the identity for an empty range.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn query(&self, l: usize, r: usize) -> T {
        assert!(l <= r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
        let (mut l, mut r) = (l + self.n, r + self.n);
        let (mut left, mut right) = (self.identity, self.identity);
        while l < r {
            if l & 1 == 1 {
                left = (self.op)(left, self.tree[l]);
                l += 1;
            }
            if r & 1 == 1 {
                r -= 1;
                right = (self.op)(self.tree[r], right);
            }
            l >>= 1;
            r >>= 1;
        }
        (self.op)(left, right)
    }

    /// `op` over the whole array (the identity if empty).
    pub fn all_prod(&self) -> T {
        self.query(0, self.n)
    }
}
