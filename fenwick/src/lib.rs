//! Fenwick Tree (Binary Indexed Tree)
//!
//! Fenwick Tree is a data structure that can efficiently update elements
//! and calculate prefix sums in a table of numbers.
//!
//! Given an array of `n` numbers, the prefix sum up to index `i` is the
//! sum of the first `i` elements. Fenwick Tree supports both update and
//! prefix-sum queries in `O(log n)` time by storing partial sums in a tree
//! whose node ranges are determined by the least significant bit of each
//! index. Building from a slice takes `O(n)`.
//!
//! [`FenwickTree::lower_bound`] walks down the tree to find, in `O(log n)`, the first prefix
//! whose sum reaches a target: "k-th smallest" queries over a 0/1 presence array.
//!
//! # Example
//!
//! ```
//! use fenwick::FenwickTree;
//!
//! let mut ft = FenwickTree::new(&[1i64, 3, 5, 7, 9]);
//! assert_eq!(ft.pref(2), 9); // 1 + 3 + 5
//! ft.update(1, 10);
//! assert_eq!(ft.sum_range(1, 3), 25); // 13 + 5 + 7
//! assert_eq!(ft.lower_bound(15), 2); // prefix sums 1, 14, 19, ...: the first >= 15 is at index 2
//! ```

use std::ops::{AddAssign, Sub};

/// Fenwick Tree (Binary Indexed Tree).
///
/// `T` is the numeric type stored in the tree. It must be `Copy`,
/// `Default` (the additive identity, typically `0`), and support
/// `+=` and `-`.
#[derive(Debug, Clone)]
pub struct FenwickTree<T> {
    size: usize,
    /// 1-based: `tree[i]` holds the sum of the `lowbit(i)` elements ending at `i`.
    tree: Vec<T>,
}

impl<T> FenwickTree<T>
where
    T: Copy + Default + AddAssign + Sub<Output = T>,
{
    /// Build a Fenwick Tree from an initial slice of values in `O(n)`.
    pub fn new(nums: &[T]) -> Self {
        let size = nums.len();
        let mut tree = vec![T::default(); size + 1];
        tree[1..].copy_from_slice(nums);
        for i in 1..=size {
            let parent = i + (i & i.wrapping_neg());
            if parent <= size {
                let v = tree[i];
                tree[parent] += v;
            }
        }
        Self { size, tree }
    }

    /// A tree of `n` zeros (`T::default()`).
    pub fn with_len(n: usize) -> Self {
        Self { size: n, tree: vec![T::default(); n + 1] }
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.size
    }

    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    /// Add `delta` to the element at `index`.
    ///
    /// # Panics
    /// Panics if `index >= len()`.
    pub fn update(&mut self, index: usize, delta: T) {
        assert!(index < self.size, "Index out of bounds");
        let mut i = index + 1;
        while i <= self.size {
            self.tree[i] += delta;
            i += i & i.wrapping_neg();
        }
    }

    /// Sum of the elements `0..=index`.
    ///
    /// # Panics
    /// Panics if `index >= len()`.
    pub fn pref(&self, index: usize) -> T {
        assert!(index < self.size, "Index out of bounds");
        self.sum_first(index + 1)
    }

    /// Sum of the first `count` elements (`0 <= count <= len()`); zero for `count == 0`.
    ///
    /// # Panics
    /// Panics if `count > len()`.
    pub fn sum_first(&self, count: usize) -> T {
        assert!(count <= self.size, "Index out of bounds");
        let mut result = T::default();
        let mut i = count;
        while i > 0 {
            result += self.tree[i];
            i &= i - 1;
        }
        result
    }

    /// Sum of the elements `left..=right` (inclusive on both ends).
    ///
    /// # Panics
    /// Panics unless `left <= right < len()`.
    pub fn sum_range(&self, left: usize, right: usize) -> T {
        assert!(left <= right && right < self.size, "Invalid range");
        self.sum_first(right + 1) - self.sum_first(left)
    }

    /// Smallest index `i` such that `pref(i) >= target`, or `len()` if no prefix reaches
    /// `target`. For `target <= 0` this is `0`. All elements must be non-negative, so that
    /// prefix sums are non-decreasing.
    ///
    /// With a 0/1 presence array this is the "`target`-th present element" query.
    pub fn lower_bound(&self, target: T) -> usize
    where
        T: PartialOrd,
    {
        if target <= T::default() {
            return 0;
        }
        let (mut pos, mut rest) = (0usize, target);
        let mut step = if self.size == 0 { 0 } else { 1usize << (usize::BITS - 1 - self.size.leading_zeros()) };
        while step > 0 {
            let next = pos + step;
            if next <= self.size && self.tree[next] < rest {
                pos = next;
                rest = rest - self.tree[next];
            }
            step >>= 1;
        }
        pos
    }
}
