//! Fenwick Tree (Binary Indexed Tree)
//!
//! Fenwick Tree is a data structure that can efficiently update elements
//! and calculate prefix sums in a table of numbers.
//!
//! Given an array of `n` numbers, the prefix sum up to index `i` is the
//! sum of the first `i` elements. Fenwick Tree supports both update and
//! prefix-sum queries in `O(log n)` time by storing partial sums in a tree
//! whose node ranges are determined by the least significant bit of each
//! index.

use std::ops::{AddAssign, Sub};

/// Fenwick Tree (Binary Indexed Tree).
///
/// `T` is the numeric type stored in the tree. It must be `Copy`,
/// `Default` (the additive identity, typically `0`), and support
/// `+=` and `-`.
#[derive(Debug, Clone)]
pub struct FenwickTree<T> {
    size: usize,
    tree: Vec<T>,
}

impl<T> FenwickTree<T>
where
    T: Copy + Default + AddAssign + Sub<Output = T>,
{
    /// Build a Fenwick Tree from an initial slice of values.
    ///
    /// Runs in `O(n log n)`. For an `O(n)` build, build an empty tree
    /// and call [`Self::update`] for each element.
    pub fn new(nums: &[T]) -> Self {
        let size = nums.len();
        let tree = vec![T::default(); size + 1];

        let mut ft = Self { size, tree };

        for (i, &val) in nums.iter().enumerate() {
            ft.update(i, val);
        }
        ft
    }

    /// Add `delta` to the element at `index`.
    ///
    /// Panics if `index >= self.size`.
    pub fn update(&mut self, index: usize, delta: T) {
        if index >= self.size {
            panic!("Index out of bounds");
        }

        let mut i = index + 1;
        loop {
            if i > self.size { break }
            self.tree[i] += delta;
            i += i & (!i + 1); // equivalent to i & -i
        }
    }

    /// Compute the prefix sum of elements `[0, index]` (inclusive).
    ///
    /// Panics if `index >= self.size`.
    pub fn pref(&self, index: usize) -> T {
        if index >= self.size {
            panic!("Index out of bounds");
        }
        let mut result = T::default();
        let mut i = index + 1;
        loop {
            if i <= 0 { break }
            result += self.tree[i];
            i -= i & (!i + 1); // equivalent to i & -i
        }
        result
    }

    /// Compute the sum of elements in the half-open range `[left, right]`.
    ///
    /// Panics if `right < left` or `right >= self.size`.
    pub fn sum_range(&self, left: usize, right: usize) -> T {
        if right < left || right >= self.size {
            panic!("Invalid range");
        }
        match left == 0 {
            true => self.pref(right),
            false => self.pref(right) - self.pref(left - 1)
        }
    }
}