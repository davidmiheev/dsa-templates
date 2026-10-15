//! A 2D Fenwick tree for point additions and rectangle sums when the points that will ever be
//! touched are known up front.
//!
//! The outer Fenwick tree runs over the distinct `x` coordinates; every outer node keeps a
//! sorted list of the `y` coordinates of the points in its range together with an inner Fenwick
//! tree over that list. Everything lives in a few flat arrays: `O(n log n)` memory for `n`
//! declared points, `O(log^2 n)` per update and query, no hashing.
//!
//! This is the *online* tool (queries may depend on earlier answers). When every update and
//! query is known in advance, [`dominance_sums`](crate::dominance_sums) is several times faster.
//!
//! # Example
//!
//! ```
//! use fenwick::Fenwick2D;
//!
//! // Declare every cell that will be updated.
//! let mut f = Fenwick2D::new(&[(1, 5), (3, 2), (3, 9)]);
//! f.add(1, 5, 10);
//! f.add(3, 2, 1);
//! f.add(3, 9, 100);
//! assert_eq!(f.prefix(4, 6), 11); // x < 4, y < 6
//! assert_eq!(f.rect(2, 4, 0, 10), 101); // 2 <= x < 4, 0 <= y < 10
//! assert_eq!(f.rect(0, 100, 0, 100), 111);
//! ```

use std::ops::{AddAssign, Sub};

/// 2D Fenwick tree over `T` values at declared `(x, y)` cells.
pub struct Fenwick2D<X, Y, T> {
    /// Sorted distinct x coordinates.
    xs: Vec<X>,
    /// Sorted distinct y coordinates; the node lists below store ranks into this array.
    all_ys: Vec<Y>,
    /// Outer node `j` (1-based) owns `ranks[start[j]..start[j + 1]]` (sorted) and the same range
    /// of `tree`, where it keeps an inner Fenwick tree (1-based inside the range).
    start: Vec<usize>,
    ranks: Vec<u32>,
    tree: Vec<T>,
}

impl<X: Ord + Copy, Y: Ord + Copy, T: Copy + Default + AddAssign + Sub<Output = T>> Fenwick2D<X, Y, T> {
    /// A tree of zeros that can be updated at the given cells (duplicates are fine).
    ///
    /// # Panics
    /// Panics if there are `2^32` or more points.
    pub fn new(points: &[(X, Y)]) -> Self {
        assert!(points.len() < u32::MAX as usize, "too many points");
        let mut xs: Vec<X> = points.iter().map(|p| p.0).collect();
        xs.sort_unstable();
        xs.dedup();
        let mut all_ys: Vec<Y> = points.iter().map(|p| p.1).collect();
        all_ys.sort_unstable();
        all_ys.dedup();
        let m = xs.len();
        let mut cells: Vec<(u32, usize)> = points
            .iter()
            .map(|&(x, y)| (all_ys.partition_point(|v| *v < y) as u32, xs.partition_point(|v| *v < x) + 1))
            .collect();
        cells.sort_unstable();
        cells.dedup();
        let mut start = vec![0usize; m + 2];
        for &(_, xi) in &cells {
            let mut j = xi;
            while j <= m {
                start[j + 1] += 1;
                j += j & j.wrapping_neg();
            }
        }
        for j in 1..=m + 1 {
            start[j] += start[j - 1];
        }
        let total = start[m + 1];
        // Filling in increasing rank keeps every node's list sorted.
        let mut fill = start.clone();
        let mut ranks = vec![0u32; total];
        for &(r, xi) in &cells {
            let mut j = xi;
            while j <= m {
                ranks[fill[j]] = r;
                fill[j] += 1;
                j += j & j.wrapping_neg();
            }
        }
        Fenwick2D { xs, all_ys, start, ranks, tree: vec![T::default(); total] }
    }

    /// Add `delta` at the cell `(x, y)`.
    ///
    /// # Panics
    /// Panics if the cell was not declared in [`new`](Self::new).
    pub fn add(&mut self, x: X, y: Y, delta: T) {
        let m = self.xs.len();
        let xi = self.xs.binary_search(&x).expect("cell was not declared in `new`") + 1;
        let r = self.all_ys.binary_search(&y).expect("cell was not declared in `new`") as u32;
        let mut j = xi;
        while j <= m {
            let (lo, hi) = (self.start[j], self.start[j + 1]);
            let k = self.ranks[lo..hi].partition_point(|v| *v < r);
            assert!(k < hi - lo && self.ranks[lo + k] == r, "cell was not declared in `new`");
            let mut i = k + 1;
            while i <= hi - lo {
                self.tree[lo + i - 1] += delta;
                i += i & i.wrapping_neg();
            }
            j += j & j.wrapping_neg();
        }
    }

    /// Sum of the inner prefix `[0, count)` of node `j`'s list restricted to ranks `< r`.
    fn node_prefix(&self, j: usize, r: u32) -> T {
        let (lo, hi) = (self.start[j], self.start[j + 1]);
        let mut i = self.ranks[lo..hi].partition_point(|v| *v < r);
        let mut sum = T::default();
        while i > 0 {
            sum += self.tree[lo + i - 1];
            i &= i - 1;
        }
        sum
    }

    /// The sum over all cells with `cx < x` and `cy < y`.
    pub fn prefix(&self, x: X, y: Y) -> T {
        let mut j = self.xs.partition_point(|v| *v < x);
        let r = self.all_ys.partition_point(|v| *v < y) as u32;
        let mut sum = T::default();
        while j > 0 {
            sum += self.node_prefix(j, r);
            j &= j - 1;
        }
        sum
    }

    /// The sum over all cells with `x_lo <= cx < x_hi` and `y_lo <= cy < y_hi` (zero if either
    /// range is empty).
    pub fn rect(&self, x_lo: X, x_hi: X, y_lo: Y, y_hi: Y) -> T {
        if x_lo >= x_hi || y_lo >= y_hi {
            return T::default();
        }
        let (mut hi, mut lo) = (self.xs.partition_point(|v| *v < x_hi), self.xs.partition_point(|v| *v < x_lo));
        let (r_lo, r_hi) = (self.all_ys.partition_point(|v| *v < y_lo) as u32, self.all_ys.partition_point(|v| *v < y_hi) as u32);
        let mut sum = T::default();
        // The Fenwick paths of `hi` and `lo` end in a common tail that cancels out.
        while hi > lo {
            sum += self.node_prefix(hi, r_hi) - self.node_prefix(hi, r_lo);
            hi &= hi - 1;
        }
        while lo > hi {
            let part = self.node_prefix(lo, r_hi) - self.node_prefix(lo, r_lo);
            sum = sum - part;
            lo &= lo - 1;
        }
        sum
    }
}
