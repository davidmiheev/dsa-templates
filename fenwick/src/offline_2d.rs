//! A 2D Fenwick tree for point additions and rectangle sums when the points that will ever be
//! touched are known up front.
//!
//! The outer Fenwick tree runs over the distinct `x` coordinates; every outer node keeps a
//! sorted list of the `y` coordinates of the points in its range together with an inner Fenwick
//! tree over that list. Everything lives in a few flat arrays: `O(n log n)` memory for `n`
//! declared points, `O(log^2 n)` per update and query, no hashing.
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
    /// Outer node `j` (1-based) owns `ys[start[j]..start[j + 1]]` (sorted) and the same range of
    /// `tree`, where it keeps an inner Fenwick tree (1-based inside the range).
    start: Vec<usize>,
    ys: Vec<Y>,
    tree: Vec<T>,
}

impl<X: Ord + Copy, Y: Ord + Copy, T: Copy + Default + AddAssign + Sub<Output = T>> Fenwick2D<X, Y, T> {
    /// A tree of zeros that can be updated at the given cells (duplicates are fine).
    pub fn new(points: &[(X, Y)]) -> Self {
        let mut xs: Vec<X> = points.iter().map(|p| p.0).collect();
        xs.sort_unstable();
        xs.dedup();
        let m = xs.len();
        let mut by_y: Vec<(Y, usize)> = points.iter().map(|&(x, y)| (y, xs.partition_point(|v| *v < x) + 1)).collect();
        by_y.sort_unstable_by_key(|p| p.0);
        by_y.dedup_by(|a, b| a.0 == b.0 && a.1 == b.1);
        let mut start = vec![0usize; m + 2];
        for &(_, xi) in &by_y {
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
        // Filling in increasing y keeps every node's list sorted.
        let mut fill = start.clone();
        let mut ys = Vec::with_capacity(total);
        if let Some(&(y0, _)) = by_y.first() {
            ys.resize(total, y0);
        }
        for &(y, xi) in &by_y {
            let mut j = xi;
            while j <= m {
                ys[fill[j]] = y;
                fill[j] += 1;
                j += j & j.wrapping_neg();
            }
        }
        Fenwick2D { xs, start, ys, tree: vec![T::default(); total] }
    }

    /// Add `delta` at the cell `(x, y)`.
    ///
    /// # Panics
    /// Panics if the cell was not declared in [`new`](Self::new).
    pub fn add(&mut self, x: X, y: Y, delta: T) {
        let m = self.xs.len();
        let xi = self.xs.binary_search(&x).expect("cell was not declared in `new`") + 1;
        let mut j = xi;
        while j <= m {
            let (lo, hi) = (self.start[j], self.start[j + 1]);
            let k = self.ys[lo..hi].partition_point(|v| *v < y);
            assert!(k < hi - lo && self.ys[lo + k] == y, "cell was not declared in `new`");
            let mut i = k + 1;
            while i <= hi - lo {
                self.tree[lo + i - 1] += delta;
                i += i & i.wrapping_neg();
            }
            j += j & j.wrapping_neg();
        }
    }

    /// The sum over all cells with `cx < x` and `cy < y`.
    pub fn prefix(&self, x: X, y: Y) -> T {
        let mut j = self.xs.partition_point(|v| *v < x);
        let mut sum = T::default();
        while j > 0 {
            let (lo, hi) = (self.start[j], self.start[j + 1]);
            let mut i = self.ys[lo..hi].partition_point(|v| *v < y);
            while i > 0 {
                sum += self.tree[lo + i - 1];
                i &= i - 1;
            }
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
        (self.prefix(x_hi, y_hi) - self.prefix(x_lo, y_hi)) - (self.prefix(x_hi, y_lo) - self.prefix(x_lo, y_lo))
    }
}
