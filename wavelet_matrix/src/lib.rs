//! Wavelet matrix over a static array: order statistics and frequencies on ranges.
//!
//! Values are first replaced by their rank among the distinct values (so any `Ord + Copy` type
//! works, negative numbers included), and the matrix has one level per bit of that rank, most
//! significant first. At each level the elements are stably partitioned by the current bit
//! (zeros first), and a rank dictionary over the bit column lets a range `[l, r)` be mapped to
//! the corresponding range of the next level in `O(1)`. All queries therefore take
//! `O(log(distinct values))`:
//!
//! * [`WaveletMatrix::kth_smallest`]: the `k`-th smallest value in `a[l..r]`,
//! * [`WaveletMatrix::count_less`] / [`count_less_eq`](WaveletMatrix::count_less_eq):
//!   how many values in `a[l..r]` are `< x` / `<= x`,
//! * [`WaveletMatrix::count_equal`] and [`WaveletMatrix::range_freq`] (values in `[lo, hi)`).
//!
//! [`WeightedWaveletMatrix`] additionally stores prefix sums of per-element weights on every
//! level, which gives "count and sum of the elements `<= x` in a range" in the same time.
//!
//! # Example
//!
//! ```
//! use wavelet_matrix::WaveletMatrix;
//!
//! let a = [5, 1, 4, 1, 5, 9, 2, 6];
//! let w = WaveletMatrix::new(&a);
//! assert_eq!(w.kth_smallest(0, 8, 0), 1); // minimum
//! assert_eq!(w.kth_smallest(2, 6, 1), 4); // sorted a[2..6] = [1, 4, 5, 9]
//! assert_eq!(w.count_less(0, 8, 5), 4);   // 1, 4, 1, 2
//! assert_eq!(w.count_equal(0, 8, 5), 2);
//! assert_eq!(w.range_freq(0, 8, 2, 6), 4); // values in [2, 6): 5, 4, 5, 2
//! ```

use std::ops::{Add, Sub};

/// Bit column with `O(1)` rank.
#[derive(Clone, Debug)]
struct BitColumn {
    words: Vec<u64>,
    /// `before[w]` = number of ones in `words[..w]`.
    before: Vec<u32>,
}

impl BitColumn {
    fn new(bits: &[bool]) -> Self {
        let mut words = vec![0u64; bits.len() / 64 + 1];
        for (i, &b) in bits.iter().enumerate() {
            if b {
                words[i >> 6] |= 1 << (i & 63);
            }
        }
        let mut before = Vec::with_capacity(words.len() + 1);
        let mut acc = 0u32;
        for w in &words {
            before.push(acc);
            acc += w.count_ones();
        }
        before.push(acc);
        BitColumn { words, before }
    }

    /// Number of ones among the first `i` bits.
    #[inline]
    fn rank1(&self, i: usize) -> usize {
        let (w, b) = (i >> 6, i & 63);
        let mask = if b == 0 { 0 } else { u64::MAX >> (64 - b) };
        self.before[w] as usize + (self.words[w] & mask).count_ones() as usize
    }
}

#[derive(Clone, Debug)]
struct Level {
    column: BitColumn,
    /// Number of zero bits (the offset at which the ones start on the next level).
    zeros: usize,
}

/// Static array with range order-statistic and frequency queries.
#[derive(Clone, Debug)]
pub struct WaveletMatrix<T> {
    n: usize,
    /// Distinct values, ascending: the value of rank `r` is `vals[r]`.
    vals: Vec<T>,
    /// Levels from the most significant bit down.
    levels: Vec<Level>,
}

impl<T: Ord + Copy> WaveletMatrix<T> {
    /// Build over `a` in `O(n log n)`.
    pub fn new(a: &[T]) -> Self {
        Self::build(a).0
    }

    /// Build the matrix and also return, for every level, the arrangement of original indices
    /// (needed by [`WeightedWaveletMatrix`]).
    fn build(a: &[T]) -> (Self, Vec<Vec<u32>>) {
        let n = a.len();
        let mut vals: Vec<T> = a.to_vec();
        vals.sort_unstable();
        vals.dedup();
        let depth = if vals.len() <= 1 { 0 } else { (usize::BITS - (vals.len() - 1).leading_zeros()) as usize };
        let rank_of = |x: &T| vals.partition_point(|v| v < x) as u32;
        let mut cur: Vec<(u32, u32)> = a.iter().enumerate().map(|(i, x)| (rank_of(x), i as u32)).collect();
        let mut levels = Vec::with_capacity(depth);
        let mut arrangements = Vec::with_capacity(depth + 1);
        arrangements.push(cur.iter().map(|e| e.1).collect::<Vec<u32>>());
        for level in (0..depth).rev() {
            let bits: Vec<bool> = cur.iter().map(|e| (e.0 >> level) & 1 == 1).collect();
            let zeros = bits.iter().filter(|&&b| !b).count();
            let (mut next, mut ones) = (Vec::with_capacity(n), Vec::new());
            for (&b, &e) in bits.iter().zip(&cur) {
                if b { ones.push(e) } else { next.push(e) }
            }
            next.extend(ones);
            levels.push(Level { column: BitColumn::new(&bits), zeros });
            cur = next;
            arrangements.push(cur.iter().map(|e| e.1).collect::<Vec<u32>>());
        }
        (WaveletMatrix { n, vals, levels }, arrangements)
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    fn check_range(&self, l: usize, r: usize) {
        assert!(l <= r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
    }

    /// The `k`-th smallest value (0-based) among `a[l..r]`.
    ///
    /// # Panics
    /// Panics unless `l < r <= len()` and `k < r - l`.
    pub fn kth_smallest(&self, mut l: usize, mut r: usize, mut k: usize) -> T {
        assert!(l < r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
        assert!(k < r - l, "k = {k} out of range for a range of length {}", r - l);
        let depth = self.levels.len();
        let mut rank = 0usize;
        for (i, level) in self.levels.iter().enumerate() {
            let (ones_l, ones_r) = (level.column.rank1(l), level.column.rank1(r));
            let (zeros_l, zeros_r) = (l - ones_l, r - ones_r);
            if k < zeros_r - zeros_l {
                l = zeros_l;
                r = zeros_r;
            } else {
                k -= zeros_r - zeros_l;
                rank |= 1 << (depth - 1 - i);
                l = level.zeros + ones_l;
                r = level.zeros + ones_r;
            }
        }
        self.vals[rank]
    }

    /// Number of elements of `a[l..r]` whose *rank* is below `rank` (`0 <= rank <= vals.len()`).
    fn count_below_rank(&self, mut l: usize, mut r: usize, rank: usize) -> usize {
        if rank == 0 {
            return 0;
        }
        if rank >= self.vals.len() {
            return r - l;
        }
        let depth = self.levels.len();
        let mut count = 0;
        for (i, level) in self.levels.iter().enumerate() {
            let (ones_l, ones_r) = (level.column.rank1(l), level.column.rank1(r));
            if (rank >> (depth - 1 - i)) & 1 == 1 {
                count += (r - ones_r) - (l - ones_l); // elements with a 0 here are smaller
                l = level.zeros + ones_l;
                r = level.zeros + ones_r;
            } else {
                l -= ones_l;
                r -= ones_r;
            }
        }
        count
    }

    /// Number of elements of `a[l..r]` that are `< x`.
    pub fn count_less(&self, l: usize, r: usize, x: T) -> usize {
        self.check_range(l, r);
        self.count_below_rank(l, r, self.vals.partition_point(|v| *v < x))
    }

    /// Number of elements of `a[l..r]` that are `<= x`.
    pub fn count_less_eq(&self, l: usize, r: usize, x: T) -> usize {
        self.check_range(l, r);
        self.count_below_rank(l, r, self.vals.partition_point(|v| *v <= x))
    }

    /// Number of elements of `a[l..r]` equal to `x`.
    pub fn count_equal(&self, mut l: usize, mut r: usize, x: T) -> usize {
        self.check_range(l, r);
        let rank = self.vals.partition_point(|v| *v < x);
        if rank >= self.vals.len() || self.vals[rank] != x {
            return 0;
        }
        // One walk down to the leaf of `rank`: what remains of the range is exactly its copies.
        let depth = self.levels.len();
        for (i, level) in self.levels.iter().enumerate() {
            let (ones_l, ones_r) = (level.column.rank1(l), level.column.rank1(r));
            if (rank >> (depth - 1 - i)) & 1 == 1 {
                l = level.zeros + ones_l;
                r = level.zeros + ones_r;
            } else {
                l -= ones_l;
                r -= ones_r;
            }
        }
        r - l
    }

    /// Number of elements of `a[l..r]` with `lo <= value < hi` (zero if `lo >= hi`).
    pub fn range_freq(&self, l: usize, r: usize, lo: T, hi: T) -> usize {
        self.check_range(l, r);
        if lo >= hi {
            return 0;
        }
        self.count_less(l, r, hi) - self.count_less(l, r, lo)
    }
}

/// A [`WaveletMatrix`] whose elements also carry a weight, so that sums can be queried too.
#[derive(Clone, Debug)]
pub struct WeightedWaveletMatrix<T, W> {
    core: WaveletMatrix<T>,
    /// `prefix[level][i]`: total weight of the first `i` elements in the arrangement *after*
    /// `level` partitions (level 0 = the original order, one more array than `levels`).
    prefix: Vec<Vec<W>>,
}

impl<T, W> WeightedWaveletMatrix<T, W>
where
    T: Ord + Copy,
    W: Copy + Default + Add<Output = W> + Sub<Output = W>,
{
    /// Build over `a` with `weight[i]` attached to `a[i]`.
    ///
    /// # Panics
    /// Panics if the two slices differ in length.
    pub fn new(a: &[T], weight: &[W]) -> Self {
        assert_eq!(a.len(), weight.len(), "values and weights must have equal length");
        let (core, arrangements) = WaveletMatrix::build(a);
        let prefix = arrangements
            .iter()
            .map(|arr| {
                let mut p = Vec::with_capacity(arr.len() + 1);
                let mut acc = W::default();
                p.push(acc);
                for &i in arr {
                    acc = acc + weight[i as usize];
                    p.push(acc);
                }
                p
            })
            .collect();
        Self { core, prefix }
    }

    /// Access to the unweighted queries.
    pub fn matrix(&self) -> &WaveletMatrix<T> {
        &self.core
    }

    /// `(count, total weight)` of the elements of `a[l..r]` whose rank is below `rank`.
    fn below_rank(&self, mut l: usize, mut r: usize, rank: usize) -> (usize, W) {
        let c = &self.core;
        if rank == 0 {
            return (0, W::default());
        }
        if rank >= c.vals.len() {
            return (r - l, self.prefix[0][r] - self.prefix[0][l]);
        }
        let depth = c.levels.len();
        let (mut count, mut sum) = (0usize, W::default());
        for (i, level) in c.levels.iter().enumerate() {
            let (ones_l, ones_r) = (level.column.rank1(l), level.column.rank1(r));
            let (zeros_l, zeros_r) = (l - ones_l, r - ones_r);
            if (rank >> (depth - 1 - i)) & 1 == 1 {
                count += zeros_r - zeros_l;
                // The zero side maps to [zeros_l, zeros_r) in the next arrangement.
                sum = sum + (self.prefix[i + 1][zeros_r] - self.prefix[i + 1][zeros_l]);
                l = level.zeros + ones_l;
                r = level.zeros + ones_r;
            } else {
                l = zeros_l;
                r = zeros_r;
            }
        }
        (count, sum)
    }

    /// `(count, sum of weights)` of the elements of `a[l..r]` with value `<= x`.
    pub fn count_sum_less_eq(&self, l: usize, r: usize, x: T) -> (usize, W) {
        self.core.check_range(l, r);
        self.below_rank(l, r, self.core.vals.partition_point(|v| *v <= x))
    }

    /// `(count, sum of weights)` of the elements of `a[l..r]` with value `< x`.
    pub fn count_sum_less(&self, l: usize, r: usize, x: T) -> (usize, W) {
        self.core.check_range(l, r);
        self.below_rank(l, r, self.core.vals.partition_point(|v| *v < x))
    }
}
