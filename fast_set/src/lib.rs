//! A set of integers in `0..n` as a 64-ary bitset tree (a van Emde Boas-like layout).
//!
//! Level 0 is a plain bitset; every higher level has one bit per 64-bit word of the level below
//! ("this word is non-empty"). Successor and predecessor queries climb until a non-empty word
//! is found and then descend with `trailing_zeros` / `leading_zeros`, so every operation takes
//! `O(log_64 n)`: at most 4 words for `n = 10^7`. Memory is `n / 64 * (1 + 1/64 + ...)` words.
//!
//! Queries that may have no answer return [`NONE`] (`usize::MAX`) instead of an `Option`.
//!
//! # Example
//!
//! ```
//! use fast_set::{FastSet, NONE};
//!
//! let mut s = FastSet::new(1000);
//! s.insert(10);
//! s.insert(500);
//! assert_eq!(s.next(11), 500); // smallest element >= 11
//! assert_eq!(s.prev(499), 10); // largest element <= 499
//! assert_eq!(s.prev(9), NONE);
//! s.remove(500);
//! assert_eq!(s.next(11), NONE);
//! ```

/// Returned by [`FastSet::next`] and [`FastSet::prev`] when there is no such element.
pub const NONE: usize = usize::MAX;

/// A set of integers in `0..n` with `O(log_64 n)` operations.
#[derive(Clone, Debug)]
pub struct FastSet {
    n: usize,
    levels: Vec<Vec<u64>>,
}

impl FastSet {
    /// An empty set over `0..n`.
    pub fn new(n: usize) -> Self {
        let mut levels = Vec::new();
        let mut size = n;
        loop {
            size = size.div_ceil(64);
            levels.push(vec![0u64; size.max(1)]);
            if size <= 1 {
                break;
            }
        }
        FastSet { n, levels }
    }

    /// The universe size.
    pub fn capacity(&self) -> usize {
        self.n
    }

    pub fn contains(&self, i: usize) -> bool {
        assert!(i < self.n, "element {i} out of range 0..{}", self.n);
        (self.levels[0][i >> 6] >> (i & 63)) & 1 == 1
    }

    /// Insert `i` (no effect if present).
    pub fn insert(&mut self, mut i: usize) {
        assert!(i < self.n, "element {i} out of range 0..{}", self.n);
        for level in &mut self.levels {
            level[i >> 6] |= 1 << (i & 63);
            i >>= 6;
        }
    }

    /// Remove `i` (no effect if absent).
    pub fn remove(&mut self, mut i: usize) {
        assert!(i < self.n, "element {i} out of range 0..{}", self.n);
        for level in &mut self.levels {
            level[i >> 6] &= !(1 << (i & 63));
            if level[i >> 6] != 0 {
                break;
            }
            i >>= 6;
        }
    }

    /// The smallest element `>= i`, or [`NONE`].
    pub fn next(&self, mut i: usize) -> usize {
        if i >= self.n {
            return NONE;
        }
        for h in 0..self.levels.len() {
            if (i >> 6) >= self.levels[h].len() {
                break;
            }
            let word = self.levels[h][i >> 6] >> (i & 63);
            if word == 0 {
                i = (i >> 6) + 1;
                continue;
            }
            i += word.trailing_zeros() as usize;
            for g in (0..h).rev() {
                i <<= 6;
                i += self.levels[g][i >> 6].trailing_zeros() as usize;
            }
            return i;
        }
        NONE
    }

    /// The largest element `<= i`, or [`NONE`].
    pub fn prev(&self, i: usize) -> usize {
        assert!(i < self.n, "element {i} out of range 0..{}", self.n);
        let mut i = i as isize;
        for h in 0..self.levels.len() {
            if i < 0 {
                break;
            }
            let word = self.levels[h][(i >> 6) as usize] << (63 - (i & 63));
            if word == 0 {
                i = (i >> 6) - 1;
                continue;
            }
            i -= word.leading_zeros() as isize;
            for g in (0..h).rev() {
                i <<= 6;
                i += 63 - self.levels[g][(i >> 6) as usize].leading_zeros() as isize;
            }
            return i as usize;
        }
        NONE
    }
}
