//! Sparse Table for idempotent range queries.
//!
//! A Sparse Table answers range queries on a static array in O(1) query time after O(n)
//! preprocessing, *as long as* the operation is associative and idempotent
//! (`op(x, x) == x`). The classic examples are `min`, `max`, `gcd`, bitwise `and`/`or`.
//!
//! # Layout
//!
//! The textbook table `st[k][i] = op(a[i..i + 2^k])` needs `n log n` memory, which at
//! `n = 5 * 10^5` is ~75 MB and makes every query two cache misses. This implementation
//! splits the array into blocks of [`BLOCK`] elements and keeps
//!
//! * `prefix[i]` / `suffix[i]`: the aggregate from the block start up to `i` / from `i` to the
//!   block end,
//! * a textbook sparse table over the *block* aggregates (`n / BLOCK` entries per level).
//!
//! A query spanning several blocks combines `suffix[l]`, one overlapped pair of table entries
//! for the blocks strictly between, and `prefix[r - 1]`. A query inside one block folds at most
//! `BLOCK` consecutive elements. Memory is about `3n`, queries stay O(1).
//!
//! # Example
//!
//! ```
//! use sparse_table::SparseTable;
//!
//! let a = vec![1, 3, 2, 5, 4];
//! let st = SparseTable::new(&a, |x, y| if x < y { x } else { y });
//! assert_eq!(st.query(1, 4), 2); // min over [1, 4) = min(3, 2, 5)
//! assert_eq!(st.query(0, 5), 1); // min over the whole array
//! ```

/// Number of elements per block (a power of two).
pub const BLOCK: usize = 32;
const SHIFT: u32 = BLOCK.trailing_zeros();

/// Sparse Table for idempotent, associative operations on a static array.
pub struct SparseTable<T, F> {
    a: Vec<T>,
    /// `prefix[i]`: `op` over `a[block_start(i)..=i]`.
    prefix: Vec<T>,
    /// `suffix[i]`: `op` over `a[i..block_end(i)]`.
    suffix: Vec<T>,
    /// `levels[k][b]`: `op` over the aggregates of blocks `b..b + 2^k`.
    levels: Vec<Vec<T>>,
    op: F,
}

impl<T: Copy, F: Fn(T, T) -> T + Copy> SparseTable<T, F> {
    /// Build a Sparse Table from a slice using `op`, which must be associative and idempotent.
    pub fn new(a: &[T], op: F) -> Self {
        let n = a.len();
        let mut prefix = a.to_vec();
        let mut suffix = a.to_vec();
        for block in prefix.chunks_mut(BLOCK) {
            for i in 1..block.len() {
                block[i] = op(block[i - 1], block[i]);
            }
        }
        for block in suffix.chunks_mut(BLOCK) {
            for i in (0..block.len().saturating_sub(1)).rev() {
                block[i] = op(block[i], block[i + 1]);
            }
        }
        // Aggregate of each block = suffix at the block start.
        let agg: Vec<T> = suffix.iter().step_by(BLOCK).copied().collect();
        let mut levels = vec![agg];
        let mut width = 1;
        while 2 * width <= levels[0].len() {
            let prev = levels.last().unwrap();
            // Level k+1 covers windows of 2 * width blocks: `prev.len() - width` of them.
            let next: Vec<T> = (0..prev.len() - width).map(|b| op(prev[b], prev[b + width])).collect();
            levels.push(next);
            width *= 2;
        }
        let _ = n;
        SparseTable { a: a.to_vec(), prefix, suffix, levels, op }
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.a.len()
    }

    pub fn is_empty(&self) -> bool {
        self.a.is_empty()
    }

    /// `op(a[l], a[l + 1], ..., a[r - 1])` over the half-open range `[l, r)`.
    ///
    /// # Panics
    /// Panics unless `l < r <= len()`.
    pub fn query(&self, l: usize, r: usize) -> T {
        assert!(l < r && r <= self.a.len(), "invalid range [{l}, {r}) for length {}", self.a.len());
        let (bl, br) = (l >> SHIFT, (r - 1) >> SHIFT);
        if bl == br {
            // Inside one block: fold at most BLOCK consecutive elements. Starting from
            // op(x, x) (= x for idempotent ops) keeps e.g. gcd(-x) = |x| for a lone element.
            let mut acc = (self.op)(self.a[l], self.a[l]);
            for &x in &self.a[l + 1..r] {
                acc = (self.op)(acc, x);
            }
            return acc;
        }
        let mut acc = self.suffix[l];
        if bl + 1 < br {
            // Blocks bl+1 .. br-1: two overlapping windows of 2^k blocks.
            let len = br - bl - 1;
            let k = (usize::BITS - 1 - len.leading_zeros()) as usize;
            let level = &self.levels[k];
            acc = (self.op)(acc, (self.op)(level[bl + 1], level[br - (1 << k)]));
        }
        (self.op)(acc, self.prefix[r - 1])
    }
}

/// `min` sparse table over any ordered type.
pub fn min_sparse_table<T: Ord + Copy>(a: &[T]) -> SparseTable<T, impl Fn(T, T) -> T + Copy> {
    SparseTable::new(a, |x, y| if y < x { y } else { x })
}

/// `max` sparse table over any ordered type.
pub fn max_sparse_table<T: Ord + Copy>(a: &[T]) -> SparseTable<T, impl Fn(T, T) -> T + Copy> {
    SparseTable::new(a, |x, y| if y > x { y } else { x })
}

/// `gcd` sparse table over `i64` (signs are ignored, `gcd(x, x) = |x|`).
pub fn gcd_sparse_table(a: &[i64]) -> SparseTable<i64, impl Fn(i64, i64) -> i64 + Copy> {
    SparseTable::new(a, |x, y| {
        let (mut a, mut b) = (x.abs(), y.abs());
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    })
}
