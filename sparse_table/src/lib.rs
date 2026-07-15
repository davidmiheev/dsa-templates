//! Sparse Table for idempotent range queries.
//!
//! A Sparse Table answers range queries on a static array in O(1) query
//! time after O(n log n) preprocessing, *as long as* the query operation
//! is idempotent and associative. The classic examples are `min`, `max`,
//! and `gcd`. (Sum works too; xor works too.)
//!
//! The trick: precompute `st[k][i] = op(a[i], a[i+1], ..., a[i+2^k-1])`
//! and overlap two intervals of length `2^k` to cover any range `[l, r]`.
//!
//! # Example
//!
//! ```
//! use sparse_table::SparseTable;
//!
//! let a = vec![1, 3, 2, 5, 4];
//! let st = SparseTable::new(&a, |x, y| if x < y { x } else { y });
//! assert_eq!(st.query(1, 4), Some(2)); // min over [1, 4) = min(3, 2, 5) = 2
//! assert_eq!(st.query(0, 5), Some(1)); // min over whole array
//! ```

/// Sparse Table for idempotent, associative operations on a static array.
pub struct SparseTable<T, F: Fn(T, T) -> T> {
    /// `st[k][i]` holds the result of `op` over the half-open window
    /// `[i, i + 2^k)` of the original array.
    st: Vec<Vec<T>>,
    /// Precomputed floor(log2) values for fast index computation.
    log: Vec<usize>,
    /// The user-supplied operation, stored for use in queries.
    op: F,
}

impl<T: Clone + Copy, F: Fn(T, T) -> T + Copy> SparseTable<T, F> {
    /// Build a Sparse Table from a slice using the given operation `op`,
    /// which must be associative and idempotent.
    pub fn new(a: &[T], op: F) -> Self {
        let n = a.len();
        if n == 0 {
            return SparseTable {
                st: Vec::new(),
                log: vec![0],
                op,
            };
        }

        let max_k = (usize::BITS - n.leading_zeros()) as usize;
        let mut st: Vec<Vec<T>> = Vec::with_capacity(max_k);
        st.push(a.to_vec());
        for k in 1..max_k {
            let prev = &st[k - 1];
            // Row k has entries for windows [i, i + 2^k) with i + 2^k <= n,
            // so the count is n - 2^k + 1 (which is >= 1 because k < max_k).
            let len = n + 1 - (1 << k);
            let mut row: Vec<T> = Vec::with_capacity(len);
            for i in 0..len {
                row.push(op(prev[i], prev[i + (1 << (k - 1))]));
            }
            st.push(row);
        }

        let mut log = vec![0usize; n + 1];
        for i in 2..=n {
            log[i] = log[i / 2] + 1;
        }

        SparseTable { st, log, op }
    }

    /// Query `op(a[l], a[l+1], ..., a[r-1])` over the half-open range
    /// `[l, r)`. Returns `None` if the range is empty or out of bounds.
    ///
    /// Two windows of length `2^k`, where `k = log2(r - l)`, exactly
    /// cover `[l, r)`. Since `op` is idempotent, overlapping them is
    /// safe.
    pub fn query(&self, l: usize, r: usize) -> Option<T> {
        let n = self.st.get(0).map(|v| v.len()).unwrap_or(0);
        if l >= r || r > n {
            return None;
        }
        let len = r - l;
        let k = self.log[len];
        let a = self.st[k][l];
        let b = self.st[k][r - (1 << k)];
        Some((self.op)(a, b))
    }
}

/// Convenience constructor for `min` queries.
pub fn min_sparse_table(a: &[i64]) -> SparseTable<i64, fn(i64, i64) -> i64> {
    SparseTable::new(a, |x, y| if x < y { x } else { y })
}

/// Convenience constructor for `max` queries.
pub fn max_sparse_table(a: &[i64]) -> SparseTable<i64, fn(i64, i64) -> i64> {
    SparseTable::new(a, |x, y| if x > y { x } else { y })
}

/// Convenience constructor for `gcd` queries.
pub fn gcd_sparse_table(a: &[i64]) -> SparseTable<i64, fn(i64, i64) -> i64> {
    SparseTable::new(a, |x, y| {
        let (mut a, mut b) = (x.abs(), y.abs());
        while b != 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        a
    })
}
