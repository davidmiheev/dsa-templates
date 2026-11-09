//! Range "add an arithmetic progression" and range minimum.
//!
//! [`LinearAddRangeMin`] keeps an `i64` array and supports
//!
//! * [`add_linear`](LinearAddRangeMin::add_linear)`(l, r, b, c)`: `a_i += b * i + c` for `i` in
//!   `[l, r)` (`i` is the absolute index),
//! * [`min`](LinearAddRangeMin::min)`(l, r)`: `min(a[l..r])`.
//!
//! The array is cut into blocks of about `sqrt(n)` elements. A block lazily remembers the
//! accumulated `(b, c)` of the updates that covered it completely, and keeps the lower envelope
//! of the lines `x -> a_i + x * i`: its minimum under pending slope `x` is a lookup in that
//! envelope (binary search, or a short walk from the previous position). A partial update
//! pushes the lazies into the block, edits it and rebuilds the envelope in `O(block)`.
//! Each operation costs `O(sqrt(n) log n)` worst case and the structure works for updates
//! of any sign.
//!
//! All arithmetic is `i64`: the values `a_i + b_total * i + c_total` must fit in `i64`, and
//! the envelope construction uses `i128` internally.
//!
//! # Example
//!
//! ```
//! use linear_add_range_min::LinearAddRangeMin;
//!
//! let mut t = LinearAddRangeMin::new(&[5, 5, 5, 5]);
//! t.add_linear(0, 4, 2, 0); // [5, 7, 9, 11]
//! t.add_linear(1, 3, -3, 1); // adds -3*i + 1 to indices 1 and 2: [5, 5, 4, 11]
//! assert_eq!(t.min(0, 4), 4);
//! assert_eq!(t.min(3, 4), 11);
//! ```

/// See the [module documentation](self).
pub struct LinearAddRangeMin {
    n: usize,
    block: usize,
    /// Values without the block lazies.
    y: Vec<i64>,
    lazy_b: Vec<i64>,
    lazy_c: Vec<i64>,
    /// Block `k`'s envelope lines (indices `i`, slope descending) live in
    /// `hull[k * block .. k * block + hull_len[k]]`.
    hull: Vec<u32>,
    hull_len: Vec<u32>,
    /// Last position found in the envelope, to start the next search from.
    pos: Vec<u32>,
}

impl LinearAddRangeMin {
    /// Build over `a` in `O(n)`.
    ///
    /// # Panics
    /// Panics if `a` is empty or has `2^32` or more elements.
    pub fn new(a: &[i64]) -> Self {
        let block = ((a.len() as f64).sqrt() * 1.2).ceil().max(1.0) as usize;
        Self::with_block_size(a, block)
    }

    /// Like [`new`](Self::new) with an explicit block size (at least 1).
    pub fn with_block_size(a: &[i64], block: usize) -> Self {
        assert!(!a.is_empty(), "the array must not be empty");
        assert!(a.len() < u32::MAX as usize && block >= 1, "invalid size or block size");
        let blocks = a.len().div_ceil(block);
        let mut t = LinearAddRangeMin {
            n: a.len(),
            block,
            y: a.to_vec(),
            lazy_b: vec![0; blocks],
            lazy_c: vec![0; blocks],
            hull: vec![0; a.len()],
            hull_len: vec![0; blocks],
            pos: vec![0; blocks],
        };
        for k in 0..blocks {
            t.rebuild(k);
        }
        t
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    fn range_of(&self, k: usize) -> (usize, usize) {
        (k * self.block, ((k + 1) * self.block).min(self.n))
    }

    /// Rebuild the lower envelope of the lines `x -> y_i + x * i` of block `k`.
    fn rebuild(&mut self, k: usize) {
        let (lo, hi) = self.range_of(k);
        let base = k * self.block;
        let mut len = 0usize;
        // Slopes `i` in decreasing order: the envelope is traversed as `x` increases.
        for i in (lo..hi).rev() {
            while len >= 2 {
                let (l1, l2) = (self.hull[base + len - 2] as usize, self.hull[base + len - 1] as usize);
                let (m1, m2, m3) = (l1 as i128, l2 as i128, i as i128);
                let (c1, c2, c3) = (self.y[l1] as i128, self.y[l2] as i128, self.y[i] as i128);
                if (c3 - c1) * (m1 - m2) <= (c2 - c1) * (m1 - m3) {
                    len -= 1;
                } else {
                    break;
                }
            }
            self.hull[base + len] = i as u32;
            len += 1;
        }
        self.hull_len[k] = len as u32;
        self.pos[k] = 0;
    }

    /// `min_i (y_i + x * i)` over block `k`.
    fn block_min(&mut self, k: usize, x: i64) -> i64 {
        let base = k * self.block;
        let len = self.hull_len[k] as usize;
        let value = |t: &Self, p: usize| {
            let i = t.hull[base + p] as usize;
            t.y[i] + x * i as i64
        };
        let mut p = (self.pos[k] as usize).min(len - 1);
        // Walk first; the envelope is unimodal along `p`, so this is correct.
        let mut steps = 0;
        loop {
            if p + 1 < len && value(self, p + 1) <= value(self, p) {
                p += 1;
            } else if p > 0 && value(self, p - 1) <= value(self, p) {
                p -= 1;
            } else {
                break;
            }
            steps += 1;
            if steps == 8 {
                // Far away: binary search for the first p whose successor is not better.
                let (mut lo, mut hi) = (0, len - 1);
                while lo < hi {
                    let mid = (lo + hi) / 2;
                    if value(self, mid + 1) <= value(self, mid) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                p = lo;
                break;
            }
        }
        self.pos[k] = p as u32;
        value(self, p)
    }

    /// Push the lazies of block `k` into `y` (the envelope stays valid up to this shift only
    /// if rebuilt, which the caller does).
    fn push(&mut self, k: usize) {
        let (b, c) = (self.lazy_b[k], self.lazy_c[k]);
        if b != 0 || c != 0 {
            let (lo, hi) = self.range_of(k);
            for i in lo..hi {
                self.y[i] += b * i as i64 + c;
            }
            self.lazy_b[k] = 0;
            self.lazy_c[k] = 0;
        }
    }

    fn check(&self, l: usize, r: usize) {
        assert!(l <= r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
    }

    /// `a_i += b * i + c` for `i` in `[l, r)`.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn add_linear(&mut self, l: usize, r: usize, b: i64, c: i64) {
        self.check(l, r);
        if l == r {
            return;
        }
        let (kl, kr) = (l / self.block, (r - 1) / self.block);
        for k in kl..=kr {
            let (lo, hi) = self.range_of(k);
            if l <= lo && hi <= r {
                self.lazy_b[k] += b;
                self.lazy_c[k] += c;
            } else {
                self.push(k);
                for i in l.max(lo)..r.min(hi) {
                    self.y[i] += b * i as i64 + c;
                }
                self.rebuild(k);
            }
        }
    }

    /// Minimum of `a[l..r]`.
    ///
    /// # Panics
    /// Panics unless `l < r <= len()`.
    pub fn min(&mut self, l: usize, r: usize) -> i64 {
        assert!(l < r && r <= self.n, "invalid non-empty range [{l}, {r}) for length {}", self.n);
        let (kl, kr) = (l / self.block, (r - 1) / self.block);
        let mut best = i64::MAX;
        for k in kl..=kr {
            let (lo, hi) = self.range_of(k);
            if l <= lo && hi <= r {
                best = best.min(self.block_min(k, self.lazy_b[k]) + self.lazy_c[k]);
            } else {
                let (b, c) = (self.lazy_b[k], self.lazy_c[k]);
                for i in l.max(lo)..r.min(hi) {
                    best = best.min(self.y[i] + b * i as i64 + c);
                }
            }
        }
        best
    }

    /// Element `i`.
    pub fn get(&mut self, i: usize) -> i64 {
        assert!(i < self.n, "index {i} out of range for length {}", self.n);
        self.min(i, i + 1)
    }
}
