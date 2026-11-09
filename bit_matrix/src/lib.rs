//! Matrices over `F_2` with the rows packed into `u64` words.
//!
//! * [`BitMatrix::mul`]: product by the method of four Russians, `O(n^3 / (8 * 64))`.
//! * [`BitMatrix::rank`], [`BitMatrix::det`], [`BitMatrix::inverse`], [`BitMatrix::solve`]:
//!   Gauss–Jordan elimination with word-parallel row operations, `O(n^3 / 64)`.
//! * [`intersect_spaces`]: a basis of the intersection of two subspaces of `F_2^30` given by
//!   bases (Zassenhaus' algorithm on `u64`s).
//!
//! Vectors are `Vec<u64>` of `ceil(len / 64)` words, bit `j` in word `j / 64`.
//!
//! # Example
//!
//! ```
//! use bit_matrix::BitMatrix;
//!
//! let mut a = BitMatrix::new(2, 2);
//! a.set(0, 0, true);
//! a.set(0, 1, true);
//! a.set(1, 1, true);
//! assert!(a.det());
//! let inv = a.inverse();
//! assert_eq!(a.mul(&inv), BitMatrix::identity(2));
//! ```

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct BitMatrix {
    rows: usize,
    cols: usize,
    wpr: usize,
    data: Vec<u64>,
}

impl BitMatrix {
    /// The zero matrix.
    pub fn new(rows: usize, cols: usize) -> Self {
        let wpr = cols.div_ceil(64);
        BitMatrix { rows, cols, wpr, data: vec![0; rows * wpr] }
    }

    pub fn identity(n: usize) -> Self {
        let mut m = Self::new(n, n);
        for i in 0..n {
            m.set(i, i, true);
        }
        m
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn get(&self, i: usize, j: usize) -> bool {
        assert!(i < self.rows && j < self.cols, "index out of range");
        self.data[i * self.wpr + j / 64] >> (j % 64) & 1 == 1
    }

    pub fn set(&mut self, i: usize, j: usize, v: bool) {
        assert!(i < self.rows && j < self.cols, "index out of range");
        let w = &mut self.data[i * self.wpr + j / 64];
        if v {
            *w |= 1 << (j % 64);
        } else {
            *w &= !(1 << (j % 64));
        }
    }

    /// Row `i` as `ceil(cols / 64)` words.
    pub fn row(&self, i: usize) -> &[u64] {
        &self.data[i * self.wpr..(i + 1) * self.wpr]
    }

    /// Mutable access to the words of row `i` (bits beyond `cols` must stay zero).
    pub fn row_mut_words(&mut self, i: usize) -> &mut [u64] {
        &mut self.data[i * self.wpr..(i + 1) * self.wpr]
    }

    /// Sets row `i` from `cols` ASCII characters `'0'` / `'1'`.
    pub fn set_row_ascii(&mut self, i: usize, s: &[u8]) {
        assert_eq!(s.len(), self.cols, "row length mismatch");
        let row = &mut self.data[i * self.wpr..(i + 1) * self.wpr];
        row.fill(0);
        for (chunk, w) in s.chunks(64).zip(row.iter_mut()) {
            let mut x = 0u64;
            for (b, &ch) in chunk.iter().enumerate() {
                x |= ((ch & 1) as u64) << b;
            }
            *w = x;
        }
    }

    /// Appends row `i` as `cols` ASCII characters `'0'` / `'1'`.
    pub fn write_row_ascii(&self, i: usize, out: &mut Vec<u8>) {
        for j in 0..self.cols {
            out.push(b'0' + (self.data[i * self.wpr + j / 64] >> (j % 64) & 1) as u8);
        }
    }

    pub fn transpose(&self) -> BitMatrix {
        let mut t = BitMatrix::new(self.cols, self.rows);
        for i in 0..self.rows {
            for (wi, &w) in self.row(i).iter().enumerate() {
                let mut x = w;
                while x != 0 {
                    let j = wi * 64 + x.trailing_zeros() as usize;
                    x &= x - 1;
                    t.data[j * t.wpr + i / 64] |= 1 << (i % 64);
                }
            }
        }
        t
    }

    /// `self * other` by the method of four Russians (tables of the 256 combinations of 8 rows
    /// of `other`).
    pub fn mul(&self, other: &BitMatrix) -> BitMatrix {
        assert_eq!(self.cols, other.rows, "dimension mismatch");
        let mut c = BitMatrix::new(self.rows, other.cols);
        let w = other.wpr;
        let mut table = vec![0u64; 256 * w];
        for g in 0..self.cols.div_ceil(8) {
            let base = g * 8;
            let width = (self.cols - base).min(8);
            // table[x] = xor of the rows base + t of `other` for the bits t of x
            for x in 1usize..1 << width {
                let low = x.trailing_zeros() as usize;
                let (prev, src) = (x & (x - 1), base + low);
                for k in 0..w {
                    table[x * w + k] = table[prev * w + k] ^ other.data[src * w + k];
                }
            }
            for i in 0..self.rows {
                let bits = (self.data[i * self.wpr + base / 64] >> (base % 64)) as usize & ((1 << width) - 1);
                if bits != 0 {
                    let dst = &mut c.data[i * w..(i + 1) * w];
                    for (d, &t) in dst.iter_mut().zip(&table[bits * w..(bits + 1) * w]) {
                        *d ^= t;
                    }
                }
            }
        }
        c
    }

    fn xor_rows(&mut self, dst: usize, src: usize, from_word: usize) {
        debug_assert!(dst != src);
        let w = self.wpr;
        let (a, b) = if dst < src {
            let (lo, hi) = self.data.split_at_mut(src * w);
            (&mut lo[dst * w..(dst + 1) * w], &hi[..w])
        } else {
            let (lo, hi) = self.data.split_at_mut(dst * w);
            (&mut hi[..w], &lo[src * w..(src + 1) * w])
        };
        for (x, &y) in a[from_word..].iter_mut().zip(&b[from_word..]) {
            *x ^= y;
        }
    }

    fn swap_rows(&mut self, a: usize, b: usize) {
        if a != b {
            for k in 0..self.wpr {
                self.data.swap(a * self.wpr + k, b * self.wpr + k);
            }
        }
    }

    fn bit(&self, i: usize, j: usize) -> bool {
        self.data[i * self.wpr + j / 64] >> (j % 64) & 1 == 1
    }

    /// Gauss–Jordan on the first `ncols` columns (eliminating above and below when `full`);
    /// returns the pivot column of each pivot row, in order.
    fn eliminate(&mut self, ncols: usize, full: bool) -> Vec<usize> {
        let mut pivots = Vec::new();
        let mut r = 0;
        for c in 0..ncols {
            if r == self.rows {
                break;
            }
            let Some(p) = (r..self.rows).find(|&i| self.bit(i, c)) else {
                continue;
            };
            self.swap_rows(r, p);
            let start = if full { 0 } else { r + 1 };
            for i in start..self.rows {
                if i != r && self.bit(i, c) {
                    self.xor_rows(i, r, c / 64);
                }
            }
            pivots.push(c);
            r += 1;
        }
        pivots
    }

    /// The rank.
    pub fn rank(&self) -> usize {
        // eliminate along the shorter dimension's rows: work with the orientation with fewer rows
        let mut m = if self.rows <= self.cols { self.clone() } else { self.transpose() };
        let nc = m.cols;
        m.eliminate(nc, false).len()
    }

    /// The determinant (1 iff invertible).
    pub fn det(&self) -> bool {
        assert_eq!(self.rows, self.cols, "determinant needs a square matrix");
        self.rank() == self.rows
    }

    /// The inverse of a square matrix with `n >= 1`, or the `0 x 0` matrix if it is singular.
    pub fn inverse(&self) -> BitMatrix {
        assert_eq!(self.rows, self.cols, "inverse needs a square matrix");
        let n = self.rows;
        assert!(n >= 1, "the matrix must be non-empty");
        let mut aug = BitMatrix::new(n, 2 * n);
        for i in 0..n {
            let (src, dst) = (self.row(i), i * aug.wpr);
            // copy the left half word by word, then set the identity bit
            aug.data[dst..dst + src.len()].copy_from_slice(src);
            let j = n + i;
            aug.data[dst + j / 64] |= 1 << (j % 64);
        }
        if aug.eliminate(n, true).len() < n {
            return BitMatrix::new(0, 0);
        }
        let mut inv = BitMatrix::new(n, n);
        for i in 0..n {
            for j in 0..n {
                if aug.bit(i, n + j) {
                    inv.set(i, j, true);
                }
            }
        }
        inv
    }

    /// Solves `self * x = b` (`b` has `rows` bits). Returns `(particular, basis)`: one solution
    /// and a basis of the solutions of the homogeneous system, as `cols`-bit vectors; the
    /// particular solution is empty when the system is inconsistent (`cols >= 1`).
    pub fn solve(&self, b: &[u64]) -> (Vec<u64>, Vec<Vec<u64>>) {
        assert!(self.cols >= 1, "the system needs at least one unknown");
        assert_eq!(b.len(), self.rows.div_ceil(64), "right-hand side length mismatch");
        let (n, m) = (self.rows, self.cols);
        let mut aug = BitMatrix::new(n, m + 1);
        for i in 0..n {
            let row = self.row(i);
            let dst = i * aug.wpr;
            aug.data[dst..dst + row.len()].copy_from_slice(row);
            if b[i / 64] >> (i % 64) & 1 == 1 {
                aug.data[dst + m / 64] |= 1 << (m % 64);
            }
        }
        let pivots = aug.eliminate(m, true);
        if (pivots.len()..n).any(|i| aug.bit(i, m)) {
            return (Vec::new(), Vec::new());
        }
        let words = m.div_ceil(64);
        let mut particular = vec![0u64; words];
        for (r, &c) in pivots.iter().enumerate() {
            if aug.bit(r, m) {
                particular[c / 64] |= 1 << (c % 64);
            }
        }
        let mut is_pivot = vec![false; m];
        for &c in &pivots {
            is_pivot[c] = true;
        }
        let mut basis = Vec::new();
        for f in (0..m).filter(|&f| !is_pivot[f]) {
            let mut v = vec![0u64; words];
            v[f / 64] |= 1 << (f % 64);
            for (r, &c) in pivots.iter().enumerate() {
                if aug.bit(r, f) {
                    v[c / 64] |= 1 << (c % 64);
                }
            }
            basis.push(v);
        }
        (particular, basis)
    }
}

/// A basis of `span(u) ∩ span(v)` where `u` and `v` are each linearly independent over `F_2`
/// and the entries are below `2^32` (Zassenhaus: eliminate the rows `(u_i | u_i)` and `(v_j | 0)`
/// on their left halves; the rows whose left half vanishes carry the intersection).
pub fn intersect_spaces(u: &[u32], v: &[u32]) -> Vec<u32> {
    // pivot[b] holds a reduced row whose highest set bit of the left half is b
    let mut pivot = [0u64; 32];
    let mut out = Vec::new();
    let rows = u.iter().map(|&x| (x as u64) << 32 | x as u64).chain(v.iter().map(|&y| (y as u64) << 32));
    for mut w in rows {
        while w >> 32 != 0 {
            let b = 63 - w.leading_zeros() as usize - 32;
            if pivot[b] == 0 {
                pivot[b] = w;
                break;
            }
            w ^= pivot[b];
        }
        if w >> 32 == 0 {
            out.push(w as u32);
        }
    }
    out
}
