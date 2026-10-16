//! A flat row-major matrix of residues modulo a prime `M < 2^31`.

use crate::pow_u32;

/// `x * f mod m` for `f < M`, with `f_shoup = floor(f * 2^32 / M)` precomputed (Shoup's trick:
/// no division, vectorizes).
#[inline(always)]
fn mul_shoup<const M: u32>(x: u32, f: u32, f_shoup: u32) -> u32 {
    let q = ((x as u64 * f_shoup as u64) >> 32) as u32;
    let r = x.wrapping_mul(f).wrapping_sub(q.wrapping_mul(M));
    if r >= M {
        r - M
    } else {
        r
    }
}

#[inline(always)]
fn shoup<const M: u32>(f: u32) -> u32 {
    (((f as u64) << 32) / M as u64) as u32
}

/// `dst -= f * src` elementwise.
#[inline(always)]
fn axpy_sub<const M: u32>(dst: &mut [u32], src: &[u32], f: u32) {
    let fs = shoup::<M>(f);
    for (d, &s) in dst.iter_mut().zip(src) {
        let t = mul_shoup::<M>(s, f, fs);
        *d = if *d >= t { *d - t } else { *d + M - t };
    }
}

/// `row *= f`.
#[inline(always)]
fn scale<const M: u32>(row: &mut [u32], f: u32) {
    let fs = shoup::<M>(f);
    for x in row.iter_mut() {
        *x = mul_shoup::<M>(*x, f, fs);
    }
}

#[inline(always)]
fn mul_mod<const M: u32>(a: u32, b: u32) -> u32 {
    (a as u64 * b as u64 % M as u64) as u32
}

#[inline(always)]
fn sub_mod<const M: u32>(a: u32, b: u32) -> u32 {
    if a >= b {
        a - b
    } else {
        a + M - b
    }
}

#[inline(always)]
fn mul_kernel<const M: u32>(a: &[u32], b: &[u32], c: &mut [u32], n: usize, m: usize, k: usize) {
    const JT: usize = 256;
    let mm = M as u64;
    // number of products that can be added to a residue without overflowing a u64
    let chunk = ((u64::MAX - mm) / ((mm - 1) * (mm - 1)).max(1)) as usize;
    let mut acc = [0u64; JT];
    for jj in (0..k).step_by(JT) {
        let jw = (k - jj).min(JT);
        for i in 0..n {
            acc[..jw].fill(0);
            let mut cnt = 0;
            for (kk, &aik) in a[i * m..(i + 1) * m].iter().enumerate() {
                if aik == 0 {
                    continue;
                }
                let brow = &b[kk * k + jj..kk * k + jj + jw];
                for (x, &y) in acc[..jw].iter_mut().zip(brow) {
                    *x += aik as u64 * y as u64;
                }
                cnt += 1;
                if cnt == chunk {
                    for x in acc[..jw].iter_mut() {
                        *x %= mm;
                    }
                    cnt = 0;
                }
            }
            for (cj, &x) in c[i * k + jj..i * k + jj + jw].iter_mut().zip(&acc[..jw]) {
                *cj = (x % mm) as u32;
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn mul_kernel_avx2<const M: u32>(a: &[u32], b: &[u32], c: &mut [u32], n: usize, m: usize, k: usize) {
    mul_kernel::<M>(a, b, c, n, m, k)
}

/// A matrix with entries in `0..M`, `M` prime below `2^31`.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Matrix<const M: u32> {
    rows: usize,
    cols: usize,
    a: Vec<u32>,
}

impl<const M: u32> Matrix<M> {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Matrix { rows, cols, a: vec![0; rows * cols] }
    }

    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m.a[i * n + i] = 1;
        }
        m
    }

    /// A matrix from its entries in row-major order (all below `M`).
    pub fn from_vec(rows: usize, cols: usize, a: Vec<u32>) -> Self {
        assert_eq!(a.len(), rows * cols, "wrong number of entries");
        assert!(a.iter().all(|&x| x < M), "entries must be reduced modulo {M}");
        Matrix { rows, cols, a }
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn get(&self, i: usize, j: usize) -> u32 {
        assert!(i < self.rows && j < self.cols, "index out of range");
        self.a[i * self.cols + j]
    }

    pub fn set(&mut self, i: usize, j: usize, v: u32) {
        assert!(i < self.rows && j < self.cols && v < M, "index or value out of range");
        self.a[i * self.cols + j] = v;
    }

    pub fn row(&self, i: usize) -> &[u32] {
        &self.a[i * self.cols..(i + 1) * self.cols]
    }

    /// All entries in row-major order.
    pub fn data(&self) -> &[u32] {
        &self.a
    }

    pub fn transpose(&self) -> Matrix<M> {
        let mut t = Matrix::zeros(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                t.a[j * self.rows + i] = self.a[i * self.cols + j];
            }
        }
        t
    }

    /// The product `self * other`.
    pub fn mul(&self, other: &Matrix<M>) -> Matrix<M> {
        assert_eq!(self.cols, other.rows, "dimension mismatch");
        let (n, m, k) = (self.rows, self.cols, other.cols);
        let mut c = Matrix::zeros(n, k);
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx2") {
                // SAFETY: the CPU supports AVX2 (checked above).
                unsafe { mul_kernel_avx2::<M>(&self.a, &other.a, &mut c.a, n, m, k) };
                return c;
            }
        }
        mul_kernel::<M>(&self.a, &other.a, &mut c.a, n, m, k);
        c
    }

    /// `self^e` by repeated squaring (`e = 0` gives the identity).
    pub fn pow(&self, mut e: u64) -> Matrix<M> {
        assert_eq!(self.rows, self.cols, "power needs a square matrix");
        let mut result = Matrix::identity(self.rows);
        let mut base = self.clone();
        while e > 0 {
            if e & 1 == 1 {
                result = result.mul(&base);
            }
            e >>= 1;
            if e > 0 {
                base = base.mul(&base);
            }
        }
        result
    }

    /// The determinant (`1` for the empty matrix).
    pub fn det(&self) -> u32 {
        assert_eq!(self.rows, self.cols, "determinant needs a square matrix");
        let n = self.rows;
        let mut a = self.a.clone();
        let mut det = 1u32;
        for c in 0..n {
            let Some(p) = (c..n).find(|&r| a[r * n + c] != 0) else {
                return 0;
            };
            if p != c {
                for j in 0..n {
                    a.swap(p * n + j, c * n + j);
                }
                det = sub_mod::<M>(0, det);
            }
            let piv = a[c * n + c];
            det = mul_mod::<M>(det, piv);
            let inv = pow_u32::<M>(piv, M as u64 - 2);
            let (head, tail) = a.split_at_mut((c + 1) * n);
            let prow = &head[c * n + c..c * n + n];
            for row in tail.chunks_exact_mut(n) {
                let x = row[c];
                if x != 0 {
                    axpy_sub::<M>(&mut row[c..], prow, mul_mod::<M>(x, inv));
                }
            }
        }
        det
    }

    /// Gauss–Jordan on the first `pivot_cols` columns of the `rows x cols` array `a`; returns the
    /// pivot column of each pivot row.
    fn rref(a: &mut [u32], rows: usize, cols: usize, pivot_cols: usize) -> Vec<usize> {
        let mut pivots = Vec::new();
        let mut r = 0;
        for c in 0..pivot_cols {
            if r == rows {
                break;
            }
            let Some(p) = (r..rows).find(|&i| a[i * cols + c] != 0) else {
                continue;
            };
            if p != r {
                for j in 0..cols {
                    a.swap(p * cols + j, r * cols + j);
                }
            }
            let inv = pow_u32::<M>(a[r * cols + c], M as u64 - 2);
            scale::<M>(&mut a[r * cols + c..(r + 1) * cols], inv);
            let (above, rest) = a.split_at_mut(r * cols);
            let (prow, below) = rest.split_at_mut(cols);
            let prow = &prow[c..];
            for row in above.chunks_exact_mut(cols).chain(below.chunks_exact_mut(cols)) {
                let x = row[c];
                if x != 0 {
                    axpy_sub::<M>(&mut row[c..], prow, x);
                }
            }
            pivots.push(c);
            r += 1;
        }
        pivots
    }

    /// The rank.
    pub fn rank(&self) -> usize {
        // eliminate below the pivot only (no back substitution)
        let (rows, cols) = (self.rows, self.cols);
        let mut a = self.a.clone();
        let mut r = 0;
        for c in 0..cols {
            if r == rows {
                break;
            }
            let Some(p) = (r..rows).find(|&i| a[i * cols + c] != 0) else {
                continue;
            };
            if p != r {
                for j in 0..cols {
                    a.swap(p * cols + j, r * cols + j);
                }
            }
            let inv = pow_u32::<M>(a[r * cols + c], M as u64 - 2);
            let (head, tail) = a.split_at_mut((r + 1) * cols);
            let prow = &head[r * cols + c..];
            for row in tail.chunks_exact_mut(cols) {
                let x = row[c];
                if x != 0 {
                    axpy_sub::<M>(&mut row[c..], prow, mul_mod::<M>(x, inv));
                }
            }
            r += 1;
        }
        r
    }

    /// The inverse of a square matrix with `n >= 1`, or the `0 x 0` matrix if it is singular.
    pub fn inverse(&self) -> Matrix<M> {
        assert_eq!(self.rows, self.cols, "inverse needs a square matrix");
        let n = self.rows;
        assert!(n >= 1, "the matrix must be non-empty");
        let w = 2 * n;
        let mut aug = vec![0u32; n * w];
        for i in 0..n {
            aug[i * w..i * w + n].copy_from_slice(&self.a[i * n..(i + 1) * n]);
            aug[i * w + n + i] = 1;
        }
        if Self::rref(&mut aug, n, w, n).len() < n {
            return Matrix::zeros(0, 0);
        }
        let mut inv = Matrix::zeros(n, n);
        for i in 0..n {
            inv.a[i * n..(i + 1) * n].copy_from_slice(&aug[i * w + n..(i + 1) * w]);
        }
        inv
    }

    /// Solves `self * x = b`: returns `(particular, basis)` with one solution and a basis of the
    /// homogeneous solutions; `particular` is empty when there is no solution (`cols >= 1`).
    pub fn solve(&self, b: &[u32]) -> (Vec<u32>, Vec<Vec<u32>>) {
        assert!(self.cols >= 1, "the system needs at least one unknown");
        assert_eq!(b.len(), self.rows, "right-hand side length mismatch");
        let (n, m) = (self.rows, self.cols);
        let w = m + 1;
        let mut aug = vec![0u32; n * w];
        for i in 0..n {
            aug[i * w..i * w + m].copy_from_slice(&self.a[i * m..(i + 1) * m]);
            aug[i * w + m] = b[i];
        }
        let pivots = Self::rref(&mut aug, n, w, m);
        if (pivots.len()..n).any(|i| aug[i * w + m] != 0) {
            return (Vec::new(), Vec::new());
        }
        let mut particular = vec![0u32; m];
        for (r, &c) in pivots.iter().enumerate() {
            particular[c] = aug[r * w + m];
        }
        let mut is_pivot = vec![false; m];
        for &c in &pivots {
            is_pivot[c] = true;
        }
        let mut basis = Vec::new();
        for f in (0..m).filter(|&f| !is_pivot[f]) {
            let mut v = vec![0u32; m];
            v[f] = 1;
            for (r, &c) in pivots.iter().enumerate() {
                v[c] = sub_mod::<M>(0, aug[r * w + f]);
            }
            basis.push(v);
        }
        (particular, basis)
    }

    /// The adjugate (`adj(A)_{ij} = (-1)^{i+j} M_{ji}`), also for singular matrices, in `O(n^3)`.
    pub fn adjugate(&self) -> Matrix<M> {
        assert_eq!(self.rows, self.cols, "adjugate needs a square matrix");
        let n = self.rows;
        if n == 1 {
            return Matrix::from_vec(1, 1, vec![1]);
        }
        let rank = self.rank();
        if rank == n {
            let det = self.det();
            let mut adj = self.inverse();
            for x in adj.a.iter_mut() {
                *x = mul_mod::<M>(*x, det);
            }
            return adj;
        }
        if rank + 2 <= n {
            return Matrix::zeros(n, n);
        }
        // rank n - 1: adj = c * x y^T with x spanning ker A and y spanning ker A^T
        let zero = vec![0u32; n];
        let x = self.solve(&zero).1.pop().expect("a rank n-1 matrix has a kernel");
        let y = self.transpose().solve(&zero).1.pop().expect("a rank n-1 matrix has a left kernel");
        let i = x.iter().position(|&v| v != 0).expect("kernel vectors are nonzero");
        let j = y.iter().position(|&v| v != 0).expect("kernel vectors are nonzero");
        // adj[i][j] = (-1)^{i+j} det(A without row j and column i)
        let mut minor = Vec::with_capacity((n - 1) * (n - 1));
        for r in (0..n).filter(|&r| r != j) {
            for c in (0..n).filter(|&c| c != i) {
                minor.push(self.a[r * n + c]);
            }
        }
        let mut d = Matrix::<M>::from_vec(n - 1, n - 1, minor).det();
        if (i + j) % 2 == 1 {
            d = sub_mod::<M>(0, d);
        }
        let c = mul_mod::<M>(d, pow_u32::<M>(mul_mod::<M>(x[i], y[j]), M as u64 - 2));
        let mut adj = Matrix::zeros(n, n);
        for (r, &xr) in x.iter().enumerate() {
            let cx = mul_mod::<M>(c, xr);
            for (s, &ys) in y.iter().enumerate() {
                adj.a[r * n + s] = mul_mod::<M>(cx, ys);
            }
        }
        adj
    }

    /// The coefficients `p_0, ..., p_n` of `det(x I - A)` (monic), via the Hessenberg form,
    /// `O(n^3)`.
    pub fn characteristic_polynomial(&self) -> Vec<u32> {
        assert_eq!(self.rows, self.cols, "characteristic polynomial needs a square matrix");
        let n = self.rows;
        let mut h = self.a.clone();
        // similarity transformations to upper Hessenberg form
        for j in 0..n.saturating_sub(2) {
            let Some(p) = (j + 1..n).find(|&i| h[i * n + j] != 0) else {
                continue;
            };
            if p != j + 1 {
                for c in 0..n {
                    h.swap(p * n + c, (j + 1) * n + c);
                }
                for r in 0..n {
                    h.swap(r * n + p, r * n + j + 1);
                }
            }
            let inv = pow_u32::<M>(h[(j + 1) * n + j], M as u64 - 2);
            for i in j + 2..n {
                let f = mul_mod::<M>(h[i * n + j], inv);
                if f == 0 {
                    continue;
                }
                // row_i -= f row_{j+1}
                let (head, tail) = h.split_at_mut(i * n);
                axpy_sub::<M>(&mut tail[j..n], &head[(j + 1) * n + j..(j + 2) * n], f);
                // col_{j+1} += f col_i
                for r in 0..n {
                    let t = mul_mod::<M>(f, h[r * n + i]);
                    let v = &mut h[r * n + j + 1];
                    *v = if *v + t >= M { *v + t - M } else { *v + t };
                }
            }
        }
        // p_{k+1} = (x - h_kk) p_k - sum_{i<k} h_ik (prod_{j=i+1}^{k} h_{j,j-1}) p_i
        let mut p: Vec<Vec<u32>> = vec![vec![1]];
        for k in 0..n {
            let mut next = vec![0u32; k + 2];
            for d in 0..=k {
                next[d + 1] = p[k][d];
                next[d] = sub_mod::<M>(next[d], mul_mod::<M>(h[k * n + k], p[k][d]));
            }
            let mut t = 1u32;
            for i in (0..k).rev() {
                t = mul_mod::<M>(t, h[(i + 1) * n + i]);
                let coef = mul_mod::<M>(t, h[i * n + k]);
                if coef != 0 {
                    for d in 0..=i {
                        next[d] = sub_mod::<M>(next[d], mul_mod::<M>(coef, p[i][d]));
                    }
                }
            }
            p.push(next);
        }
        p.pop().expect("n + 1 polynomials were built")
    }

    /// The Pfaffian of an alternating matrix of even order (`A^T = -A`, zero diagonal), `O(n^3)`.
    pub fn pfaffian(&self) -> u32 {
        assert_eq!(self.rows, self.cols, "the Pfaffian needs a square matrix");
        let n = self.rows;
        assert!(n.is_multiple_of(2), "the Pfaffian needs an even order");
        for i in 0..n {
            for j in 0..n {
                assert_eq!(self.a[i * n + j], sub_mod::<M>(0, self.a[j * n + i]), "the matrix is not alternating");
            }
        }
        let mut a = self.a.clone();
        let mut pf = 1u32;
        for k in (0..n).step_by(2) {
            let Some(p) = (k + 1..n).find(|&p| a[k * n + p] != 0) else {
                return 0;
            };
            if p != k + 1 {
                for c in 0..n {
                    a.swap(p * n + c, (k + 1) * n + c);
                }
                for r in 0..n {
                    a.swap(r * n + p, r * n + k + 1);
                }
                pf = sub_mod::<M>(0, pf);
            }
            let piv = a[k * n + k + 1];
            pf = mul_mod::<M>(pf, piv);
            let inv = pow_u32::<M>(piv, M as u64 - 2);
            // congruence by I - f_j e_j e_{k+1}^T for j >= k+2 clears row k:
            // a[j][l] += f_l b_j - f_j b_l with b = row k+1 and f_j = a[k][j] / a[k][k+1]
            let f: Vec<u32> = (0..n).map(|j| if j >= k + 2 { mul_mod::<M>(a[k * n + j], inv) } else { 0 }).collect();
            let b: Vec<u32> = (0..n).map(|j| a[(k + 1) * n + j]).collect();
            for j in k + 2..n {
                for l in k + 2..n {
                    let t = sub_mod::<M>(mul_mod::<M>(f[l], b[j]), mul_mod::<M>(f[j], b[l]));
                    let v = &mut a[j * n + l];
                    *v = if *v + t >= M { *v + t - M } else { *v + t };
                }
            }
        }
        pf
    }
}
