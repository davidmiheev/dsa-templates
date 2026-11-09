//! The hafnian by the subset recursion over pairs of vertices (Björklund et al.).

use crate::Matrix;

struct State {
    deg: usize,
    /// `mats[level]` holds the lower-triangular matrix of polynomials of that recursion level
    mats: Vec<Vec<u32>>,
}

#[inline(always)]
fn tri(i: usize) -> usize {
    i * i.wrapping_sub(1) / 2
}

impl State {
    /// `out[k+1] += sum_{a+b=k} (p[a] q[b] + r[a] s[b])` for the polynomials of length `deg`.
    #[inline(always)]
    fn add_products<const M: u32>(deg: usize, out: &mut [u32], p: &[u32], q: &[u32], rs: Option<(&[u32], &[u32])>) {
        let k_lim = 8 * (M as u64) * (M as u64);
        for k in 0..deg - 1 {
            let mut acc = 0u64;
            for a in 0..=k {
                acc += p[a] as u64 * q[k - a] as u64;
                acc = acc.min(acc.wrapping_sub(k_lim));
            }
            if let Some((r, s)) = rs {
                for a in 0..=k {
                    acc += r[a] as u64 * s[k - a] as u64;
                    acc = acc.min(acc.wrapping_sub(k_lim));
                }
            }
            out[k + 1] = ((out[k + 1] as u64 + acc) % M as u64) as u32;
        }
    }

    /// The polynomial of the matrix held in `mats[level][..tri(size)]`.
    fn f<const M: u32>(&mut self, level: usize, size: usize) -> Vec<u32> {
        let deg = self.deg;
        let mut ret = vec![0u32; deg];
        if size == 0 {
            ret[0] = 1;
            return ret;
        }
        let n = size - 2;
        let mut cur = std::mem::take(&mut self.mats[level]);
        let prefix = tri(n) * deg;
        // exclude the pair
        self.mats[level + 1][..prefix].copy_from_slice(&cur[..prefix]);
        let zero = self.f::<M>(level + 1, n);
        // include it: merge the two vertices into every pair of the others
        {
            let (top, rest) = cur.split_at_mut(prefix);
            let (y, x) = rest.split_at(n * deg);
            for i in 0..n {
                let (xi, yi) = (&x[i * deg..(i + 1) * deg], &y[i * deg..(i + 1) * deg]);
                for j in 0..i {
                    let (xj, yj) = (&x[j * deg..(j + 1) * deg], &y[j * deg..(j + 1) * deg]);
                    let cell = &mut top[(tri(i) + j) * deg..(tri(i) + j + 1) * deg];
                    Self::add_products::<M>(deg, cell, xi, yj, Some((yi, xj)));
                }
            }
        }
        self.mats[level + 1][..prefix].copy_from_slice(&cur[..prefix]);
        let all = self.f::<M>(level + 1, n);
        let xn = &cur[prefix + n * deg + n * deg..prefix + n * deg + (n + 1) * deg];
        Self::add_products::<M>(deg, &mut ret, xn, &all, None);
        for k in 0..deg {
            let v = ret[k] as u64 + all[k] as u64 + M as u64 - zero[k] as u64;
            ret[k] = (v % M as u64) as u32;
        }
        self.mats[level] = cur;
        ret
    }
}

/// The hafnian of a symmetric matrix of even order `n` with zero diagonal, i.e. the sum over the
/// perfect matchings of the product of the edge weights, in `O(2^{n/2} n^2 (n/2)^2)`.
/// The hafnian of the empty matrix is 1.
pub fn hafnian<const M: u32>(a: &Matrix<M>) -> u32 {
    assert!((M as u64) < 1 << 30, "the modulus must be below 2^30");
    let n = a.rows();
    assert_eq!(a.cols(), n, "the hafnian needs a square matrix");
    assert!(n.is_multiple_of(2), "the hafnian needs an even order");
    for i in 0..n {
        assert_eq!(a.get(i, i), 0, "the diagonal must be zero");
        for j in 0..i {
            assert_eq!(a.get(i, j), a.get(j, i), "the matrix must be symmetric");
        }
    }
    if n == 0 {
        return 1;
    }
    let deg = n / 2 + 1;
    let mut mats: Vec<Vec<u32>> = (0..=n / 2).map(|d| vec![0u32; tri(n - 2 * d) * deg]).collect();
    for i in 0..n {
        for j in 0..i {
            mats[0][(tri(i) + j) * deg] = a.get(i, j);
        }
    }
    let mut st = State { deg, mats };
    st.f::<M>(0, n)[deg - 1]
}
