//! Determinant modulo an arbitrary integer.

/// `x mod m` for `x < 2^64` and `2 <= m < 2^32` by Barrett reduction.
#[derive(Clone, Copy)]
struct Barrett {
    m: u64,
    inv: u64,
}

impl Barrett {
    fn new(m: u64) -> Self {
        Barrett { m, inv: ((1u128 << 64) / m as u128) as u64 }
    }

    #[inline(always)]
    fn reduce(&self, x: u64) -> u64 {
        let q = ((x as u128 * self.inv as u128) >> 64) as u64;
        let r = x - q * self.m;
        if r >= self.m {
            r - self.m
        } else {
            r
        }
    }
}

/// The determinant of the `n x n` matrix `a` (row-major, entries below `m`) modulo `m`, for any
/// `1 <= m < 2^32`: Gaussian elimination where a column is cleared with Euclid's algorithm on
/// the rows instead of a division. `O(n^3 + n^2 log m)` operations.
pub fn determinant_mod(n: usize, a: &[u32], m: u32) -> u32 {
    assert_eq!(a.len(), n * n, "wrong number of entries");
    assert!(m >= 1, "the modulus must be positive");
    if m == 1 {
        return 0;
    }
    assert!(a.iter().all(|&x| x < m), "entries must be reduced modulo m");
    let bar = Barrett::new(m as u64);
    let mut a: Vec<u64> = a.iter().map(|&x| x as u64).collect();
    let mut negate = false;
    let mut det = 1u64;
    for c in 0..n {
        for r in c + 1..n {
            while a[r * n + c] != 0 {
                let q = a[c * n + c] / a[r * n + c];
                if q != 0 {
                    let (head, tail) = a.split_at_mut(r * n);
                    let (crow, rrow) = (&mut head[c * n + c..c * n + n], &tail[c..n]);
                    for (x, &y) in crow.iter_mut().zip(rrow) {
                        let t = bar.reduce(q * y);
                        *x = if *x >= t { *x - t } else { *x + m as u64 - t };
                    }
                }
                for j in c..n {
                    a.swap(c * n + j, r * n + j);
                }
                negate = !negate;
            }
        }
        det = bar.reduce(det * a[c * n + c]);
        if det == 0 {
            return 0;
        }
    }
    let det = if negate { (m as u64 - det) % m as u64 } else { det };
    det as u32
}
