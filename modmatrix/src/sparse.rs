//! Wiedemann's algorithm for the determinant of a sparse matrix.

use linear_recurrence::berlekamp_massey;
use modint::StaticModInt;

struct SplitMix(u64);

impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// The determinant of the `n x n` matrix with the given nonzero entries `(row, col, value)`
/// (distinct positions, values below `M`) in `O(n (n + nnz))`: with a random diagonal `D`, the
/// minimal polynomial of `A D` found by Berlekamp–Massey on `u^T (A D)^t v` gives `det(A D)`.
/// Randomized (fixed seeds): a nonsingular matrix is reported as singular with probability about
/// `n^2 / M` per attempt, three attempts are made.
pub fn sparse_determinant<const M: u32>(n: usize, entries: &[(usize, usize, u32)]) -> u32 {
    type Z<const M: u32> = StaticModInt<M>;
    if n == 0 {
        return 1;
    }
    assert!(entries.iter().all(|&(i, j, c)| i < n && j < n && c < M), "entry out of range");
    let mut rng = SplitMix(0x1234_5678_9ABC_DEF0);
    for _ in 0..3 {
        let rand = |rng: &mut SplitMix| Z::<M>::new(1 + rng.next() % (M as u64 - 1));
        let d: Vec<Z<M>> = (0..n).map(|_| rand(&mut rng)).collect();
        let u: Vec<Z<M>> = (0..n).map(|_| rand(&mut rng)).collect();
        let mut w: Vec<Z<M>> = (0..n).map(|_| rand(&mut rng)).collect();
        // B = A D has entries c d_j at (i, j)
        let b: Vec<(usize, usize, Z<M>)> = entries.iter().map(|&(i, j, c)| (i, j, Z::<M>::new(c as u64) * d[j])).collect();
        let mut seq = Vec::with_capacity(2 * n);
        for _ in 0..2 * n {
            seq.push(u.iter().zip(&w).fold(Z::<M>::zero(), |s, (&x, &y)| s + x * y));
            let mut next = vec![Z::<M>::zero(); n];
            for &(i, j, c) in &b {
                next[i] += c * w[j];
            }
            w = next;
        }
        let rec = berlekamp_massey(&seq);
        if rec.len() < n {
            continue;
        }
        // p(x) = x^n - c_1 x^{n-1} - ... - c_n and det(-B) = p(0)
        let mut det_b = rec[n - 1];
        if n.is_multiple_of(2) {
            det_b = -det_b;
        }
        let det_d = d.iter().fold(Z::<M>::one(), |s, &x| s * x);
        return (det_b / det_d).val();
    }
    0
}
