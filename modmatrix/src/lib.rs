//! Dense matrices over `Z/MZ` (prime `M`).
//!
//! * [`determinant`]: Gaussian elimination on a `Vec<Vec<StaticModInt<M>>>`, `O(n^3)`.
//! * [`Matrix`]: a flat row-major matrix of residues with product (cache-blocked, AVX2 when
//!   available), power, determinant, rank, inverse, adjugate (also for singular matrices),
//!   linear systems, characteristic polynomial (Hessenberg form) and Pfaffian.
//! * [`hafnian`]: the hafnian of a symmetric matrix with zero diagonal in `O(2^{n/2} n^2 d^2)`.
//! * [`sparse_determinant`]: Wiedemann's algorithm for a sparse matrix.
//! * [`determinant_mod`]: determinant modulo an arbitrary integer (Euclidean elimination).
//!
//! # Example
//!
//! ```
//! use modint::ModInt998244353 as Mint;
//! use modmatrix::determinant;
//!
//! let a = vec![vec![Mint::new(1), Mint::new(2)], vec![Mint::new(3), Mint::new(4)]];
//! assert_eq!(determinant(a), -Mint::new(2));
//! ```

mod arbitrary;
mod hafnian;
mod matrix;
mod sparse;

pub use arbitrary::determinant_mod;
pub use hafnian::hafnian;
pub use matrix::Matrix;
pub use sparse::sparse_determinant;

use modint::StaticModInt;

/// `base^e mod M`.
pub(crate) fn pow_u32<const M: u32>(base: u32, mut e: u64) -> u32 {
    let (mut r, mut b) = (1u64, base as u64 % M as u64);
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % M as u64;
        }
        b = b * b % M as u64;
        e >>= 1;
    }
    r as u32
}

/// Determinant of a square matrix (consumed, as it is eliminated in place). The determinant of
/// the `0 x 0` matrix is 1.
///
/// # Panics
/// Panics if the matrix is not square.
pub fn determinant<const M: u32>(mut a: Vec<Vec<StaticModInt<M>>>) -> StaticModInt<M> {
    let n = a.len();
    assert!(a.iter().all(|row| row.len() == n), "determinant needs a square matrix");
    let mut det = StaticModInt::<M>::one();
    for col in 0..n {
        let Some(pivot) = (col..n).find(|&r| a[r][col].val() != 0) else {
            return StaticModInt::zero();
        };
        if pivot != col {
            a.swap(pivot, col);
            det = -det;
        }
        let p = a[col][col];
        det *= p;
        let inv = p.inv();
        let (head, tail) = a.split_at_mut(col + 1);
        let prow = &head[col];
        for row in tail.iter_mut() {
            if row[col].val() == 0 {
                continue;
            }
            let f = row[col] * inv;
            for j in col..n {
                row[j] -= f * prow[j];
            }
        }
    }
    det
}
