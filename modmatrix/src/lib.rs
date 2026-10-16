//! Dense matrices over `Z/MZ` (prime `M`) as `Vec<Vec<StaticModInt<M>>>`.
//!
//! [`determinant`] is Gaussian elimination in `O(n^3)`.
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

use modint::StaticModInt;

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
