//! Berlekamp–Massey over `Z/MZ` (prime `M`).
//!
//! [`berlekamp_massey`] returns the shortest `c_1, ..., c_L` with
//! `s_n = c_1 s_{n-1} + ... + c_L s_{n-L}` for all `n >= L`, in `O(n L)`.
//!
//! # Example
//!
//! ```
//! use linear_recurrence::berlekamp_massey;
//! use modint::ModInt998244353 as Mint;
//!
//! // Fibonacci numbers
//! let s: Vec<Mint> = [0u64, 1, 1, 2, 3, 5, 8, 13, 21, 34].iter().map(|&x| Mint::new(x)).collect();
//! assert_eq!(berlekamp_massey(&s), vec![Mint::new(1), Mint::new(1)]);
//! ```

use modint::StaticModInt;

/// The shortest linear recurrence `s_n = sum_{i=1}^{L} c_i s_{n-i}` (`n >= L`) satisfied by `s`,
/// returned as `[c_1, ..., c_L]`; empty for the zero sequence.
pub fn berlekamp_massey<const M: u32>(s: &[StaticModInt<M>]) -> Vec<StaticModInt<M>> {
    type Z<const M: u32> = StaticModInt<M>;
    let zero = Z::<M>::zero();
    // c: connection polynomial 1 + c[1] x + ... ; b: its value before the last length change
    let mut c: Vec<Z<M>> = vec![Z::one()];
    let mut b: Vec<Z<M>> = vec![Z::one()];
    let (mut l, mut m) = (0usize, 1usize);
    let mut bb = Z::<M>::one();
    for i in 0..s.len() {
        let mut d = s[i];
        for j in 1..=l.min(c.len() - 1) {
            d += c[j] * s[i - j];
        }
        if d == zero {
            m += 1;
            continue;
        }
        let coef = d / bb;
        let prev = if 2 * l <= i { Some(c.clone()) } else { None };
        if c.len() < b.len() + m {
            c.resize(b.len() + m, zero);
        }
        for (j, &bj) in b.iter().enumerate() {
            c[j + m] -= coef * bj;
        }
        if let Some(prev) = prev {
            l = i + 1 - l;
            b = prev;
            bb = d;
            m = 1;
        } else {
            m += 1;
        }
    }
    c.resize(l + 1, zero);
    c[1..].iter().map(|&x| -x).collect()
}
