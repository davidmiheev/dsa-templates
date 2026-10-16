//! Series with a handful of non-zero coefficients (`K` terms): every function runs in
//! `O(N K)` by the linear recurrence that the defining differential equation gives, instead of
//! the `O(N log N)` Newton iterations. Input is the length `n` and the terms `(index, value)`
//! in increasing index order with non-zero values.

use modint::StaticModInt;
use number_theory::{mod_sqrt, NONE};

type Mint<const M: u32> = StaticModInt<M>;

fn inverses<const M: u32>(n: usize) -> Vec<Mint<M>> {
    crate::fps::inverses::<M>(n + 1)
}

/// `1 / f` modulo `x^n`; the term at index 0 must exist.
pub fn inv_sparse<const M: u32>(n: usize, terms: &[(usize, Mint<M>)]) -> Vec<Mint<M>> {
    assert!(!terms.is_empty() && terms[0].0 == 0, "the constant term is required");
    let a0_inv = terms[0].1.inv();
    let mut g = vec![Mint::<M>::zero(); n];
    for i in 0..n {
        let mut s = if i == 0 { Mint::one() } else { Mint::zero() };
        for &(k, v) in &terms[1..] {
            if k > i {
                break;
            }
            s -= v * g[i - k];
        }
        g[i] = s * a0_inv;
    }
    g
}

/// `log f` modulo `x^n` where the constant term of `f` is 1 (given as the term at index 0).
pub fn log_sparse<const M: u32>(n: usize, terms: &[(usize, Mint<M>)]) -> Vec<Mint<M>> {
    assert!(!terms.is_empty() && terms[0] == (0, Mint::one()), "log needs constant term 1");
    let inv = inverses::<M>(n);
    // n g_n = n f_n - sum_{k=1}^{n-1} f_k (n - k) g_{n-k}
    let mut g = vec![Mint::<M>::zero(); n];
    for i in 1..n {
        let mut s = Mint::<M>::zero();
        for &(k, v) in &terms[1..] {
            if k > i {
                break;
            }
            if k == i {
                s += v * Mint::from(i);
            } else {
                s -= v * Mint::from(i - k) * g[i - k];
            }
        }
        g[i] = s * inv[i];
    }
    g
}

/// `exp f` modulo `x^n` (no constant term).
pub fn exp_sparse<const M: u32>(n: usize, terms: &[(usize, Mint<M>)]) -> Vec<Mint<M>> {
    assert!(terms.iter().all(|&(i, _)| i > 0), "exp needs constant term 0");
    let inv = inverses::<M>(n);
    // n g_n = sum_k k f_k g_{n-k}
    let mut g = vec![Mint::<M>::zero(); n];
    if n > 0 {
        g[0] = Mint::one();
    }
    for i in 1..n {
        let mut s = Mint::<M>::zero();
        for &(k, v) in terms {
            if k > i {
                break;
            }
            s += v * Mint::from(k) * g[i - k];
        }
        g[i] = s * inv[i];
    }
    g
}

/// `h^e` for `h` with constant term 1 by `h g' = e h' g`; `e` is an element of the field.
fn pow_one_sparse<const M: u32>(n: usize, terms: &[(usize, Mint<M>)], e: Mint<M>) -> Vec<Mint<M>> {
    let inv = inverses::<M>(n);
    // n g_n = sum_{j>=1} (e j - (n - j)) h_j g_{n-j}
    let mut g = vec![Mint::<M>::zero(); n];
    if n > 0 {
        g[0] = Mint::one();
    }
    for i in 1..n {
        let mut s = Mint::<M>::zero();
        for &(j, v) in terms {
            if j > i {
                break;
            }
            s += (e * Mint::from(j) - Mint::from(i - j)) * v * g[i - j];
        }
        g[i] = s * inv[i];
    }
    g
}

/// `f^e` modulo `x^n` for any exponent `e >= 0`.
pub fn pow_sparse<const M: u32>(n: usize, terms: &[(usize, Mint<M>)], e: u64) -> Vec<Mint<M>> {
    let mut result = vec![Mint::<M>::zero(); n];
    if n == 0 {
        return result;
    }
    if e == 0 {
        result[0] = Mint::one();
        return result;
    }
    let Some(&(k, c)) = terms.first() else { return result };
    if (k as u128) * (e as u128) >= n as u128 {
        return result;
    }
    let shift = k * e as usize;
    let c_inv = c.inv();
    let h: Vec<(usize, Mint<M>)> = terms.iter().map(|&(i, v)| (i - k, v * c_inv)).collect();
    let len = n - shift;
    let g = pow_one_sparse(len, &h, Mint::from(e % M as u64));
    let scale = c.pow(e % (M as u64 - 1));
    for i in 0..len {
        result[shift + i] = g[i] * scale;
    }
    result
}

/// A square root of `f` modulo `x^n`, or an **empty vector** if none exists (`n >= 1`).
pub fn sqrt_sparse<const M: u32>(n: usize, terms: &[(usize, Mint<M>)]) -> Vec<Mint<M>> {
    assert!(n >= 1, "need at least one term");
    let Some(&(k, c)) = terms.first() else { return vec![Mint::zero(); n] };
    if k >= n {
        return vec![Mint::zero(); n];
    }
    if k % 2 == 1 {
        return Vec::new();
    }
    let r = mod_sqrt(c.val() as u64, M as u64);
    if r == NONE {
        return Vec::new();
    }
    let c_inv = c.inv();
    let h: Vec<(usize, Mint<M>)> = terms.iter().map(|&(i, v)| (i - k, v * c_inv)).collect();
    let len = n - k / 2;
    let half = Mint::<M>::new(2).inv();
    let g = pow_one_sparse(len, &h, half);
    let mut result = vec![Mint::zero(); n];
    for i in 0..len {
        result[k / 2 + i] = g[i] * Mint::new(r);
    }
    result
}
