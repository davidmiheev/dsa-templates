//! Products of polynomials in `K` variables stored as dense arrays in mixed radix: the entry at
//! `i = i_1 + i_2 n_1 + i_3 n_1 n_2 + ...` is the coefficient of `x_1^{i_1} x_2^{i_2} ...`.
//!
//! * [`multivariate_convolution`] (modulo 998244353): the product modulo `x_j^{n_j}` for every
//!   `j`. Kronecker substitution with the radices of the *input* layout would let digit carries
//!   pollute the result; the number of carries is detected by the grading
//!   `chi(x) = sum_{j < K-1} floor(x / (n_1 ... n_j))`, which changes by exactly the number of
//!   "carried" prefix sums (`0..K-1`), so splitting the arrays into `K` classes by `chi mod K`,
//!   convolving class pairs cyclically in the class index, and reading each position from its
//!   own class keeps exactly the carry-free terms. `K` transforms for each input, `K^2` pointwise
//!   products.
//! * [`multivariate_convolution_cyclic`] (modulo a prime `p` with every `n_j | p - 1`): the product
//!   modulo `x_j^{n_j} - 1`, by a length-`n_j` DFT along every axis (primitive roots of unity
//!   exist because `n_j | p - 1`); axes longer than 64 use Bluestein's chirp-z transform with an
//!   arbitrary-modulus convolution.
//!
//! # Example
//!
//! ```
//! use multivariate_convolution::{multivariate_convolution, multivariate_convolution_cyclic};
//!
//! // (1 + x)(1 + x) mod x^2 = 1 + 2x
//! assert_eq!(multivariate_convolution(&[2], &[1, 1], &[1, 1]), vec![1, 2]);
//! // (1 + x)(1 + x) mod (x^2 - 1) = 2 + 2x, over F_3
//! assert_eq!(multivariate_convolution_cyclic(3, &[2], &[1, 1], &[1, 1]), vec![2, 2]);
//! ```

use convolution::{convolution, convolution_arbitrary_mod, intt_unscaled, ntt, P1};
use number_theory::{mod_inverse, pow_mod, primitive_root};

const MOD: u64 = P1 as u64;

/// Product of `f` and `g` modulo `x_j^{shape[j]}` for all `j`, coefficients modulo 998244353.
/// An empty `shape` means constants (length-1 arrays).
pub fn multivariate_convolution(shape: &[usize], f: &[u32], g: &[u32]) -> Vec<u32> {
    let k = shape.len();
    let total: usize = shape.iter().product();
    assert!(f.len() == total && g.len() == total, "arrays must have length {total}");
    if k == 0 {
        return vec![(f[0] as u64 * g[0] as u64 % MOD) as u32];
    }
    if k == 1 {
        let mut c = convolution(f, g);
        c.truncate(total);
        return c;
    }
    // chi(x) = sum over prefix products B_1 .. B_{K-1} of floor(x / B)
    let mut bases = Vec::with_capacity(k - 1);
    let mut b = 1usize;
    for &s in &shape[..k - 1] {
        b *= s;
        bases.push(b);
    }
    let chi: Vec<usize> = (0..total).map(|x| bases.iter().map(|&b| x / b).sum::<usize>() % k).collect();
    let size = (2 * total).next_power_of_two();
    let split = |v: &[u32]| -> Vec<Vec<u32>> {
        let mut parts = vec![vec![0u32; size]; k];
        for (x, &c) in v.iter().enumerate() {
            parts[chi[x]][x] = c;
        }
        for p in parts.iter_mut() {
            ntt::<P1>(p);
        }
        parts
    };
    let (fa, fb) = (split(f), split(g));
    let mut out = vec![vec![0u32; size]; k];
    for (l, x) in fa.iter().enumerate() {
        for (r, y) in fb.iter().enumerate() {
            let lr = if l + r >= k { l + r - k } else { l + r };
            for ((o, &u), &v) in out[lr].iter_mut().zip(x).zip(y) {
                let s = *o as u64 + u as u64 * v as u64 % MOD;
                *o = (if s >= MOD { s - MOD } else { s }) as u32;
            }
        }
    }
    for o in out.iter_mut() {
        intt_unscaled::<P1>(o);
    }
    let inv = pow_mod(size as u64, MOD - 2, MOD);
    (0..total).map(|x| (out[chi[x]][x] as u64 * inv % MOD) as u32).collect()
}

/// DFT of length `n` modulo the prime `p` with the primitive `n`-th root of unity `w`:
/// `X[k] = sum_j x[j] w^{jk}`.
fn dft_line(x: &[u64], w: u64, p: u64) -> Vec<u64> {
    let n = x.len();
    if n <= 64 {
        let mut pw = vec![1u64; n];
        for i in 1..n {
            pw[i] = pw[i - 1] * w % p;
        }
        return (0..n)
            .map(|k| (0..n).fold(0u64, |s, j| (s + x[j] * pw[j * k % n]) % p))
            .collect();
    }
    // Bluestein: jk = C(j + k, 2) - C(j, 2) - C(k, 2)
    let winv = mod_inverse(w, p);
    let tri = |t: usize| -> u64 { (t as u64 * (t as u64).wrapping_sub(1) / 2) % n as u64 }; // C(t,2) mod n
    let mut wp = vec![1u64; n];
    let mut wi = vec![1u64; n];
    for i in 1..n {
        wp[i] = wp[i - 1] * w % p;
        wi[i] = wi[i - 1] * winv % p;
    }
    let u: Vec<u32> = (0..n).rev().map(|j| (x[j] * wi[tri(j) as usize] % p) as u32).collect();
    let v: Vec<u32> = (0..2 * n - 1).map(|t| wp[tri(t) as usize] as u32).collect();
    let conv = convolution_arbitrary_mod(&u, &v, p as u32);
    (0..n).map(|k| conv[n - 1 + k] as u64 * wi[tri(k) as usize] % p).collect()
}

/// Product of `f` and `g` modulo `x_j^{shape[j]} - 1` for all `j`, coefficients modulo the prime
/// `p` (below `2^31`). Requires `shape[j] | p - 1`.
pub fn multivariate_convolution_cyclic(p: u64, shape: &[usize], f: &[u32], g: &[u32]) -> Vec<u32> {
    let total: usize = shape.iter().product();
    assert!(f.len() == total && g.len() == total, "arrays must have length {total}");
    if shape.is_empty() {
        return vec![(f[0] as u64 * g[0] as u64 % p) as u32];
    }
    let root = primitive_root(p);
    let mut a: Vec<u64> = f.iter().map(|&x| x as u64).collect();
    let mut b: Vec<u64> = g.iter().map(|&x| x as u64).collect();
    let transform = |data: &mut Vec<u64>, inverse: bool| {
        let mut stride = 1usize;
        for &n in shape {
            assert!((p - 1) % n as u64 == 0, "axis length {n} does not divide p - 1 = {}", p - 1);
            let mut w = pow_mod(root, (p - 1) / n as u64, p);
            if inverse {
                w = mod_inverse(w, p);
            }
            let block = stride * n;
            for start in (0..total).step_by(block) {
                for off in 0..stride {
                    let line: Vec<u64> = (0..n).map(|j| data[start + off + j * stride]).collect();
                    let out = dft_line(&line, w, p);
                    for (j, &v) in out.iter().enumerate() {
                        data[start + off + j * stride] = v;
                    }
                }
            }
            stride = block;
        }
    };
    transform(&mut a, false);
    transform(&mut b, false);
    for (x, &y) in a.iter_mut().zip(&b) {
        *x = *x * y % p;
    }
    transform(&mut a, true);
    let scale = mod_inverse(total as u64 % p, p);
    a.iter().map(|&x| (x * scale % p) as u32).collect()
}
