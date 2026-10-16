//! Composition of power series in `O(n log^2 n)` (Kinoshita–Li) and the compositional inverse.
//!
//! **Power projection.** For `g` with `g(0) = 0` the map `w -> p`, `p_k = sum_m w_m [x^m] g^k`, is
//! `p(y) = [x^{n-1}] P_0(x) / (1 - y g(x))` with `P_0 = reverse(w)`. Bostan–Mori halves the target
//! index at every step: multiply numerator and denominator by `Q(-x, y)`, keep the parity class of the
//! target index. Over `log n` steps the `x`-size halves while the `y`-degree doubles, so every step
//! costs one bivariate product of total size `O(n)`.
//!
//! **Composition** `f(g)` is the *transpose* of power projection (`c_m = sum_k a_k [x^m] g^k` is the
//! same matrix applied the other way). Transposing the algorithm step by step — selection of a parity
//! class becomes an embedding, a product with a fixed polynomial becomes a correlation (a "middle
//! product"), steps in reverse order — gives an algorithm with the same complexity that never forms
//! the powers of `g`. The denominators `Q_t` do not depend on the input vector, so they are computed
//! forwards once and stored.
//!
//! **Compositional inverse** by Newton's iteration `g <- g - (f(g) - x) / f'(g)`, where one
//! composition per doubling suffices because `f'(g) = (f(g))' / g'`.

use crate::fps::{derivative, inv, mul, mul_trunc, resized};
use modint::StaticModInt;

type Mint<const M: u32> = StaticModInt<M>;

/// The denominators `Q_0, ..., Q_T` of the Bostan–Mori recursion for `1 - y g(x)`: level `t` has
/// `m_t + 1` rows (`m_t = (n - 1) >> t`) and `2^t + 1` columns, row-major.
fn denominators<const M: u32>(g: &[Mint<M>], n: usize) -> Vec<Vec<Mint<M>>> {
    let levels = usize::BITS as usize - (n - 1).leading_zeros() as usize; // T = bit length of n - 1
    let mut q0 = vec![Mint::<M>::zero(); n * 2];
    for i in 0..n {
        q0[2 * i + 1] = -if i < g.len() { g[i] } else { Mint::zero() };
    }
    q0[0] += Mint::one();
    let mut qs = vec![q0];
    for t in 0..levels {
        let m_t = (n - 1) >> t;
        let m_next = m_t >> 1;
        let cols = (1usize << t) + 1;
        let stride = (1usize << (t + 1)) + 1; // y-degree of the product is 2^(t+1)
        let q = &qs[t];
        let spread = |negate_odd: bool| -> Vec<Mint<M>> {
            let mut v = vec![Mint::<M>::zero(); m_t * stride + cols];
            for i in 0..=m_t {
                for y in 0..cols {
                    let x = q[i * cols + y];
                    v[i * stride + y] = if negate_odd && i % 2 == 1 { -x } else { x };
                }
            }
            v
        };
        let c = mul(&spread(false), &spread(true));
        let mut next = vec![Mint::<M>::zero(); (m_next + 1) * stride];
        for j in 0..=m_next {
            for y in 0..stride {
                let pos = 2 * j * stride + y;
                if pos < c.len() {
                    next[j * stride + y] = c[pos];
                }
            }
        }
        qs.push(next);
    }
    qs
}

/// The coefficients `[x^0 .. x^{n-1}]` of `f(g(x))` where `g(0) = 0`. `f` and `g` are truncated or
/// zero-padded to `n` terms.
///
/// # Panics
/// Panics if `g` has a non-zero constant term.
pub fn compose<const M: u32>(f: &[Mint<M>], g: &[Mint<M>], n: usize) -> Vec<Mint<M>> {
    assert!(g.is_empty() || g[0].val() == 0, "composition needs g(0) = 0");
    if n == 0 {
        return Vec::new();
    }
    let qs = denominators(g, n);
    let levels = qs.len() - 1;
    // transposed pass; z has (m_t + 1) rows and 2^t columns at level t
    let mut z = vec![Mint::<M>::zero(); 1 << levels];
    for (k, slot) in z.iter_mut().enumerate().take(n.min(f.len())) {
        *slot = f[k];
    }
    for t in (0..levels).rev() {
        let m_t = (n - 1) >> t;
        let b_t = m_t & 1;
        let m_next = m_t >> 1;
        let w = 1usize << (t + 1); // columns of z (level t + 1) and flattening stride
        let cols_q = (1usize << t) + 1;
        // R: rows 2j + b_t carry row j of z
        let mut r = vec![Mint::<M>::zero(); (m_t + 1) * w];
        for j in 0..=m_next {
            let row = 2 * j + b_t;
            r[row * w..row * w + w].copy_from_slice(&z[j * w..(j + 1) * w]);
        }
        // Q~ rows i = 0..m_t, flattened with stride w, then reversed for a correlation
        let q = &qs[t];
        let len_q = m_t * w + cols_q;
        let mut qt = vec![Mint::<M>::zero(); len_q];
        for i in 0..=m_t {
            for y in 0..cols_q {
                let x = q[i * cols_q + y];
                qt[i * w + y] = if i % 2 == 1 { -x } else { x };
            }
        }
        qt.reverse();
        let conv = mul(&r, &qt);
        let cols_new = 1usize << t;
        let mut next = vec![Mint::<M>::zero(); (m_t + 1) * cols_new];
        for i in 0..=m_t {
            for u in 0..cols_new {
                let pos = i * w + u + len_q - 1;
                if pos < conv.len() {
                    next[i * cols_new + u] = conv[pos];
                }
            }
        }
        z = next;
    }
    // level 0: n rows, one column; c_m = z[n - 1 - m]
    (0..n).map(|m| z[n - 1 - m]).collect()
}

/// The compositional inverse of `f` (`f(0) = 0`, `f'(0) != 0`) modulo `x^n`: the series `g` with
/// `f(g(x)) = x`.
pub fn compositional_inverse<const M: u32>(f: &[Mint<M>], n: usize) -> Vec<Mint<M>> {
    assert!(f.len() >= 2 || n <= 1, "need the linear term");
    assert!(f.is_empty() || f[0].val() == 0, "f(0) must be 0");
    if n <= 1 {
        return vec![Mint::zero(); n];
    }
    assert!(f[1].val() != 0, "f'(0) must be invertible");
    let mut g = vec![Mint::<M>::zero(), f[1].inv()];
    let mut m = 2;
    while m < n {
        let target = (2 * m).min(n.next_power_of_two().max(2 * m));
        let c = compose(f, &resized(&g, target), target);
        // E = f(g) - x, which is divisible by x^m
        let mut e = c.clone();
        e[1] -= Mint::one();
        // f'(g) = (f(g))' / g' modulo x^m
        let gp = derivative(&resized(&g, target));
        let fp = mul_trunc(&derivative(&c), &inv(&gp, m), m);
        let step = mul_trunc(&e, &inv(&fp, m), target);
        let mut next = resized(&g, target);
        for i in 0..target {
            next[i] -= step[i];
        }
        g = next;
        m = target;
    }
    g.truncate(n);
    g
}
