//! Exact graph colouring for small graphs (`n <= ~22`) by inclusion–exclusion over vertex subsets.
//!
//! Let `A_S(y) = sum over independent sets I inside S of y^|I|`. Then for every `k`
//!
//! ```text
//!   sum_S (-1)^(n-|S|) * A_S(y)^k  =  sum over k-tuples of independent sets with union V of y^(total size),
//! ```
//!
//! and the coefficient of `y^n` counts the tuples that are also disjoint, i.e. the proper
//! `k`-colourings `P(k)`. Putting `y = 1` instead counts covers, whose positivity decides the
//! chromatic number.
//!
//! * [`chromatic_number`] — `O(2^n * chi)`.
//! * [`chromatic_polynomial`] — `O(2^n * n^3)` worst case: with `A = 1 + u`, `P(k) = sum_i C(k, i) T_i`
//!   where `T_i = sum_S sign * [y^n] u_S^i`, and all `[y^n] u^i` come out of the recurrence
//!   `g_m(z) = z * sum_r u_r g_{m-r}(z)` for `g_m(z) = [y^m] 1 / (1 - z u(y))`.
//!
//! # Example
//!
//! ```
//! use graph_coloring::{chromatic_number, chromatic_polynomial};
//! use modint::ModInt998244353 as Mint;
//!
//! // A triangle needs 3 colours and has P(x) = x (x-1) (x-2) = x^3 - 3x^2 + 2x.
//! let tri = [(0, 1), (1, 2), (2, 0)];
//! assert_eq!(chromatic_number(3, &tri), 3);
//! let p = chromatic_polynomial::<998_244_353>(3, &tri);
//! assert_eq!(p, vec![Mint::new(0), Mint::new(2), -Mint::new(3), Mint::new(1)]);
//! ```

use modint::StaticModInt;

const P: u64 = 998_244_353;

fn adjacency(n: usize, edges: &[(usize, usize)]) -> (Vec<u32>, u32) {
    assert!(n <= 26, "subset DP is limited to 26 vertices, got {n}");
    let mut adj = vec![0u32; n];
    let mut loops = 0u32;
    for &(u, v) in edges {
        assert!(u < n && v < n, "edge ({u}, {v}) has an endpoint outside 0..{n}");
        if u == v {
            loops |= 1 << u;
        } else {
            adj[u] |= 1 << v;
            adj[v] |= 1 << u;
        }
    }
    (adj, loops)
}

/// Chromatic number of the simple graph `edges` on `0..n` (`n >= 1`, parallel edges allowed).
///
/// # Panics
/// Panics on a self-loop (no proper colouring exists), on `n == 0` or `n > 26`.
pub fn chromatic_number(n: usize, edges: &[(usize, usize)]) -> usize {
    assert!(n >= 1, "the empty graph has no vertices to colour");
    let (adj, loops) = adjacency(n, edges);
    assert!(loops == 0, "a graph with a self-loop has no proper colouring");
    // independent[S] = number of independent subsets of S
    let mut independent = vec![0u32; 1 << n];
    independent[0] = 1;
    for s in 1usize..1 << n {
        let v = s.trailing_zeros() as usize;
        let rest = s & (s - 1);
        independent[s] = independent[rest] + independent[rest & !(adj[v] as usize)];
    }
    // sign = (-1)^(n - |S|)
    let sign = |s: usize| (n - s.count_ones() as usize) & 1;
    let mut pw_p = vec![1u64; 1 << n];
    let mut pw_w = vec![1u64; 1 << n];
    for k in 1..=n {
        let (mut acc_p, mut acc_w) = (0u64, 0u64);
        for s in 0..1usize << n {
            let c = independent[s] as u64;
            pw_p[s] = pw_p[s] * c % P;
            pw_w[s] = pw_w[s].wrapping_mul(c);
            if sign(s) == 0 {
                acc_p += pw_p[s];
                acc_w = acc_w.wrapping_add(pw_w[s]);
            } else {
                acc_p += P - pw_p[s];
                acc_w = acc_w.wrapping_sub(pw_w[s]);
            }
            if acc_p >= 1 << 62 {
                acc_p %= P;
            }
        }
        // The cover count is a positive integer: it vanishes only if both residues do.
        if acc_p % P != 0 || acc_w != 0 {
            return k;
        }
    }
    n
}

/// Coefficients `p_0, ..., p_n` of the chromatic polynomial of the multigraph `edges` on
/// `0..n` modulo the prime `M` (a self-loop makes it identically zero).
///
/// # Panics
/// Panics if `n == 0` or `n > 22` (memory is `2^n * (n + 1)` words), or `M` is not the prime
/// 998 244 353 used by the accumulators.
pub fn chromatic_polynomial<const M: u32>(n: usize, edges: &[(usize, usize)]) -> Vec<StaticModInt<M>> {
    assert!(M as u64 == P, "the accumulators are specialised to 998244353");
    assert!(n >= 1 && n <= 22, "need 1..=22 vertices, got {n}");
    let (adj, loops) = adjacency(n, edges);
    let w = n + 1;
    // ranked[S * w + r] = number of independent r-subsets of S
    let mut ranked = vec![0u32; (1 << n) * w];
    ranked[0] = 1;
    for s in 1usize..1 << n {
        let v = s.trailing_zeros() as usize;
        let rest = s & (s - 1);
        let without = rest & !(adj[v] as usize);
        for r in 0..w {
            let mut c = ranked[rest * w + r];
            if loops >> v & 1 == 0 && r > 0 {
                c += ranked[without * w + r - 1];
            }
            ranked[s * w + r] = c;
        }
    }
    // T_i = sum_S sign(S) [y^n] u_S^i, accumulated as positive and negative parts.
    let (mut pos, mut neg) = (vec![0u64; w], vec![0u64; w]);
    let mut g = vec![0u64; w * w]; // g[m * w + i] = [y^m z^i] 1/(1 - z u(y))
    g[0] = 1;
    for s in 1usize..1 << n {
        let u = &ranked[s * w..s * w + w];
        let alpha = (0..w).rev().find(|&r| u[r] != 0).unwrap();
        if alpha == 0 {
            continue; // only loops: A_S = 1
        }
        for m in 1..=n {
            let mut acc = [0u64; 27];
            for r in 1..=m.min(alpha) {
                let ur = u[r] as u64;
                let src = &g[(m - r) * w..(m - r) * w + (m - r + 1)];
                for (i, &x) in src.iter().enumerate() {
                    acc[i + 1] += ur * x;
                }
            }
            for i in 0..=m {
                g[m * w + i] = acc[i] % P;
            }
        }
        let target = if (n - s.count_ones() as usize) & 1 == 0 { &mut pos } else { &mut neg };
        for i in 0..w {
            target[i] += g[n * w + i];
        }
    }
    type Mint<const M: u32> = StaticModInt<M>;
    let t: Vec<Mint<M>> = (0..w).map(|i| Mint::new(pos[i] % P) - Mint::new(neg[i] % P)).collect();
    // P(x) = sum_i T_i * C(x, i); expand the falling factorials.
    let mut poly = vec![Mint::<M>::zero(); w];
    let mut falling = vec![Mint::<M>::one()]; // prod_{j < i} (x - j)
    let mut fact_inv = Mint::<M>::one();
    for i in 0..w {
        if i > 0 {
            fact_inv /= Mint::from(i);
        }
        for (d, &c) in falling.iter().enumerate() {
            poly[d] += t[i] * fact_inv * c;
        }
        // falling *= (x - i)
        let mut next = vec![Mint::<M>::zero(); falling.len() + 1];
        for (d, &c) in falling.iter().enumerate() {
            next[d + 1] += c;
            next[d] -= c * Mint::from(i);
        }
        falling = next;
    }
    poly
}
