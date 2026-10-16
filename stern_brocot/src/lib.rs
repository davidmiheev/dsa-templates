//! The Stern–Brocot tree: every positive rational `a/b` (in lowest terms) is a node, `1/1` is the
//! root, and the path to a node is a string of `L` / `R` moves, stored run-length encoded as
//! `(is_right, count)` pairs with alternating directions. All operations are `O(log)` in the
//! numerators and denominators (the paths themselves can be exponentially long).
//!
//! * [`encode_path`] / [`decode_path`]: the run-length path of a fraction and back.
//! * [`lca`], [`ancestor`] (the ancestor at a given depth; `(0, 0)` if the node is not that deep)
//!   and [`descendant_range`] (the open interval holding the subtree).
//! * [`best_approximations`]: the closest fractions to `x/y` from below and above with
//!   numerator and denominator at most `n`.
//!
//! A fraction is a pair `(numerator, denominator)`; `(1, 0)` stands for `+infinity` and `(0, 1)`
//! for `0` where they appear as interval ends.
//!
//! # Example
//!
//! ```
//! use stern_brocot::*;
//!
//! // 1/1 -R-> 2/1 -L-> 3/2 -R-> 5/3
//! assert_eq!(encode_path(5, 3), vec![(true, 1), (false, 1), (true, 1)]);
//! assert_eq!(decode_path(&[(true, 1), (false, 1), (true, 1)]), (5, 3));
//! assert_eq!(lca((5, 3), (3, 2)), (3, 2));
//! assert_eq!(ancestor(1, (5, 3)), (2, 1));
//! assert_eq!(best_approximations(30, 31_415, 10_000), ((25, 8), (22, 7)));
//! ```

/// A run of `count` moves, to the right child if `.0`, otherwise to the left.
pub type Run = (bool, u64);

/// The run-length encoded path from `1/1` to `a/b` (`gcd(a, b) = 1`, `a, b >= 1`).
pub fn encode_path(mut a: u64, mut b: u64) -> Vec<Run> {
    assert!(a >= 1 && b >= 1, "not a positive fraction");
    let mut runs = Vec::new();
    while !(a == 1 && b == 1) {
        if a > b {
            let n = (a - 1) / b;
            runs.push((true, n));
            a -= n * b;
        } else {
            let n = (b - 1) / a;
            assert!(n >= 1, "the fraction is not in lowest terms");
            runs.push((false, n));
            b -= n * a;
        }
    }
    runs
}

/// The two ends `(lo, hi)` of the open interval that holds the subtree of the node reached by
/// `runs`; the node itself is the mediant of the ends.
fn bounds(runs: &[Run]) -> ((u64, u64), (u64, u64)) {
    let (mut lo, mut hi) = ((0u64, 1u64), (1u64, 0u64));
    for &(right, n) in runs {
        if right {
            lo = (lo.0 + n * hi.0, lo.1 + n * hi.1);
        } else {
            hi = (hi.0 + n * lo.0, hi.1 + n * lo.1);
        }
    }
    (lo, hi)
}

/// The fraction at the end of the path.
pub fn decode_path(runs: &[Run]) -> (u64, u64) {
    let (lo, hi) = bounds(runs);
    (lo.0 + hi.0, lo.1 + hi.1)
}

/// The lowest common ancestor of two nodes.
pub fn lca(u: (u64, u64), v: (u64, u64)) -> (u64, u64) {
    let (pu, pv) = (encode_path(u.0, u.1), encode_path(v.0, v.1));
    let mut common: Vec<Run> = Vec::new();
    for (&(cu, nu), &(cv, nv)) in pu.iter().zip(&pv) {
        if cu != cv {
            break;
        }
        common.push((cu, nu.min(nv)));
        if nu != nv {
            break;
        }
    }
    decode_path(&common)
}

/// The ancestor of `u` at depth `k` (the root has depth 0; `u` is its own ancestor at its own
/// depth), or `(0, 0)` if `u` is shallower than `k`.
pub fn ancestor(k: u64, u: (u64, u64)) -> (u64, u64) {
    let mut left = k;
    let mut kept: Vec<Run> = Vec::new();
    for (right, n) in encode_path(u.0, u.1) {
        let take = n.min(left);
        if take > 0 {
            kept.push((right, take));
        }
        left -= take;
        if left == 0 {
            break;
        }
    }
    if left > 0 {
        return (0, 0);
    }
    decode_path(&kept)
}

/// The infimum `f/g` and supremum `h/k` of the descendants of `u`, as `((f, g), (h, k))`.
pub fn descendant_range(u: (u64, u64)) -> ((u64, u64), (u64, u64)) {
    bounds(&encode_path(u.0, u.1))
}

/// For `n, x, y >= 1`: the largest fraction `p/q <= x/y` and the smallest `p/q >= x/y` with
/// `1 <= p, q <= n`, as `((p, q), (r, s))`. `(0, 1)` / `(1, 0)` stand for "none" below / above.
pub fn best_approximations(n: u64, x: u64, y: u64) -> ((u64, u64), (u64, u64)) {
    assert!(n >= 1 && x >= 1 && y >= 1, "arguments must be positive");
    let g = gcd(x, y);
    let (x, y) = (x / g, y / g);
    if x <= n && y <= n {
        return ((x, y), (x, y));
    }
    let (n, x, y) = (n as i128, x as i128, y as i128);
    let (mut lo, mut hi) = ((0i128, 1i128), (1i128, 0i128));
    loop {
        // steps to the right: lo += j hi while lo + j hi < x/y and the coordinates stay <= n
        let num = x * lo.1 - lo.0 * y; // > 0
        let den = hi.0 * y - x * hi.1; // > 0
        let mut j = num / den;
        if hi.0 > 0 {
            j = j.min((n - lo.0) / hi.0);
        }
        if hi.1 > 0 {
            j = j.min((n - lo.1) / hi.1);
        }
        let right = j > 0;
        if right {
            lo = (lo.0 + j * hi.0, lo.1 + j * hi.1);
        }
        // steps to the left: hi += j lo while hi + j lo > x/y and the coordinates stay <= n
        let num = hi.0 * y - x * hi.1;
        let den = x * lo.1 - lo.0 * y;
        let mut j = num / den;
        if lo.0 > 0 {
            j = j.min((n - hi.0) / lo.0);
        }
        if lo.1 > 0 {
            j = j.min((n - hi.1) / lo.1);
        }
        let left = j > 0;
        if left {
            hi = (hi.0 + j * lo.0, hi.1 + j * lo.1);
        }
        if !right && !left {
            return ((lo.0 as u64, lo.1 as u64), (hi.0 as u64, hi.1 as u64));
        }
    }
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}
