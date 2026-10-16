//! Gaussian integers: gcd and the representations as a sum of two squares.

use super::{factorize, mod_sqrt, NONE};

type G = (i128, i128);

fn mul(a: G, b: G) -> G {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}

/// `round(a / b)` componentwise (`b != 0`), the quotient of the Euclidean division.
fn div_round(a: G, b: G) -> G {
    let n = b.0 * b.0 + b.1 * b.1;
    let (re, im) = mul(a, (b.0, -b.1));
    let r = |v: i128| (2 * v + n).div_euclid(2 * n);
    (r(re), r(im))
}

fn gcd_inner(mut a: G, mut b: G) -> G {
    while b != (0, 0) {
        let q = div_round(a, b);
        let qb = mul(q, b);
        (a, b) = (b, (a.0 - qb.0, a.1 - qb.1));
    }
    a
}

/// A greatest common divisor of `a0 + a1 i` and `b0 + b1 i` in `Z[i]` (unique up to a unit),
/// by the Euclidean algorithm. `gcd(0, 0) = 0`.
pub fn gaussian_gcd(a: (i64, i64), b: (i64, i64)) -> (i64, i64) {
    let g = gcd_inner((a.0 as i128, a.1 as i128), (b.0 as i128, b.1 as i128));
    (g.0 as i64, g.1 as i64)
}

/// All ordered pairs `(a, b)` of non-negative integers with `a^2 + b^2 = n`, sorted.
/// Expected time `O(n^{1/4})` (factorisation) plus the output size.
pub fn two_squares(n: u64) -> Vec<(u64, u64)> {
    if n == 0 {
        return vec![(0, 0)];
    }
    // all Gaussian integers of norm n, up to units, as products over the prime factors
    let mut elems: Vec<G> = vec![(1, 0)];
    for (p, e) in factorize(n) {
        let e = e as usize;
        if p == 2 {
            let mut f = (1i128, 0i128);
            for _ in 0..e {
                f = mul(f, (1, 1));
            }
            elems = elems.iter().map(|&z| mul(z, f)).collect();
        } else if p % 4 == 3 {
            if e % 2 == 1 {
                return Vec::new();
            }
            let mut f = 1i128;
            for _ in 0..e / 2 {
                f *= p as i128;
            }
            elems = elems.iter().map(|&z| (z.0 * f, z.1 * f)).collect();
        } else {
            let r = mod_sqrt(p - 1, p);
            assert!(r != NONE);
            let pi = gcd_inner((p as i128, 0), (r as i128, 1));
            let conj = (pi.0, -pi.1);
            // pi^j conj^(e - j) for j = 0..=e
            let mut pows_pi = vec![(1i128, 0i128)];
            let mut pows_c = vec![(1i128, 0i128)];
            for _ in 0..e {
                pows_pi.push(mul(*pows_pi.last().unwrap(), pi));
                pows_c.push(mul(*pows_c.last().unwrap(), conj));
            }
            let mut next = Vec::with_capacity(elems.len() * (e + 1));
            for &z in &elems {
                for j in 0..=e {
                    next.push(mul(z, mul(pows_pi[j], pows_c[e - j])));
                }
            }
            elems = next;
        }
    }
    let mut out = Vec::new();
    for z in elems {
        let mut w = z;
        for _ in 0..4 {
            if w.0 >= 0 && w.1 >= 0 {
                out.push((w.0 as u64, w.1 as u64));
            }
            w = (-w.1, w.0); // times i
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}
