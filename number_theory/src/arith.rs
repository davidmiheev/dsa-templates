//! Floor sums, minimum of a linear function modulo `m`, Euler's totient and tetration.

use super::factorize;

/// `sum_{i=0}^{n-1} floor((a i + b) / m)`, computed modulo `2^128` (so exact whenever the result
/// fits), in `O(log m)`. Arguments must be below `2^63`.
pub fn floor_sum(n: u64, m: u64, a: u64, b: u64) -> u128 {
    assert!(m > 0, "m must be positive");
    assert!(n < 1 << 63 && m < 1 << 63 && a < 1 << 63 && b < 1 << 63, "arguments too large");
    let (mut n, mut m, mut a, mut b) = (n as u128, m as u128, a as u128, b as u128);
    let mut ans = 0u128;
    loop {
        if a >= m {
            ans = ans.wrapping_add((n * n.saturating_sub(1) / 2).wrapping_mul(a / m));
            a %= m;
        }
        if b >= m {
            ans = ans.wrapping_add(n.wrapping_mul(b / m));
            b %= m;
        }
        let y = a * n + b;
        if y < m {
            return ans;
        }
        n = y / m;
        b = y % m;
        std::mem::swap(&mut m, &mut a);
    }
}

/// `min over 0 <= x < n of (a x + b) mod m`, in `O(log m)`.
///
/// # Panics
/// Panics unless `n >= 1`, `m >= 1`.
pub fn min_mod_linear(n: u64, m: u64, a: u64, b: u64) -> u64 {
    assert!(n >= 1 && m >= 1, "n and m must be positive");
    let (mut n, mut m, mut a, mut b) = (n as u128, m as u128, a as u128 % m as u128, b as u128 % m as u128);
    let mut best = u128::MAX;
    loop {
        if a == 0 {
            return best.min(b) as u64;
        }
        if 2 * a > m {
            // run the progression backwards: slope m - a, starting at the last term
            b = (a * (n - 1) + b) % m;
            a = m - a;
        }
        best = best.min(b);
        // the first element of every block between two wraps is a candidate
        let k = (a * (n - 1) + b) / m;
        if k == 0 {
            return best as u64;
        }
        // values (b - m - t m) mod a for t = 0..k
        let (a2, b2) = ((a - m % a) % a, (b % a + a - m % a) % a);
        (n, m, a, b) = (k, a, a2, b2);
    }
}

/// Euler's totient function.
pub fn euler_phi(n: u64) -> u64 {
    factorize(n).into_iter().fold(n, |r, (p, _)| r / p * (p - 1))
}

/// `a ↑↑ b = a^(a^(...^a)) (b times) mod m` for `m >= 1` (`0^0 = 1`).
pub fn tetration_mod(a: u64, b: u64, m: u64) -> u64 {
    assert!(m >= 1, "modulus must be positive");
    if m == 1 {
        return 0;
    }
    if b == 0 {
        return 1;
    }
    if a == 0 {
        // 0, 1, 0, 1, ... for b = 1, 2, 3, ...
        return (b % 2) ^ 1;
    }
    if a == 1 {
        return 1;
    }
    tower(a, b, m) % m
}

/// `a ↑↑ b` for `a >= 2`, as `v` if `v < m` and `v mod m + m` otherwise.
fn tower(a: u64, b: u64, m: u64) -> u64 {
    if m == 1 {
        return 1; // the value is at least 1
    }
    if b == 0 {
        return 1;
    }
    let e = tower(a, b - 1, euler_phi(m));
    pow_capped(a, e, m)
}

fn cap(v: u128, m: u64) -> u64 {
    if v < m as u128 {
        v as u64
    } else {
        (v % m as u128) as u64 + m
    }
}

/// `a^e` capped like [`tower`].
fn pow_capped(a: u64, mut e: u64, m: u64) -> u64 {
    let (mut r, mut x) = (cap(1, m), cap(a as u128, m));
    while e > 0 {
        if e & 1 == 1 {
            r = cap(r as u128 * x as u128, m);
        }
        x = cap(x as u128 * x as u128, m);
        e >>= 1;
    }
    r
}
