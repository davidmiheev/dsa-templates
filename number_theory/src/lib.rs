//! 64-bit number theory: modular arithmetic helpers, deterministic Miller–Rabin, Pollard's rho
//! factorisation and primitive roots.
//!
//! * [`mul_mod`], [`pow_mod`], [`gcd`], [`mod_inverse`] (asserts invertibility).
//! * [`is_prime`]: deterministic for every `u64` (bases `2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37`).
//! * [`factorize`]: prime factorisation with multiplicities, sorted, via trial division by small
//!   primes and Brent's variant of Pollard's rho; expected `O(n^{1/4})`.
//! * [`primitive_root`]: a generator of the multiplicative group modulo a prime.
//! * [`mod_sqrt`]: a square root modulo a prime (Tonelli–Shanks), [`NONE`] if there is none.
//!
//! # Example
//!
//! ```
//! use number_theory::{factorize, is_prime, primitive_root};
//!
//! assert!(is_prime(998_244_353));
//! assert_eq!(factorize(600_851_475_143), vec![(71, 1), (839, 1), (1471, 1), (6857, 1)]);
//! assert_eq!(primitive_root(998_244_353), 3);
//! ```

/// `a * b mod m` without overflow (`m > 0`).
pub fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    (a as u128 * b as u128 % m as u128) as u64
}

/// `base^exp mod m` (`m > 0`; the result for `m = 1` is 0).
pub fn pow_mod(base: u64, mut exp: u64, m: u64) -> u64 {
    let mut result = 1 % m;
    let mut b = base % m;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mul_mod(result, b, m);
        }
        b = mul_mod(b, b, m);
        exp >>= 1;
    }
    result
}

/// Greatest common divisor (`gcd(0, 0) = 0`).
pub fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// The inverse of `a` modulo `m` in `0..m`.
///
/// # Panics
/// Panics if `gcd(a, m) != 1`.
pub fn mod_inverse(a: u64, m: u64) -> u64 {
    let (mut r0, mut r1) = (m as i128, (a % m) as i128);
    let (mut t0, mut t1) = (0i128, 1i128);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (t0, t1) = (t1, t0 - q * t1);
    }
    assert!(r0 == 1, "{a} is not invertible modulo {m}");
    t0.rem_euclid(m as i128) as u64
}

/// Deterministic primality test for any `u64`.
pub fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for &p in &BASES {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    'bases: for &a in &BASES {
        let mut x = pow_mod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mul_mod(x, x, n);
            if x == n - 1 {
                continue 'bases;
            }
        }
        return false;
    }
    true
}

/// A nontrivial divisor of the composite `n` (Brent's cycle detection with batched gcds).
fn rho(n: u64) -> u64 {
    if n.is_multiple_of(2) {
        return 2;
    }
    let mut c = 1u64;
    loop {
        let f = |x: u64| (mul_mod(x, x, n) + c) % n;
        let (mut x, mut y, mut q, mut g) = (2u64, 2u64, 1u64, 1u64);
        let (mut r, m) = (1usize, 128usize);
        let mut ys = y;
        while g == 1 {
            x = y;
            for _ in 0..r {
                y = f(y);
            }
            let mut k = 0;
            while k < r && g == 1 {
                ys = y;
                for _ in 0..m.min(r - k) {
                    y = f(y);
                    q = mul_mod(q, x.abs_diff(y), n);
                }
                g = gcd(q, n);
                k += m;
            }
            r *= 2;
        }
        if g == n {
            g = 1;
            while g == 1 {
                ys = f(ys);
                g = gcd(x.abs_diff(ys), n);
            }
        }
        if g != n {
            return g;
        }
        c += 1;
    }
}

fn split(n: u64, out: &mut Vec<u64>) {
    if n == 1 {
        return;
    }
    if is_prime(n) {
        out.push(n);
        return;
    }
    let d = rho(n);
    split(d, out);
    split(n / d, out);
}

/// Prime factorisation `[(p, e), ...]` with increasing primes; empty for `n = 1`.
///
/// # Panics
/// Panics for `n = 0`.
pub fn factorize(mut n: u64) -> Vec<(u64, u32)> {
    assert!(n >= 1, "cannot factorise 0");
    let mut primes = Vec::new();
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] {
        while n.is_multiple_of(p) {
            primes.push(p);
            n /= p;
        }
    }
    split(n, &mut primes);
    primes.sort_unstable();
    let mut result: Vec<(u64, u32)> = Vec::new();
    for p in primes {
        match result.last_mut() {
            Some((q, e)) if *q == p => *e += 1,
            _ => result.push((p, 1)),
        }
    }
    result
}

/// The smallest primitive root of the prime `p`.
///
/// # Panics
/// Panics if `p` is not prime.
pub fn primitive_root(p: u64) -> u64 {
    assert!(is_prime(p), "{p} is not prime");
    if p == 2 {
        return 1;
    }
    let factors = factorize(p - 1);
    (2..p)
        .find(|&g| factors.iter().all(|&(q, _)| pow_mod(g, (p - 1) / q, p) != 1))
        .expect("every prime has a primitive root")
}

/// Returned by [`mod_sqrt`] when the argument is a non-residue.
pub const NONE: u64 = u64::MAX;

/// A square root of `a` modulo the odd or even prime `p` (Tonelli–Shanks), or [`NONE`] if `a` is
/// not a quadratic residue. The root returned is one of the two (`r` and `p - r`).
///
/// # Panics
/// Panics if `p` is not prime (checked only cheaply, via `p >= 2`).
pub fn mod_sqrt(a: u64, p: u64) -> u64 {
    assert!(p >= 2, "modulus must be a prime");
    let a = a % p;
    if a == 0 || p == 2 {
        return a;
    }
    if pow_mod(a, (p - 1) / 2, p) != 1 {
        return NONE;
    }
    let (mut q, mut s) = (p - 1, 0u32);
    while q.is_multiple_of(2) {
        q /= 2;
        s += 1;
    }
    let z = (2..p).find(|&z| pow_mod(z, (p - 1) / 2, p) == p - 1).expect("a non-residue exists");
    let (mut m, mut c) = (s, pow_mod(z, q, p));
    let (mut t, mut r) = (pow_mod(a, q, p), pow_mod(a, q.div_ceil(2), p));
    while t != 1 {
        let mut i = 0;
        let mut x = t;
        while x != 1 {
            x = mul_mod(x, x, p);
            i += 1;
        }
        let b = pow_mod(c, 1u64 << (m - i - 1), p);
        m = i;
        c = mul_mod(b, b, p);
        t = mul_mod(t, c, p);
        r = mul_mod(r, b, p);
    }
    r
}
