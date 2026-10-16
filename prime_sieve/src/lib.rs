//! Sieves of Eratosthenes.
//!
//! * [`primes_up_to`]: all primes `<= n` with a cache-blocked odd-only segmented sieve,
//!   about `0.4` s for `n = 5 * 10^8`.
//! * [`mobius_table`]: the Moebius function alone, 1.2 bytes per number.
//! * [`mobius_and_prefix`]: the Moebius function and its prefix sums (the Mertens function) up to
//!   `n` by a linear sieve.
//! * [`totient_table`]: Euler's phi for every `i <= n` by a linear sieve.
//!
//! # Example
//!
//! ```
//! use prime_sieve::*;
//!
//! assert_eq!(primes_up_to(30), vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
//! let (mu, mertens) = mobius_and_prefix(10);
//! assert_eq!(mu, vec![0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1]);
//! assert_eq!(mertens[10], -1);
//! assert_eq!(totient_table(10), vec![0, 1, 1, 2, 2, 4, 2, 6, 4, 6, 4]);
//! ```

/// All primes `<= n` in increasing order.
///
/// # Panics
/// Panics if `n >= 2^32`.
pub fn primes_up_to(n: u64) -> Vec<u32> {
    assert!(n < 1 << 32, "n must be below 2^32");
    let n = n as usize;
    if n < 2 {
        return Vec::new();
    }
    let mut primes: Vec<u32> = vec![2];
    if n < 3 {
        return primes;
    }
    primes.reserve(estimate_pi(n));
    // odd numbers 2i + 1, i = 0..=(n - 1) / 2; segments of SEG odd numbers
    const SEG: usize = 1 << 17;
    let half = (n - 1) / 2 + 1; // number of odd numbers 1, 3, ..., <= n
    let root = (n as f64).sqrt() as usize + 1;
    // base primes up to sqrt(n), found by a plain sieve
    let mut base_sieve = vec![true; root / 2 + 1];
    let mut base: Vec<usize> = Vec::new();
    for i in 1..=root / 2 {
        let p = 2 * i + 1;
        if p > root {
            break;
        }
        if base_sieve[i] {
            base.push(p);
            let mut j = (p * p - 1) / 2;
            while j <= root / 2 {
                base_sieve[j] = false;
                j += p;
            }
        }
    }
    // next[k]: the first odd index >= the current segment start that is a multiple of base[k]
    let mut next: Vec<usize> = base.iter().map(|&p| (p * p - 1) / 2).collect();
    let mut seg = vec![true; SEG];
    let mut lo = 0usize;
    while lo < half {
        let hi = (lo + SEG).min(half);
        let len = hi - lo;
        seg[..len].fill(true);
        for (k, &p) in base.iter().enumerate() {
            let mut j = next[k];
            if j >= hi {
                continue;
            }
            while j < hi {
                seg[j - lo] = false;
                j += p;
            }
            next[k] = j;
        }
        let start = if lo == 0 { 1 } else { 0 }; // index 0 is the number 1
        for (t, &is_p) in seg[..len].iter().enumerate().skip(start) {
            if is_p {
                primes.push((2 * (lo + t) + 1) as u32);
            }
        }
        lo = hi;
    }
    primes
}

fn estimate_pi(n: usize) -> usize {
    let x = n as f64;
    (1.26 * x / x.ln()) as usize + 16
}

/// The Moebius function `mu[i]` for `i` in `0..=n` (`mu[0] = 0`) by a linear sieve, with a bitset
/// for the composite flags.
pub fn mobius_table(n: usize) -> Vec<i8> {
    let mut mu = vec![0i8; n + 1];
    let mut composite = vec![0u64; n / 64 + 1];
    let mut primes: Vec<u32> = Vec::new();
    if n >= 1 {
        mu[1] = 1;
    }
    for i in 2..=n {
        if composite[i / 64] >> (i % 64) & 1 == 0 {
            primes.push(i as u32);
            mu[i] = -1;
        }
        for &p in &primes {
            let ip = i * p as usize;
            if ip > n {
                break;
            }
            composite[ip / 64] |= 1 << (ip % 64);
            if i % p as usize == 0 {
                mu[ip] = 0;
                break;
            }
            mu[ip] = -mu[i];
        }
    }
    mu
}

/// `(mu, mertens)` for `0..=n`: `mu[i]` is the Moebius function (`mu[0] = 0`) and
/// `mertens[i] = mu[1] + ... + mu[i]`.
pub fn mobius_and_prefix(n: usize) -> (Vec<i8>, Vec<i32>) {
    let mut mu = vec![0i8; n + 1];
    let mut composite = vec![false; n + 1];
    let mut primes: Vec<u32> = Vec::new();
    if n >= 1 {
        mu[1] = 1;
    }
    for i in 2..=n {
        if !composite[i] {
            primes.push(i as u32);
            mu[i] = -1;
        }
        for &p in &primes {
            let ip = i * p as usize;
            if ip > n {
                break;
            }
            composite[ip] = true;
            if i % p as usize == 0 {
                mu[ip] = 0;
                break;
            }
            mu[ip] = -mu[i];
        }
    }
    let mut mertens = vec![0i32; n + 1];
    for i in 1..=n {
        mertens[i] = mertens[i - 1] + mu[i] as i32;
    }
    (mu, mertens)
}

/// Euler's totient `phi[i]` for `i` in `0..=n` (`phi[0] = 0`).
pub fn totient_table(n: usize) -> Vec<u32> {
    let mut phi = vec![0u32; n + 1];
    let mut composite = vec![false; n + 1];
    let mut primes: Vec<u32> = Vec::new();
    if n >= 1 {
        phi[1] = 1;
    }
    for i in 2..=n {
        if !composite[i] {
            primes.push(i as u32);
            phi[i] = i as u32 - 1;
        }
        for &p in &primes {
            let ip = i * p as usize;
            if ip > n {
                break;
            }
            composite[ip] = true;
            if i % p as usize == 0 {
                phi[ip] = phi[i] * p;
                break;
            }
            phi[ip] = phi[i] * (p - 1);
        }
    }
    phi
}
