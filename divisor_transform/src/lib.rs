//! Transforms on arrays indexed by `1..=n` (slot 0 is unused) along the divisibility order, each
//! in `O(n log log n)` by one pass per prime.
//!
//! * [`zeta_divisors`] / [`mobius_divisors`]: `b[k] = sum_{d | k} a[d]` and its inverse.
//! * [`zeta_multiples`] / [`mobius_multiples`]: `b[k] = sum_{k | m} a[m]` and its inverse.
//! * [`gcd_convolution`]: `c[k] = sum_{gcd(i, j) = k} a[i] b[j]` (multiples transform).
//! * [`lcm_convolution`]: `c[k] = sum_{lcm(i, j) = k} a[i] b[j]` (divisors transform; products
//!   of pairs with `lcm > n` are dropped).
//!
//! # Example
//!
//! ```
//! use divisor_transform::{gcd_convolution, lcm_convolution};
//!
//! const M: u32 = 998_244_353;
//! // a, b indexed 1..=4 (index 0 unused): all ones
//! let ones = [0, 1, 1, 1, 1];
//! assert_eq!(gcd_convolution::<M>(&ones, &ones), vec![0, 11, 3, 1, 1]);
//! assert_eq!(lcm_convolution::<M>(&ones, &ones), vec![0, 1, 3, 3, 5]);
//! ```

fn primes_up_to(n: usize) -> Vec<usize> {
    let mut sieve = vec![true; n + 1];
    let mut primes = Vec::new();
    for i in 2..=n {
        if sieve[i] {
            primes.push(i);
            let mut j = i * i;
            while j <= n {
                sieve[j] = false;
                j += i;
            }
        }
    }
    primes
}

#[inline(always)]
fn add<const M: u32>(a: u32, b: u32) -> u32 {
    let s = a + b;
    if s >= M {
        s - M
    } else {
        s
    }
}

#[inline(always)]
fn sub<const M: u32>(a: u32, b: u32) -> u32 {
    if a >= b {
        a - b
    } else {
        a + M - b
    }
}

/// `a[k] <- sum_{d | k} a[d]` for `1 <= k < a.len()`.
pub fn zeta_divisors<const M: u32>(a: &mut [u32]) {
    let n = a.len() - 1;
    for p in primes_up_to(n) {
        for i in 1..=n / p {
            a[i * p] = add::<M>(a[i * p], a[i]);
        }
    }
}

/// Inverse of [`zeta_divisors`].
pub fn mobius_divisors<const M: u32>(a: &mut [u32]) {
    let n = a.len() - 1;
    for p in primes_up_to(n) {
        for i in (1..=n / p).rev() {
            a[i * p] = sub::<M>(a[i * p], a[i]);
        }
    }
}

/// `a[k] <- sum_{k | m <= n} a[m]`.
pub fn zeta_multiples<const M: u32>(a: &mut [u32]) {
    let n = a.len() - 1;
    for p in primes_up_to(n) {
        for i in (1..=n / p).rev() {
            a[i] = add::<M>(a[i], a[i * p]);
        }
    }
}

/// Inverse of [`zeta_multiples`].
pub fn mobius_multiples<const M: u32>(a: &mut [u32]) {
    let n = a.len() - 1;
    for p in primes_up_to(n) {
        for i in 1..=n / p {
            a[i] = sub::<M>(a[i], a[i * p]);
        }
    }
}

fn pointwise<const M: u32>(a: &mut [u32], b: &[u32]) {
    for (x, &y) in a.iter_mut().zip(b) {
        *x = (*x as u64 * y as u64 % M as u64) as u32;
    }
}

/// `c[k] = sum_{gcd(i, j) = k} a[i] b[j]` for `1 <= k <= n`; slot 0 of the inputs is ignored and
/// of the output is 0. Both slices have length `n + 1`.
pub fn gcd_convolution<const M: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len(), "length mismatch");
    let (mut fa, mut fb) = (a.to_vec(), b.to_vec());
    fa[0] = 0;
    fb[0] = 0;
    zeta_multiples::<M>(&mut fa);
    zeta_multiples::<M>(&mut fb);
    pointwise::<M>(&mut fa, &fb);
    mobius_multiples::<M>(&mut fa);
    fa
}

/// `c[k] = sum_{lcm(i, j) = k} a[i] b[j]` for `1 <= k <= n`.
pub fn lcm_convolution<const M: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len(), "length mismatch");
    let (mut fa, mut fb) = (a.to_vec(), b.to_vec());
    fa[0] = 0;
    fb[0] = 0;
    zeta_divisors::<M>(&mut fa);
    zeta_divisors::<M>(&mut fb);
    pointwise::<M>(&mut fa, &fb);
    mobius_divisors::<M>(&mut fa);
    fa
}
