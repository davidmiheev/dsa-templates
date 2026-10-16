//! Stress test: primality and factorisation against trial division, primitive roots by order.

use number_theory::{factorize, gcd, is_prime, mod_inverse, mul_mod, pow_mod, primitive_root};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn trial_factor(mut n: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut p = 2;
    while p * p <= n {
        if n % p == 0 {
            let mut e = 0;
            while n % p == 0 {
                n /= p;
                e += 1;
            }
            out.push((p, e));
        }
        p += 1;
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

#[test]
fn small_numbers_match_trial_division() {
    let mut st = Stress::new("number_theory/small");
    for n in 1..20_000u64 {
        let f = trial_factor(n);
        assert_eq!(factorize(n), f, "n={n}");
        assert_eq!(is_prime(n), f.len() == 1 && f[0].1 == 1, "n={n}");
        st.ops(1);
    }
    st.case();
    st.done();
}

#[test]
fn large_numbers() {
    let mut st = Stress::new("number_theory/large");
    let mut rng = Rng::new(1);
    for case in 0..300 * scale() {
        // products of random primes near 2^30..2^31 and arbitrary 64-bit numbers
        let n = match case % 3 {
            0 => rng.next_u64() | 1,
            1 => {
                let mut p = || loop {
                    let c = (rng.next_u64() >> 33) | (1 << 30);
                    if is_prime(c) {
                        return c;
                    }
                };
                p() * p()
            }
            _ => 1 + rng.next_u64() % 1_000_000_000_000,
        };
        let f = factorize(n);
        let prod = f.iter().fold(1u128, |acc, &(p, e)| acc * (p as u128).pow(e));
        assert_eq!(prod, n as u128);
        assert!(f.iter().all(|&(p, _)| is_prime(p)));
        assert!(f.windows(2).all(|w| w[0].0 < w[1].0));
        st.case();
    }
    // known values
    assert!(is_prime(18446744073709551557)); // largest prime below 2^64
    assert!(!is_prime(3215031751)); // strong pseudoprime to bases 2, 3, 5, 7
    assert_eq!(factorize(u64::MAX), vec![(3, 1), (5, 1), (17, 1), (257, 1), (641, 1), (65537, 1), (6700417, 1)]);
    st.done();
}

#[test]
fn modular_helpers_and_primitive_roots() {
    let mut rng = Rng::new(2);
    for _ in 0..1000 {
        let m = 2 + rng.next_u64() % 1_000_000_007;
        let a = rng.next_u64() % m;
        if gcd(a, m) == 1 {
            assert_eq!(mul_mod(a, mod_inverse(a, m), m), 1 % m);
        }
        assert_eq!(pow_mod(a, 0, m), 1 % m);
    }
    for p in (2..400u64).filter(|&p| is_prime(p)) {
        let g = primitive_root(p);
        let mut x = 1;
        let mut seen = std::collections::HashSet::new();
        for _ in 0..p - 1 {
            x = x * g % p;
            seen.insert(x);
        }
        assert_eq!(seen.len() as u64, p - 1, "p={p} g={g}");
    }
    assert_eq!(primitive_root(998_244_353), 3);
    assert_eq!(primitive_root(1_000_000_007), 5);
}
