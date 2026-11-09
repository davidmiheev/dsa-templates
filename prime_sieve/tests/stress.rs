//! Stress test: the sieves against trial division and the definitions.

use prime_sieve::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn is_prime(n: u64) -> bool {
    n >= 2 && (2..).take_while(|d| d * d <= n).all(|d| !n.is_multiple_of(d))
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

#[test]
fn sieves_match_definitions() {
    let mut st = Stress::new("prime_sieve");
    for n in 0..300usize {
        let want: Vec<u32> = (0..=n as u64).filter(|&i| is_prime(i)).map(|i| i as u32).collect();
        assert_eq!(primes_up_to(n as u64), want, "n={n}");
        let phi = totient_table(n);
        let (mu, mertens) = mobius_and_prefix(n);
        assert_eq!(mobius_table(n), mu);
        for (i, &ph) in phi.iter().enumerate().skip(1) {
            assert_eq!(ph as usize, (1..=i).filter(|&j| gcd(i, j) == 1).count());
            // sum of mu over divisors is [i == 1]
            let s: i32 = (1..=i).filter(|d| i % d == 0).map(|d| mu[d] as i32).sum();
            assert_eq!(s, (i == 1) as i32);
        }
        assert_eq!(mertens.last().copied().unwrap_or(0), (1..=n).map(|i| mu[i] as i32).sum::<i32>());
        st.case();
        st.ops(n);
    }
    // segment boundaries and a few large limits against known prime counts
    for (n, pi) in [(100u64, 25usize), (1_000, 168), (131_071, 12_251), (1_000_000, 78_498), (10_000_000, 664_579)] {
        let p = primes_up_to(n);
        assert_eq!(p.len(), pi, "n={n}");
        assert!(p.windows(2).all(|w| w[0] < w[1]));
        assert!(p.last().is_none_or(|&l| l as u64 <= n));
    }
    let mut rng = Rng::new(1);
    for _ in 0..30 * scale() {
        let n = rng.range(2, 3_000_000) as u64;
        let p = primes_up_to(n);
        let q = rng.below(p.len());
        assert!(is_prime(p[q] as u64));
        // consecutive primes have no prime between them
        if q + 1 < p.len() {
            assert!((p[q] as u64 + 1..p[q + 1] as u64).all(|x| !is_prime(x)));
        }
        st.case();
    }
    assert_eq!(primes_up_to((1 << 32) - 1).len(), 203_280_221);
    st.done();
}
