//! Stress tests: every routine against the direct computation for small `n`, and known values.

use prime_sieve::primes_up_to;
use sublinear::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const M: u64 = 469_762_049;

#[test]
fn lucy_and_min25_match_direct_sums() {
    let mut st = Stress::new("sublinear/lucy_min25");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = if case < 300 { 1 + case as u64 } else { rng.range(1, 200_000) as u64 };
        let primes = primes_up_to(n);
        let pi = prime_pi(n);
        let (cnt, sum) = prime_count_and_sum::<M>(n);
        for v in quotient_values(n) {
            let c = primes.iter().take_while(|&&p| p as u64 <= v).count() as u64;
            let s = primes.iter().take_while(|&&p| p as u64 <= v).fold(0u64, |a, &p| (a + p as u64) % M);
            assert_eq!(pi.get(v), c, "pi({v}) for n={n}");
            assert_eq!(cnt.get(v), c);
            assert_eq!(sum.get(v), s, "sum for n={n} v={v}");
        }
        assert_eq!(pi.ascending().len(), quotient_values(n).len());
        // multiplicative function with f(p^e) = a e + b p
        let (a, b) = (rng.below(M as usize) as u64, rng.below(M as usize) as u64);
        let pf = Quotients::from_fn(n, |x| (a * cnt.get(x) + b * sum.get(x)) % M);
        let f = |p: u64, e: u32| (a * e as u64 + b * p) % M;
        let got = multiplicative_prefix_sum::<M>(n, &primes, &pf, f);
        // direct: factorise every x <= n
        let mut want = 0u64;
        for x in 1..=n.min(20_000) {
            let mut y = x;
            let mut val = 1u64;
            let mut d = 2;
            while d * d <= y {
                if y % d == 0 {
                    let mut e = 0;
                    while y % d == 0 {
                        y /= d;
                        e += 1;
                    }
                    val = val * f(d, e) % M;
                }
                d += 1;
            }
            if y > 1 {
                val = val * f(y, 1) % M;
            }
            want = (want + val) % M;
        }
        if n <= 20_000 {
            assert_eq!(got, want, "n={n} a={a} b={b}");
        }
        st.case();
        st.ops(n as usize);
    }
    st.done();
}

#[test]
fn prime_counts_at_known_limits() {
    assert_eq!(prime_pi(1_000_000_000).get(1_000_000_000), 50_847_534);
    assert_eq!(prime_pi(10_000_000_000).get(10_000_000_000), 455_052_511);
    assert_eq!(prime_pi(100_000_000_000).get(100_000_000_000), 4_118_054_813);
}

const P: u64 = 998_244_353;

fn prefix_table(n: u64, f: &[u64]) -> Quotients<u64> {
    // f[i] for i in 1..=n (f[0] unused)
    let pre: Vec<u64> = std::iter::once(0).chain(f[1..].iter().scan(0u64, |a, &x| { *a = (*a + x) % P; Some(*a) })).collect();
    Quotients::from_fn(n, |v| pre[v as usize])
}

#[test]
fn dirichlet_product_and_inverse() {
    let mut st = Stress::new("sublinear/dirichlet");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(50_000 + case as u64);
        let n = if case < 100 { 1 + case as u64 } else { rng.range(1, 3000) as u64 };
        let nu = n as usize;
        let f: Vec<u64> = (0..=nu).map(|i| if i == 0 { 0 } else { 1 + rng.below(P as usize - 1) as u64 }).collect();
        let g: Vec<u64> = (0..=nu).map(|i| if i == 0 { 0 } else { rng.below(P as usize) as u64 }).collect();
        let mut h = vec![0u64; nu + 1];
        for i in 1..=nu {
            for j in 1..=nu / i {
                h[i * j] = (h[i * j] + f[i] * g[j]) % P;
            }
        }
        let (tf, tg, th) = (prefix_table(n, &f), prefix_table(n, &g), prefix_table(n, &h));
        let prod = dirichlet_product::<P>(&tf, &tg);
        for v in quotient_values(n) {
            assert_eq!(prod.get(v), th.get(v), "product n={n} v={v}");
        }
        assert_eq!(prod.ascending(), th.ascending());
        // the inverse of f: f * inverse = identity, whose prefix sums are all 1
        let inv = dirichlet_inverse::<P>(&tf);
        let check = dirichlet_product::<P>(&tf, &inv);
        assert!(quotient_values(n).iter().all(|&v| check.get(v) == 1), "inverse n={n}");
        let rebuilt = Quotients::from_ascending(n, &prod.ascending());
        assert_eq!(&rebuilt.small[1..], &prod.small[1..]);
        assert_eq!(&rebuilt.large[1..], &prod.large[1..]);
        st.case();
        st.ops(nu);
    }
    st.done();
}

#[test]
fn division_sums() {
    let mut st = Stress::new("sublinear/division");
    let phi = prime_sieve::totient_table(3000);
    let (mu, _) = prime_sieve::mobius_and_prefix(3000);
    let mut sq = vec![0u64; 3001];
    for i in 1..=3000usize {
        sq[i] = sq[i - 1] + (mu_free(i) as u64);
    }
    for n in 1..=3000u64 {
        let want = phi[1..=n as usize].iter().fold(0u64, |a, &x| (a + x as u64) % P);
        assert_eq!(totient_sum::<998_244_353>(n), want, "totient n={n}");
        assert_eq!(squarefree_count(n), sq[n as usize], "squarefree n={n}");
        st.case();
    }
    let _ = mu;
    // known values
    assert_eq!(squarefree_count(1_000_000_000), 607_927_124);
    // the plain formula sum_d mu(d) floor(n / d^2) for large n (the split at L is exercised)
    for n in [10_000_000_000_000u64, 99_999_999_999_999, 1_234_567_890_123_456] {
        let s = isqrt(n);
        let mu = prime_sieve::mobius_table(s as usize);
        let direct: i64 = (1..=s).map(|d| mu[d as usize] as i64 * (n / (d * d)) as i64).sum();
        assert_eq!(squarefree_count(n), direct as u64, "n={n}");
    }
    assert_eq!(totient_sum::<998_244_353>(1_000_000_000), 303_963_551_173_008_414 % P);
    let table: Vec<i32> = prime_sieve::mobius_and_prefix(1000).1;
    assert_eq!(mertens(1_000_000, &table), 212);
    st.done();
}

fn mu_free(i: usize) -> bool {
    (2..).take_while(|d| d * d <= i).all(|d| !i.is_multiple_of(d * d))
}
