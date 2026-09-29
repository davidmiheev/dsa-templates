//! Stress test: `ModBinomial` against Pascal's triangle and combinatorial identities.

use math::ModBinomial;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const PRIMES: [i64; 5] = [13, 101, 10_007, 998_244_353, 1_000_000_007];

fn pascal(max_n: usize, p: u64) -> Vec<Vec<u64>> {
    let mut t = vec![vec![0u64; max_n + 1]; max_n + 1];
    for n in 0..=max_n {
        t[n][0] = 1 % p;
        for r in 1..=n {
            t[n][r] = (t[n - 1][r - 1] + if r <= n - 1 { t[n - 1][r] } else { 0 }) % p;
        }
    }
    t
}

#[test]
fn matches_pascal_triangle() {
    let mut st = Stress::new("math/pascal");
    for case in 0..60 * scale() {
        let mut rng = Rng::new(case as u64);
        let p = PRIMES[case % PRIMES.len()];
        let max_n = (rng.range(1, 250) as usize).min(p as usize - 1); // table needs max_n < p
        let mb = ModBinomial::new(max_n, p);
        let t = pascal(max_n, p as u64);
        st.case();
        for n in 0..=max_n {
            for r in 0..=max_n {
                let want = if r <= n { t[n][r] as i64 } else { 0 };
                assert_eq!(mb.ncr(n, r), want, "C({n},{r}) mod {p} (max_n {max_n})");
            }
            st.ops(max_n + 1);
        }
    }
    st.done();
}

#[test]
fn identities_hold_for_large_n() {
    let mut st = Stress::new("math/identities");
    let p = 998_244_353i64;
    let max_n = 1_000_000;
    let mb = ModBinomial::new(max_n, p);
    let mut rng = Rng::new(21);
    st.case();
    for _ in 0..200_000 * scale() {
        let n = rng.range(1, max_n as i64) as usize;
        let r = rng.below(n + 1);
        assert_eq!(mb.ncr(n, r), mb.ncr(n, n - r), "symmetry n={n} r={r}");
        assert_eq!(mb.ncr(n, 0), 1);
        assert_eq!(mb.ncr(n, n), 1);
        assert_eq!(mb.ncr(n, 1), n as i64 % p);
        if r >= 1 && r <= n - 1 {
            // Pascal's rule
            assert_eq!(mb.ncr(n, r), (mb.ncr(n - 1, r - 1) + mb.ncr(n - 1, r)) % p, "pascal n={n} r={r}");
        }
        assert_eq!(mb.ncr(n, n + 1 + rng.below(5)), 0, "r > n must be 0");
        st.ops(6);
    }
    // Row sums: sum_r C(n, r) = 2^n
    for n in [0usize, 1, 2, 10, 1000, 5000] {
        let sum = (0..=n).fold(0i64, |acc, r| (acc + mb.ncr(n, r)) % p);
        let mut pow = 1i64;
        for _ in 0..n {
            pow = pow * 2 % p;
        }
        assert_eq!(sum, pow, "row sum n={n}");
        st.ops(n + 1);
    }
    st.done();
}
