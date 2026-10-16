//! Stress tests against dynamic programming / brute force definitions.
#![allow(clippy::needless_range_loop, clippy::int_plus_one)]

use combinatorics::*;
use modint::StaticModInt;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const P: u32 = 998_244_353;
type Z = StaticModInt<P>;

fn vals(v: &[Z]) -> Vec<u64> {
    v.iter().map(|x| x.val() as u64).collect()
}

#[test]
fn block_factorials() {
    let mut st = Stress::new("combinatorics/factorial");
    for block in [256usize, 1024] {
        let table = BlockFactorials::<P>::new(block);
        let mut rng = Rng::new(block as u64);
        for _ in 0..10 * scale() {
            let n = rng.below(2_000_000) as u64;
            let want = (1..=n).fold(Z::one(), |a, i| a * Z::new(i));
            assert_eq!(table.fact(n), want, "block {block} n={n}");
            st.case();
        }
        // Wilson: (p-1)! = -1, (p-2)! = 1, (p-3)! = (p-1)/2... = 2^{-1} * ... check (p-1)!/(p-1) = 1
        assert_eq!(table.fact(P as u64 - 1), -Z::one());
        assert_eq!(table.fact(P as u64 - 2), Z::one());
        assert_eq!(table.fact(P as u64 - 3), Z::new(P as u64 - 2).inv());
        // around block boundaries and the end
        for n in [0u64, 1, block as u64 - 1, block as u64, block as u64 + 1, 3 * block as u64 / 2, P as u64 - 1 - block as u64] {
            let direct = if n < 3_000_000 { (1..=n).fold(Z::one(), |a, i| a * Z::new(i)) } else { table.fact(n) };
            assert_eq!(table.fact(n), direct);
        }
        // n! (n+1) = (n+1)!
        for _ in 0..200 {
            let n = rng.below(P as usize - 2) as u64;
            assert_eq!(table.fact(n) * Z::new(n + 1), table.fact(n + 1));
        }
    }
    st.done();
}

fn pascal_mod(n: usize, m: u64) -> Vec<Vec<u64>> {
    let mut t = vec![vec![0u64; n + 1]; n + 1];
    for i in 0..=n {
        t[i][0] = 1 % m;
        for j in 1..=i {
            t[i][j] = (t[i - 1][j - 1] + if j <= i - 1 { t[i - 1][j] } else { 0 }) % m;
        }
    }
    t
}

#[test]
fn binomials() {
    let mut st = Stress::new("combinatorics/binomial");
    let tri = pascal_mod(120, u64::MAX / 4);
    // exact values in u128 for n <= 120 (C(120, 60) < 2^116)
    let mut exact = vec![vec![0u128; 121]; 121];
    for i in 0..=120 {
        exact[i][0] = 1;
        for j in 1..=i {
            exact[i][j] = exact[i - 1][j - 1] + if j < i { exact[i - 1][j] } else { 0 };
        }
    }
    let _ = tri;
    let mut ms: Vec<u64> = (1..=70).collect();
    ms.extend([97, 128, 243, 625, 1000, 4096, 6561, 9973, 30_030, 65_536, 510_510, 999_983, 1_000_000]);
    for &m in &ms {
        let b = BinomialMod::new(m);
        for n in (0..=120).step_by(7).chain([119, 120]) {
            for k in 0..=n + 1 {
                let want = if k > n { 0 } else { (exact[n][k] % m as u128) as u64 };
                assert_eq!(b.binom(n as u64, k as u64), want, "C({n},{k}) mod {m}");
            }
        }
        st.case();
    }
    // Lucas / Kummer checks for huge n: C(2^60, 2^59) is divisible by 2 exactly once... use identities
    let b = BinomialMod::new(1_000_000);
    let mut rng = Rng::new(3);
    for _ in 0..2000 * scale() {
        let n = 1 + (rng.next_u64() >> 2) % 1_000_000_000_000_000_000;
        let k = rng.next_u64() % (n + 1);
        // Pascal: C(n, k) = C(n-1, k-1) + C(n-1, k)
        if k >= 1 && k <= n - 1 {
            assert_eq!(b.binom(n, k), (b.binom(n - 1, k - 1) + b.binom(n - 1, k)) % 1_000_000, "n={n} k={k}");
        }
        assert_eq!(b.binom(n, k), b.binom(n, n - k));
    }
    // prime modulus tables
    for &m in &[2u64, 3, 7, 101, 998_244_353, 1_073_741_789] {
        let limit = 60.min(m as usize - 1);
        let t = BinomialPrime::new(m, limit);
        for n in 0..=limit {
            for k in 0..=limit {
                let want = if k > n { 0 } else { (exact[n][k] % m as u128) as u64 };
                assert_eq!(t.binom(n, k), want);
            }
        }
    }
    st.done();
}

#[test]
fn q_binomials() {
    let mut st = Stress::new("combinatorics/qbinomial");
    for &m in &[2u64, 3, 5, 7, 11, 13, 29, 31, 101] {
        // q-Pascal: [n,k] = [n-1,k-1] + q^k [n-1,k]
        let nmax = (m as usize - 1).min(90);
        for q in 0..m {
            let mut t = vec![vec![0u64; nmax + 2]; nmax + 1];
            for n in 0..=nmax {
                t[n][0] = 1;
                for k in 1..=n {
                    let qk = (0..k).fold(1u64, |a, _| a * q % m);
                    t[n][k] = (t[n - 1][k - 1] + qk * if k < n { t[n - 1][k] } else { 0 }) % m;
                }
            }
            let tab = QBinomialPrime::new(m, q, nmax);
            for n in 0..=nmax {
                for k in 0..=n + 1 {
                    let want = if k > n { 0 } else { t[n][k] };
                    assert_eq!(tab.binom(n, k), want, "[{n} {k}]_{q} mod {m}");
                }
            }
            st.case();
        }
    }
    st.done();
}

fn stirling_tables(n: usize, p: u64) -> (Vec<Vec<u64>>, Vec<Vec<u64>>) {
    // signed first kind, second kind modulo p
    let mut s1 = vec![vec![0u64; n + 2]; n + 1];
    let mut s2 = vec![vec![0u64; n + 2]; n + 1];
    s1[0][0] = 1 % p;
    s2[0][0] = 1 % p;
    for i in 1..=n {
        for k in 1..=i {
            s1[i][k] = (s1[i - 1][k - 1] + (p - (i as u64 - 1) % p) * s1[i - 1][k]) % p;
            s2[i][k] = (s2[i - 1][k - 1] + k as u64 % p * s2[i - 1][k]) % p;
        }
    }
    (s1, s2)
}

#[test]
fn stirling_numbers() {
    let mut st = Stress::new("combinatorics/stirling");
    let (s1, s2) = stirling_tables(70, P as u64);
    for n in [0usize, 1, 2, 3, 5, 8, 17, 31, 64, 70] {
        assert_eq!(vals(&stirling1_row::<P>(n)), s1[n][..=n].to_vec(), "s1 row {n}");
        assert_eq!(vals(&stirling2_row::<P>(n)), s2[n][..=n].to_vec(), "s2 row {n}");
        for k in 0..=n {
            let a = stirling1_fixed_k::<P>(n, k);
            let b = stirling2_fixed_k::<P>(n, k);
            for m in k..=n {
                assert_eq!(a[m - k].val() as u64, s1[m][k], "s1({m},{k})");
                assert_eq!(b[m - k].val() as u64, s2[m][k], "S({m},{k})");
            }
        }
        st.case();
    }
    // small primes, all (n, k) with n up to 140
    for &p in &[2u64, 3, 5, 7, 11, 13, 17, 31] {
        let (t1, t2) = stirling_tables(140, p);
        let queries: Vec<(u64, u64)> = (0..=140u64).flat_map(|n| (0..=n).map(move |k| (n, k))).collect();
        let a1 = stirling1_small_prime(p, &queries);
        let a2 = stirling2_small_prime(p, &queries);
        for (i, &(n, k)) in queries.iter().enumerate() {
            assert_eq!(a1[i], t1[n as usize][k as usize], "s1({n},{k}) mod {p}");
            assert_eq!(a2[i], t2[n as usize][k as usize], "S({n},{k}) mod {p}");
        }
        st.case();
    }
    // huge n: periodicity checks. S(n, k) mod p is a function of n - k modulo p - 1 once n - k is
    // large enough only through binomials, so compare n and n + p (p - 1)... via a structural
    // identity: S(n + 1, k) = S(n, k - 1) + k S(n, k)
    let p = 4999u64;
    let mut rng = Rng::new(8);
    for _ in 0..300 {
        let k = 1 + rng.next_u64() % 1_000_000_000_000_000;
        let n = k + rng.next_u64() % 1_000_000_000_000_000;
        let qs = [(n + 1, k), (n, k - 1), (n, k)];
        let r = stirling2_small_prime(p, &qs);
        assert_eq!(r[0], (r[1] + k % p * r[2]) % p, "S recurrence n={n} k={k}");
        // s(n+1, k) = s(n, k-1) - n s(n, k)
        let r = stirling1_small_prime(p, &qs);
        assert_eq!(r[0], (r[1] + p - n % p * r[2] % p) % p, "s recurrence n={n} k={k}");
    }
    st.done();
}

#[test]
fn sequences() {
    let mut st = Stress::new("combinatorics/sequences");
    // Bell numbers by the triangle
    let n = 40;
    let mut tri = vec![vec![0u64; n + 2]; n + 2];
    tri[0][0] = 1;
    for i in 1..=n {
        tri[i][0] = tri[i - 1][i - 1];
        for j in 1..=i {
            tri[i][j] = (tri[i][j - 1] + tri[i - 1][j - 1]) % P as u64;
        }
    }
    let bell: Vec<u64> = (0..=n).map(|i| tri[i][0]).collect();
    assert_eq!(vals(&bell_numbers::<P>(n)), bell);
    // partitions by the pentagonal recurrence-free DP (coin change)
    let mut p = vec![0u64; 301];
    p[0] = 1;
    for part in 1..=300 {
        for s in part..=300 {
            p[s] = (p[s] + p[s - part]) % P as u64;
        }
    }
    assert_eq!(vals(&partition_numbers::<P>(300)), p);
    assert_eq!(partition_numbers::<P>(0).len(), 1);
    for case in 0..300 * scale() {
        let mut rng = Rng::new(case as u64);
        // subset sums by knapsack
        let target = 1 + rng.below(40);
        let cnt = rng.below(25);
        let s: Vec<usize> = (0..cnt).map(|_| 1 + rng.below(target)).collect();
        let mut dp = vec![0u64; target + 1];
        dp[0] = 1;
        for &x in &s {
            for t in (x..=target).rev() {
                dp[t] = (dp[t] + dp[t - x]) % P as u64;
            }
        }
        assert_eq!(vals(&subset_sum_counts::<P>(&s, target)), dp[1..].to_vec());
        // derangements
        let m = 1 + rng.next_u64() % 1_000_000_000;
        let d = derangements_mod(10, m);
        let exact = [0u64, 1, 2, 9, 44, 265, 1854, 14833, 133496, 1334961];
        assert_eq!(d, exact.iter().map(|&x| x % m).collect::<Vec<_>>());
        // distinct subsequences
        let len = rng.below(14);
        let a: Vec<u32> = (0..len).map(|_| rng.below(4) as u32).collect();
        let mut set = std::collections::HashSet::new();
        for mask in 1u32..1 << len {
            set.insert((0..len).filter(|&i| mask >> i & 1 == 1).map(|i| a[i]).collect::<Vec<_>>());
        }
        assert_eq!(distinct_subsequences::<P>(&a), set.len() as u64);
        st.case();
    }
    assert_eq!(derangements_mod(5, 1), vec![0; 5]);
    st.done();
}

#[test]
fn increasing_sequences_between_bounds() {
    let mut st = Stress::new("combinatorics/lattice");
    for case in 0..4000 * scale() {
        let mut rng = Rng::new(7_000 + case as u64);
        let n = 1 + rng.below(if case % 10 == 0 { 40 } else { 8 });
        let m = 1 + rng.below(if case % 10 == 0 { 40 } else { 9 });
        let mut a = vec![0usize; n];
        let mut b = vec![0usize; n];
        for i in 0..n {
            let x = rng.below(m);
            let y = x + 1 + rng.below(m - x);
            a[i] = if rng.chance(1, 4) { 0 } else { x };
            b[i] = if rng.chance(1, 4) { m } else { y.max(a[i] + 1).min(m) };
            if b[i] <= a[i] {
                a[i] = b[i] - 1;
            }
        }
        // naive dp over the values of the sequence
        let mut f = vec![0u64; m];
        for h in a[0]..b[0] {
            f[h] = 1;
        }
        for i in 1..n {
            let mut g = vec![0u64; m];
            let mut run = 0u64;
            for h in 0..m {
                run = (run + f[h]) % P as u64;
                if h >= a[i] && h < b[i] {
                    g[h] = run;
                }
            }
            f = g;
        }
        let want = f.iter().fold(0u64, |s, &x| (s + x) % P as u64);
        for base in [1usize, 2, 6, 1024] {
            assert_eq!(count_with_base::<P>(m, &a, &b, base) as u64, want, "n={n} m={m} a={a:?} b={b:?} base={base}");
        }
        st.case();
        st.ops(n * m);
    }
    // a wide instance exercises the convolution blocks: unconstrained sequences give C(n + m - 1, n)
    let (n, m) = (300, 200);
    let got = count_increasing_sequences::<P>(m, &vec![0; n], &vec![m; n]);
    let t = BinomialPrime::new(P as u64, n + m);
    assert_eq!(got as u64, t.binom(n + m - 1, n));
    st.done();
}
