//! Stress tests: `steady_ant` against the definition (min-plus of distribution matrices) and
//! `range_lis` against patience sorting, with duplicates, strict and non-strict.

use seaweed::{range_lis, steady_ant};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn random_perm(rng: &mut Rng, n: usize) -> Vec<u32> {
    let mut p: Vec<u32> = (0..n as u32).collect();
    for i in (1..n).rev() {
        p.swap(i, rng.below(i + 1));
    }
    p
}

/// Σ(i, j) = #{ r >= i : p[r] < j } for i, j in 0..=n.
fn distribution(p: &[u32]) -> Vec<Vec<i32>> {
    let n = p.len();
    let mut s = vec![vec![0; n + 1]; n + 1];
    for i in (0..n).rev() {
        for j in 0..=n {
            s[i][j] = s[i + 1][j] + (p[i] < j as u32) as i32;
        }
    }
    s
}

fn brute_product(p: &[u32], q: &[u32]) -> Vec<u32> {
    let n = p.len();
    let (a, b) = (distribution(p), distribution(q));
    let mut c = vec![vec![0; n + 1]; n + 1];
    for i in 0..=n {
        for k in 0..=n {
            c[i][k] = (0..=n).map(|j| a[i][j] + b[j][k]).min().unwrap();
        }
    }
    (0..n)
        .map(|r| (0..n).find(|&col| c[r][col + 1] - c[r][col] - c[r + 1][col + 1] + c[r + 1][col] == 1).unwrap() as u32)
        .collect()
}

fn brute_lis(a: &[i64], strict: bool) -> usize {
    let mut tails: Vec<i64> = Vec::new();
    for &x in a {
        let pos = if strict { tails.partition_point(|&t| t < x) } else { tails.partition_point(|&t| t <= x) };
        if pos == tails.len() {
            tails.push(x);
        } else {
            tails[pos] = x;
        }
    }
    tails.len()
}

#[test]
fn steady_ant_matches_definition() {
    let mut st = Stress::new("seaweed/steady_ant");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 30) as usize;
        let (p, q) = (random_perm(&mut rng, n), random_perm(&mut rng, n));
        assert_eq!(steady_ant(&p, &q), brute_product(&p, &q), "p={p:?} q={q:?}");
        st.case();
        st.ops(1);
    }
    // Structured inputs.
    for n in [1usize, 2, 3, 8, 17, 64] {
        let id: Vec<u32> = (0..n as u32).collect();
        let rev: Vec<u32> = id.iter().rev().copied().collect();
        for (p, q) in [(&id, &id), (&id, &rev), (&rev, &id), (&rev, &rev)] {
            assert_eq!(steady_ant(p, q), brute_product(p, q), "n={n}");
        }
    }
    assert!(steady_ant(&[], &[]).is_empty());
    st.done();
}

#[test]
fn product_is_associative_with_identity() {
    let mut rng = Rng::new(99);
    for _ in 0..200 * scale() {
        let n = rng.range(1, 40) as usize;
        let (a, b, c) = (random_perm(&mut rng, n), random_perm(&mut rng, n), random_perm(&mut rng, n));
        assert_eq!(steady_ant(&steady_ant(&a, &b), &c), steady_ant(&a, &steady_ant(&b, &c)));
        let id: Vec<u32> = (0..n as u32).collect();
        assert_eq!(steady_ant(&a, &id), a);
        assert_eq!(steady_ant(&id, &a), a);
    }
}

#[test]
fn range_lis_matches_patience_sorting() {
    let mut st = Stress::new("seaweed/range_lis");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(7000 + case as u64);
        let n = rng.range(0, 40) as usize;
        let values = *[2i64, 5, 1000].get(case % 3).unwrap();
        let a: Vec<i64> = if case % 3 == 2 { random_perm(&mut rng, n).iter().map(|&x| x as i64).collect() } else { (0..n).map(|_| rng.range(0, values)).collect() };
        let queries: Vec<(usize, usize)> = (0..30)
            .map(|_| {
                let l = rng.below(n + 1);
                (l, l + rng.below(n + 1 - l))
            })
            .collect();
        for strict in [true, false] {
            let got = range_lis(&a, &queries, strict);
            for (i, &(l, r)) in queries.iter().enumerate() {
                assert_eq!(got[i], brute_lis(&a[l..r], strict), "query [{l}, {r}) strict={strict} a={a:?}");
            }
        }
        st.case();
        st.ops(queries.len());
    }
    st.done();
}

#[test]
fn range_lis_large() {
    let mut st = Stress::new("seaweed/range_lis-large");
    let n = 100_000;
    let mut rng = Rng::new(5);
    let a: Vec<i64> = random_perm(&mut rng, n).iter().map(|&x| x as i64).collect();
    let queries: Vec<(usize, usize)> = (0..100_000 * scale())
        .map(|_| {
            let l = rng.below(n);
            (l, l + 1 + rng.below(n - l))
        })
        .collect();
    st.case();
    let got = range_lis(&a, &queries, true);
    for (i, &(l, r)) in queries.iter().enumerate().filter(|(_, &(l, r))| r - l <= 3000).take(40) {
        assert_eq!(got[i], brute_lis(&a[l..r], true), "query {i} = [{l}, {r})");
    }
    // Sorted and reversed sequences have trivial answers.
    let sorted: Vec<i64> = (0..n as i64).collect();
    assert_eq!(range_lis(&sorted, &[(0, n), (10, 20)], true), vec![n, 10]);
    let reversed: Vec<i64> = (0..n as i64).rev().collect();
    assert_eq!(range_lis(&reversed, &[(0, n), (10, 20)], true), vec![1, 1]);
    st.ops(queries.len());
    st.done();
}

#[test]
#[should_panic(expected = "invalid query")]
fn out_of_range_query_panics() {
    range_lis(&[1, 2, 3], &[(2, 4)], true);
}

#[test]
fn prefix_substring_lcs_matches_dp() {
    use seaweed::prefix_substring_lcs;
    let mut st = Stress::new("seaweed/lcs");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(40_000 + case as u64);
        let (m, n) = (1 + rng.below(12), 1 + rng.below(12));
        let sigma = 1 + rng.below(4);
        let s: Vec<u8> = (0..m).map(|_| rng.below(sigma) as u8).collect();
        let t: Vec<u8> = (0..n).map(|_| rng.below(sigma) as u8).collect();
        let mut queries = Vec::new();
        for a in 0..=m {
            for b in 0..=n {
                for c in b..=n {
                    queries.push((a, b, c));
                }
            }
        }
        let got = prefix_substring_lcs(&s, &t, &queries);
        for (&(a, b, c), &g) in queries.iter().zip(&got) {
            // plain dynamic programming
            let (x, y) = (&s[..a], &t[b..c]);
            let mut dp = vec![vec![0usize; y.len() + 1]; x.len() + 1];
            for i in 0..x.len() {
                for j in 0..y.len() {
                    dp[i + 1][j + 1] = if x[i] == y[j] { dp[i][j] + 1 } else { dp[i][j + 1].max(dp[i + 1][j]) };
                }
            }
            assert_eq!(g, dp[x.len()][y.len()], "s={s:?} t={t:?} query=({a},{b},{c})");
        }
        st.case();
        st.ops(queries.len());
    }
    st.done();
}
