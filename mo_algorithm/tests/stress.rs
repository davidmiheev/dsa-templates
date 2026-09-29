//! Stress test: Mo's algorithm against naive per-query scans, with window-invariant
//! checks on the add/remove callbacks, and an offline BIT reference at scale.

use mo_algorithm::{count_distinct_in_range, solve_blocked, Query};
use std::collections::HashSet;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn random_queries(rng: &mut Rng, n: usize, q: usize) -> Vec<Query> {
    (0..q)
        .map(|_| {
            let a = rng.below(n + 1);
            let b = rng.below(n + 1);
            Query { lo: a.min(b), hi: a.max(b) } // empty ranges (lo == hi) included
        })
        .collect()
}

#[test]
fn count_distinct_matches_naive() {
    let mut st = Stress::new("mo/distinct");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 60) as usize;
        let values = *[2usize, 5, 60].get(case % 3).unwrap();
        let a: Vec<usize> = (0..n).map(|_| rng.below(values)).collect();
        let q = rng.below(80);
        let queries = random_queries(&mut rng, n, q);
        let got = count_distinct_in_range(&a, &queries);
        assert_eq!(got.len(), queries.len(), "one answer per query, case {case}");
        for (i, q) in queries.iter().enumerate() {
            let want = a[q.lo..q.hi].iter().collect::<HashSet<_>>().len();
            assert_eq!(got[i], want, "query {i} = [{}, {}) case {case} n={n}", q.lo, q.hi);
        }
        st.case();
        st.ops(queries.len());
    }
    st.done();
}

#[test]
fn window_callbacks_are_consistent() {
    // A custom state: sum of the window, plus a membership bitmap that catches double-adds,
    // removes of absent elements, and a wrong window at `visit` time.
    struct Window {
        sum: i64,
        inside: Vec<bool>,
        visited: Vec<u32>,
    }
    let mut st = Stress::new("mo/window");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(7000 + case as u64);
        let n = rng.range(1, 50) as usize;
        let a: Vec<i64> = (0..n).map(|_| rng.range(-100, 100)).collect();
        let q = rng.below(70);
        let queries = random_queries(&mut rng, n, q);
        let mut w = Window { sum: 0, inside: vec![false; n], visited: vec![0; queries.len()] };
        let mut calls = 0usize;
        let mut wrong = None;
        solve_blocked(
            n,
            &queries,
            &mut w,
            |w, i| {
                assert!(!w.inside[i], "add({i}) while already inside, case {case}");
                w.inside[i] = true;
                w.sum += a[i];
            },
            |w, i| {
                assert!(w.inside[i], "remove({i}) while not inside, case {case}");
                w.inside[i] = false;
                w.sum -= a[i];
            },
            |w, qi| {
                let q = queries[qi];
                let ok = (0..n).all(|i| w.inside[i] == (q.lo <= i && i < q.hi));
                if !ok || w.sum != a[q.lo..q.hi].iter().sum::<i64>() {
                    wrong = Some(qi);
                }
                w.visited[qi] += 1;
                calls += 1;
            },
        );
        assert_eq!(wrong, None, "window content at visit time, case {case}");
        assert!(w.visited.iter().all(|&c| c == 1), "each query visited exactly once, case {case}");
        assert_eq!(calls, queries.len());
        st.case();
        st.ops(queries.len());
    }
    st.done();
}

#[test]
fn no_queries_means_no_callbacks() {
    let mut touched = false;
    solve_blocked(5, &[], &mut touched, |t, _| *t = true, |t, _| *t = true, |t, _| *t = true);
    assert!(!touched);
    assert!(count_distinct_in_range(&[1, 2, 3], &[]).is_empty());
}

/// Offline reference: sort queries by `hi`, keep a BIT with a 1 at the last occurrence of each value.
fn distinct_offline(a: &[usize], queries: &[Query]) -> Vec<usize> {
    let n = a.len();
    let mut bit = vec![0i64; n + 2];
    let add = |bit: &mut Vec<i64>, mut i: usize, d: i64| {
        i += 1;
        while i <= n + 1 {
            bit[i] += d;
            i += i & i.wrapping_neg();
        }
    };
    let pref = |bit: &Vec<i64>, mut i: usize| {
        let mut s = 0;
        while i > 0 {
            s += bit[i];
            i -= i & i.wrapping_neg();
        }
        s
    };
    let mut order: Vec<usize> = (0..queries.len()).collect();
    order.sort_by_key(|&i| queries[i].hi);
    let mut last = std::collections::HashMap::new();
    let (mut ans, mut pos) = (vec![0usize; queries.len()], 0usize);
    for qi in order {
        while pos < queries[qi].hi {
            if let Some(&p) = last.get(&a[pos]) {
                add(&mut bit, p, -1);
            }
            add(&mut bit, pos, 1);
            last.insert(a[pos], pos);
            pos += 1;
        }
        let q = queries[qi];
        ans[qi] = (pref(&bit, q.hi) - pref(&bit, q.lo)) as usize;
    }
    ans
}

#[test]
fn large_matches_offline_bit_reference() {
    let mut st = Stress::new("mo/large");
    let mut rng = Rng::new(31);
    // Debug builds run the same checks on a smaller input to keep `cargo test` quick.
    let (n, q) = if cfg!(debug_assertions) { (20_000, 20_000) } else { (200_000, 200_000) };
    for values in [10usize, 1000, n] {
        let a: Vec<usize> = (0..n).map(|_| rng.below(values)).collect();
        let queries = random_queries(&mut rng, n, q);
        assert_eq!(count_distinct_in_range(&a, &queries), distinct_offline(&a, &queries), "values={values}");
        st.case();
        st.ops(q);
    }
    st.done();
}
