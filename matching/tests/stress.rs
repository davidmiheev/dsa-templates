//! Stress test: the three matchings against exhaustive search / bitmask DP.

use matching::{blossom, hopcroft_karp, hungarian, NONE};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

/// Max matching size of a general graph by bitmask DP.
fn brute_matching(n: usize, edges: &[(usize, usize)]) -> usize {
    let mut adj = vec![0u32; n];
    for &(u, v) in edges {
        if u != v {
            adj[u] |= 1 << v;
            adj[v] |= 1 << u;
        }
    }
    let mut best = vec![0usize; 1 << n];
    for mask in 1u32..1 << n {
        let v = mask.trailing_zeros() as usize;
        let rest = mask & (mask - 1);
        let mut b = best[rest as usize]; // v unmatched
        let mut cand = adj[v] & rest;
        while cand != 0 {
            let w = cand.trailing_zeros();
            cand &= cand - 1;
            b = b.max(1 + best[(rest & !(1 << w)) as usize]);
        }
        best[mask as usize] = b;
    }
    best[(1 << n) - 1]
}

#[test]
fn bipartite_matches_general_brute_force() {
    let mut st = Stress::new("matching/hopcroft_karp");
    for case in 0..2000 * scale() {
        let mut rng = Rng::new(case as u64);
        let (l, r) = (rng.range(1, 6) as usize, rng.range(1, 6) as usize);
        let edges: Vec<(usize, usize)> = (0..rng.range(0, 14)).map(|_| (rng.below(l), rng.below(r))).collect();
        let m = hopcroft_karp(l, r, &edges);
        let as_general: Vec<(usize, usize)> = edges.iter().map(|&(a, b)| (a, l + b)).collect();
        assert_eq!(m.size, brute_matching(l + r, &as_general), "case {case}: {edges:?}");
        // Consistency of the arrays.
        let mut count = 0;
        for a in 0..l {
            if m.left_to_right[a] != NONE {
                count += 1;
                assert_eq!(m.right_to_left[m.left_to_right[a]], a);
                assert!(edges.contains(&(a, m.left_to_right[a])), "matched pair is an edge");
            }
        }
        assert_eq!(count, m.size);
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn blossom_matches_bitmask_dp() {
    let mut st = Stress::new("matching/blossom");
    for case in 0..2000 * scale() {
        let mut rng = Rng::new(case as u64 + 7);
        let n = rng.range(1, 11) as usize;
        let edges: Vec<(usize, usize)> = (0..rng.range(0, 22)).map(|_| (rng.below(n), rng.below(n))).collect();
        let m = blossom(n, &edges);
        assert_eq!(m.size, brute_matching(n, &edges), "case {case}: n={n} {edges:?}");
        for v in 0..n {
            if m.mate[v] != NONE {
                assert_eq!(m.mate[m.mate[v]], v);
                assert!(edges.contains(&(v, m.mate[v])) || edges.contains(&(m.mate[v], v)), "matched pair is an edge");
                assert_ne!(m.mate[v], v);
            }
        }
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn hungarian_matches_permutation_search() {
    fn permutations(n: usize, f: &mut impl FnMut(&[usize])) {
        let mut p: Vec<usize> = (0..n).collect();
        fn rec(k: usize, p: &mut Vec<usize>, f: &mut impl FnMut(&[usize])) {
            if k == p.len() {
                f(p);
                return;
            }
            for i in k..p.len() {
                p.swap(k, i);
                rec(k + 1, p, f);
                p.swap(k, i);
            }
        }
        rec(0, &mut p, f);
    }
    let mut st = Stress::new("matching/hungarian");
    for case in 0..800 * scale() {
        let mut rng = Rng::new(case as u64 + 99);
        let n = rng.range(1, 7) as usize;
        let amp = *[3i64, 1000, 1_000_000_000].get(case % 3).unwrap();
        let cost: Vec<Vec<i64>> = (0..n).map(|_| (0..n).map(|_| rng.range(-amp, amp)).collect()).collect();
        let mut best = i64::MAX;
        permutations(n, &mut |p| best = best.min((0..n).map(|i| cost[i][p[i]]).sum()));
        let (total, assign) = hungarian(&cost);
        assert_eq!(total, best, "case {case}");
        let mut seen = vec![false; n];
        for &c in &assign {
            assert!(!seen[c]);
            seen[c] = true;
        }
        assert_eq!((0..n).map(|i| cost[i][assign[i]]).sum::<i64>(), total);
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn large_inputs() {
    let mut st = Stress::new("matching/large");
    let mut rng = Rng::new(4);
    st.case();
    // Bipartite: a random graph with a planted perfect matching.
    let (l, m) = (100_000, 200_000);
    let mut edges: Vec<(usize, usize)> = (0..l).map(|i| (i, (i * 7919) % l)).collect();
    edges.extend((0..m - l).map(|_| (rng.below(l), rng.below(l))));
    assert_eq!(hopcroft_karp(l, l, &edges).size, l);
    // General: n = 500 complete graph has a perfect matching; the Hungarian algorithm at n = 500.
    let n = 500;
    let mut kn = Vec::new();
    for a in 0..n {
        for b in a + 1..n {
            kn.push((a, b));
        }
    }
    assert_eq!(blossom(n, &kn).size, n / 2);
    let cost: Vec<Vec<i64>> = (0..n).map(|_| (0..n).map(|_| rng.range(-1_000_000_000, 1_000_000_000)).collect()).collect();
    let (total, assign) = hungarian(&cost);
    assert_eq!(total, (0..n).map(|i| cost[i][assign[i]]).sum::<i64>());
    st.ops(l + n * n);
    st.done();
}

#[test]
fn weighted_blossom_matches_bitmask_dp() {
    use matching::max_weight_matching;
    let mut st = Stress::new("matching/weighted");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(case as u64 + 9001);
        let n = rng.range(1, 9) as usize;
        let density = rng.range(2, 10) as u64;
        let mut w = vec![vec![0i64; n]; n];
        let mut edges = Vec::new();
        for u in 0..n {
            for v in u + 1..n {
                if rng.chance(density, 10) {
                    let weight = if rng.chance(1, 2) { rng.range(1, 3) } else { rng.range(1, 20) };
                    w[u][v] = weight;
                    w[v][u] = weight;
                    edges.push((u, v, weight));
                }
            }
        }
        // best[mask] = best weight of a matching using only the vertices of `mask`
        let mut best = vec![0i64; 1 << n];
        for mask in 1u32..1 << n {
            let v = mask.trailing_zeros() as usize;
            let rest = mask & (mask - 1);
            let mut b = best[rest as usize];
            for u in 0..n {
                if rest >> u & 1 == 1 && w[v][u] > 0 {
                    b = b.max(w[v][u] + best[(rest & !(1 << u)) as usize]);
                }
            }
            best[mask as usize] = b;
        }
        let m = max_weight_matching(n, &edges);
        assert_eq!(m.weight, best[(1 << n) - 1], "case {case}: {edges:?}");
        // consistency
        let mut total = 0;
        let mut count = 0;
        for v in 0..n {
            if m.mate[v] != NONE {
                assert_eq!(m.mate[m.mate[v]], v);
                assert!(w[v][m.mate[v]] > 0, "matched pair is an edge");
                if v < m.mate[v] {
                    total += w[v][m.mate[v]];
                    count += 1;
                }
            }
        }
        assert_eq!((total, count), (m.weight, m.size));
        st.case();
    }
    st.done();
}
