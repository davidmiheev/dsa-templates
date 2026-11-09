//! Stress tests for `tree_algos` against brute-force definitions.

use std::collections::HashMap;
use tree_algos::{cartesian_tree, rooted_isomorphism_classes, tree_diameter, NONE};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn cartesian_tree_matches_recursive_definition() {
    fn build(a: &[u32], lo: usize, hi: usize, parent_of: usize, out: &mut [usize]) {
        if lo >= hi {
            return;
        }
        // Leftmost minimum.
        let m = (lo..hi).min_by_key(|&i| (a[i], i)).unwrap();
        out[m] = parent_of;
        build(a, lo, m, m, out);
        build(a, m + 1, hi, m, out);
    }
    let mut st = Stress::new("tree_algos/cartesian");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(0, 30) as usize;
        let a: Vec<u32> = (0..n).map(|_| rng.below(if case % 2 == 0 { 5 } else { 1000 }) as u32).collect();
        let mut want = vec![NONE; n];
        build(&a, 0, n, NONE, &mut want);
        assert_eq!(cartesian_tree(&a), want, "a={a:?}");
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn diameter_matches_all_pairs() {
    let mut st = Stress::new("tree_algos/diameter");
    for case in 0..1000 * scale() {
        let mut rng = Rng::new(case as u64 + 5);
        let n = rng.range(1, 25) as usize;
        let edges: Vec<(usize, usize, u64)> = (1..n).map(|i| (rng.below(i), i, rng.below(10) as u64 + 1)).collect();
        // Floyd-Warshall.
        let mut d = vec![vec![u64::MAX / 4; n]; n];
        for i in 0..n {
            d[i][i] = 0;
        }
        for &(u, v, w) in &edges {
            d[u][v] = w;
            d[v][u] = w;
        }
        for k in 0..n {
            for i in 0..n {
                for j in 0..n {
                    d[i][j] = d[i][j].min(d[i][k] + d[k][j]);
                }
            }
        }
        let best = d.iter().flatten().copied().max().unwrap();
        let (len, path) = tree_diameter(n, &edges);
        assert_eq!(len, best, "n={n} edges={edges:?}");
        assert_eq!(d[path[0]][*path.last().unwrap()], len);
        // The path is a real path whose edge weights add up.
        let w: HashMap<(usize, usize), u64> = edges.iter().flat_map(|&(u, v, w)| [((u, v), w), ((v, u), w)]).collect();
        let total: u64 = path.windows(2).map(|p| w[&(p[0], p[1])]).sum();
        assert_eq!(total, len);
        st.case();
        st.ops(1);
    }
    st.done();
}

fn canonical(parent: &[usize], children: &[Vec<usize>], v: usize) -> String {
    let mut parts: Vec<String> = children[v].iter().map(|&c| canonical(parent, children, c)).collect();
    parts.sort();
    format!("({})", parts.concat())
}

#[test]
fn isomorphism_classes_match_canonical_strings() {
    let mut st = Stress::new("tree_algos/isomorphism");
    for case in 0..1000 * scale() {
        let mut rng = Rng::new(case as u64 + 11);
        let n = rng.range(1, 30) as usize;
        let shape = case % 3;
        let mut parent = vec![0usize; n];
        for i in 1..n {
            parent[i] = match shape {
                0 => rng.below(i),
                1 => i.saturating_sub(1 + rng.below(2)),
                _ => (i - 1) / 2,
            };
        }
        // Shuffle labels so the root is not always 0 (keep the root at a random label).
        let perm = {
            let mut p: Vec<usize> = (0..n).collect();
            for i in (1..n).rev() {
                p.swap(i, rng.below(i + 1));
            }
            p
        };
        let mut par = vec![NONE; n];
        for i in 1..n {
            par[perm[i]] = perm[parent[i]];
        }
        let root = perm[0];
        let mut children = vec![vec![]; n];
        for v in 0..n {
            if v != root {
                children[par[v]].push(v);
            }
        }
        let canon: Vec<String> = (0..n).map(|v| canonical(&par, &children, v)).collect();
        let (k, class) = rooted_isomorphism_classes(&par, root);
        assert_eq!(k, class.iter().max().unwrap() + 1);
        for u in 0..n {
            for v in 0..n {
                assert_eq!(class[u] == class[v], canon[u] == canon[v], "case {case} u={u} v={v}");
            }
        }
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn large_inputs() {
    let mut st = Stress::new("tree_algos/large");
    let n = 500_000;
    let mut rng = Rng::new(1);
    st.case();
    let a: Vec<u32> = (0..n).map(|_| rng.next_u64() as u32).collect();
    let p = cartesian_tree(&a);
    assert_eq!(p.iter().filter(|&&x| x == NONE).count(), 1);
    let edges: Vec<(usize, usize, u64)> = (1..n).map(|i| (i - 1, i, 1)).collect();
    assert_eq!(tree_diameter(n, &edges).0, (n - 1) as u64);
    let parent: Vec<usize> = (0..n).map(|i| if i == 0 { NONE } else { rng.below(i) }).collect();
    let (k, _) = rooted_isomorphism_classes(&parent, 0);
    assert!(k > 1);
    let chain: Vec<usize> = (0..n).map(|i| if i == 0 { NONE } else { i - 1 }).collect();
    assert_eq!(rooted_isomorphism_classes(&chain, 0).0, n);
    st.ops(n);
    st.done();
}

#[test]
#[should_panic(expected = "not connected")]
fn disconnected_diameter_panics() {
    tree_diameter(4, &[(0, 1, 1u32), (0, 1, 1), (2, 3, 1)]);
}

mod reroot {
    use super::*;
    use tree_algos::rerooting;

    const P: u64 = 998_244_353;

    #[test]
    fn eccentricity_matches_all_pairs_distances() {
        let mut st = Stress::new("tree_algos/rerooting-ecc");
        for case in 0..1000 * scale() {
            let mut rng = Rng::new(case as u64 + 77);
            let n = rng.range(1, 25) as usize;
            let mut edges: Vec<(usize, usize)> = Vec::new();
            for i in 1..n {
                let p = rng.below(i);
                edges.push(if rng.chance(1, 2) { (i, p) } else { (p, i) });
            }
            let mut d = vec![vec![usize::MAX / 4; n]; n];
            for i in 0..n {
                d[i][i] = 0;
            }
            for &(u, v) in &edges {
                d[u][v] = 1;
                d[v][u] = 1;
            }
            for k in 0..n {
                for i in 0..n {
                    for j in 0..n {
                        d[i][j] = d[i][j].min(d[i][k] + d[k][j]);
                    }
                }
            }
            let got = rerooting(n, &edges, 0usize, |a, b| *a.max(b), |&x, _, _, _| x + 1, |&x, _| x);
            for r in 0..n {
                assert_eq!(got[r], *d[r].iter().max().unwrap(), "case {case} root {r}");
            }
            st.case();
            st.ops(1);
        }
        st.done();
    }

    /// Σ_y P(x, y) with P(x, y) = f_{e0}(f_{e1}(... f_{ek}(a_y))): a *non-commutative-looking*
    /// DP whose per-edge map is affine; checked against enumerating all paths.
    #[test]
    fn affine_path_sums_match_enumeration() {
        let mut st = Stress::new("tree_algos/rerooting-affine");
        for case in 0..600 * scale() {
            let mut rng = Rng::new(case as u64 + 9);
            let n = rng.range(1, 15) as usize;
            let edges: Vec<(usize, usize)> = (1..n).map(|i| (rng.below(i), i)).collect();
            let a: Vec<u64> = (0..n).map(|_| rng.below(P as usize) as u64).collect();
            let f: Vec<(u64, u64)> = (0..n.saturating_sub(1)).map(|_| (1 + rng.below(P as usize - 1) as u64, rng.below(P as usize) as u64)).collect();
            // Brute force by DFS from every x.
            let mut adj = vec![vec![]; n];
            for (e, &(u, v)) in edges.iter().enumerate() {
                adj[u].push((v, e));
                adj[v].push((u, e));
            }
            fn dfs(adj: &[Vec<(usize, usize)>], f: &[(u64, u64)], a: &[u64], v: usize, p: usize) -> u64 {
                // Sum over y in the component of v (away from p) of P(v, y).
                let mut s = a[v];
                for &(u, e) in &adj[v] {
                    if u != p {
                        // P(v, y) = f_e(P(u, y)) summed over y: needs the count too.
                        let (su, cu) = pair(adj, f, a, u, v);
                        s = (s + f[e].0 * su % P + f[e].1 * cu) % P;
                    }
                }
                s
            }
            fn pair(adj: &[Vec<(usize, usize)>], f: &[(u64, u64)], a: &[u64], v: usize, p: usize) -> (u64, u64) {
                let mut cnt = 1u64;
                let mut s = a[v];
                for &(u, e) in &adj[v] {
                    if u != p {
                        let (su, cu) = pair(adj, f, a, u, v);
                        s = (s + f[e].0 * su % P + f[e].1 * cu) % P;
                        cnt += cu;
                    }
                }
                (s, cnt)
            }
            let got = rerooting(
                n,
                &edges,
                (0u64, 0u64),
                |x, y| ((x.0 + y.0) % P, x.1 + y.1),
                |&(s, k), e, _, _| ((f[e].0 * s + f[e].1 * (k % P)) % P, k),
                |&(s, k), v| ((s + a[v]) % P, k + 1),
            );
            for x in 0..n {
                assert_eq!(got[x].0, dfs(&adj, &f, &a, x, usize::MAX), "case {case} root {x}");
            }
            st.case();
            st.ops(1);
        }
        st.done();
    }

    #[test]
    fn large_path_and_star() {
        let n = 500_000;
        let path: Vec<(usize, usize)> = (1..n).map(|i| (i - 1, i)).collect();
        let ecc = rerooting(n, &path, 0usize, |a, b| *a.max(b), |&x, _, _, _| x + 1, |&x, _| x);
        assert_eq!((ecc[0], ecc[n / 2], ecc[n - 1]), (n - 1, n / 2, n - 1));
        let star: Vec<(usize, usize)> = (1..n).map(|i| (0, i)).collect();
        let ecc = rerooting(n, &star, 0usize, |a, b| *a.max(b), |&x, _, _, _| x + 1, |&x, _| x);
        assert_eq!((ecc[0], ecc[1]), (1, 2));
    }
}

mod schedule {
    use super::*;
    use tree_algos::min_inversion_order;

    fn cost(order: &[usize], c: &[u64], d: &[u64]) -> u128 {
        let (mut sd, mut x) = (0u128, 0u128);
        for &v in order {
            x += c[v] as u128 * sd;
            sd += d[v] as u128;
        }
        x
    }

    /// Minimum cost over all orders with parents before children (enumerated).
    fn brute(parent: &[usize], c: &[u64], d: &[u64]) -> u128 {
        fn rec(parent: &[usize], placed: &mut Vec<bool>, order: &mut Vec<usize>, c: &[u64], d: &[u64], best: &mut u128) {
            let n = parent.len();
            if order.len() == n {
                *best = (*best).min(cost(order, c, d));
                return;
            }
            for v in 0..n {
                if !placed[v] && (parent[v] == NONE || placed[parent[v]]) {
                    placed[v] = true;
                    order.push(v);
                    rec(parent, placed, order, c, d, best);
                    order.pop();
                    placed[v] = false;
                }
            }
        }
        let mut best = u128::MAX;
        rec(parent, &mut vec![false; parent.len()], &mut vec![], c, d, &mut best);
        best
    }

    #[test]
    fn matches_exhaustive_search() {
        let mut st = Stress::new("tree_algos/schedule");
        for case in 0..1500 * scale() {
            let mut rng = Rng::new(case as u64 + 31);
            let n = rng.range(1, 8) as usize;
            let amp = *[1usize, 3, 20].get(case % 3).unwrap();
            let mut parent: Vec<usize> = vec![NONE; n];
            for i in 1..n {
                parent[i] = rng.below(i);
            }
            let (c, d): (Vec<u64>, Vec<u64>) = (0..n).map(|_| (rng.below(amp + 1) as u64, rng.below(amp + 1) as u64)).unzip();
            let (order, got) = min_inversion_order(&parent, 0, &c, &d);
            assert_eq!(got, cost(&order, &c, &d));
            assert_eq!(got, brute(&parent, &c, &d), "parent={parent:?} c={c:?} d={d:?}");
            // Valid topological order.
            let mut pos = vec![0; n];
            for (i, &v) in order.iter().enumerate() {
                pos[v] = i;
            }
            assert!((1..n).all(|v| pos[parent[v]] < pos[v]));
            st.case();
            st.ops(1);
        }
        st.done();
    }

    #[test]
    fn large_inputs_run_fast() {
        let n = 200_000;
        let mut rng = Rng::new(2);
        let parent: Vec<usize> = (0..n).map(|i| if i == 0 { NONE } else { rng.below(i) }).collect();
        let c: Vec<u64> = (0..n).map(|_| rng.below(10_000) as u64).collect();
        let d: Vec<u64> = (0..n).map(|_| rng.below(10_000) as u64).collect();
        let (order, x) = min_inversion_order(&parent, 0, &c, &d);
        assert_eq!(order.len(), n);
        assert_eq!(x, cost(&order, &c, &d));
        let chain: Vec<usize> = (0..n).map(|i| if i == 0 { NONE } else { i - 1 }).collect();
        assert_eq!(min_inversion_order(&chain, 0, &c, &d).0, (0..n).collect::<Vec<_>>());
    }
}
