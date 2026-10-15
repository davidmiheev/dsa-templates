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
