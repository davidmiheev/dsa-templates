//! Stress test: the diameter against the best over all spanning trees (enumerating subsets of n-1 edges).

use graph_center::minimum_diameter_spanning_tree;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn tree_diameter(n: usize, tree: &[(usize, usize, u64)]) -> u64 {
    // all-pairs via Floyd on the tree
    let inf = u64::MAX / 4;
    let mut d = vec![vec![inf; n]; n];
    for i in 0..n {
        d[i][i] = 0;
    }
    for &(u, v, w) in tree {
        d[u][v] = d[u][v].min(w);
        d[v][u] = d[v][u].min(w);
    }
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                d[i][j] = d[i][j].min(d[i][k] + d[k][j]);
            }
        }
    }
    let mut best = 0;
    for i in 0..n {
        for j in 0..n {
            assert!(d[i][j] < inf, "not connected");
            best = best.max(d[i][j]);
        }
    }
    best
}

#[test]
fn matches_exhaustive_spanning_trees() {
    let mut st = Stress::new("graph_center");
    for case in 0..2500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 6) as usize;
        let mut edges: Vec<(usize, usize, u64)> = (1..n).map(|v| (rng.below(v), v, rng.range(0, 9) as u64)).collect();
        for _ in 0..rng.range(0, 5) {
            edges.push((rng.below(n), rng.below(n), rng.range(0, 9) as u64));
        }
        let m = edges.len();
        let mut best = u64::MAX;
        for mask in 0u32..1 << m {
            if mask.count_ones() as usize != n - 1 {
                continue;
            }
            let chosen: Vec<(usize, usize, u64)> = (0..m).filter(|&i| mask >> i & 1 == 1).map(|i| edges[i]).collect();
            let mut comp: Vec<usize> = (0..n).collect();
            let mut ok = true;
            for &(u, v, _) in &chosen {
                if comp[u] == comp[v] {
                    ok = false;
                    break;
                }
                let (a, b) = (comp[u], comp[v]);
                for c in comp.iter_mut() {
                    if *c == b {
                        *c = a;
                    }
                }
            }
            if ok {
                best = best.min(tree_diameter(n, &chosen));
            }
        }
        let (d, tree) = minimum_diameter_spanning_tree(n, &edges);
        assert_eq!(d, best, "case {case}: n={n} {edges:?}");
        assert_eq!(tree.len(), n - 1);
        let chosen: Vec<(usize, usize, u64)> = tree.iter().map(|&i| edges[i]).collect();
        assert_eq!(tree_diameter(n, &chosen), d, "case {case}: reported diameter differs from the tree's");
        st.case();
    }
    st.done();
}
