//! Stress test: min over tree nodes of (cut + weight) equals the minimum over all vertex subsets.

use extreme_sets::{extreme_set_tree, NONE};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn tree_nodes_contain_the_minimum() {
    let mut st = Stress::new("extreme_sets");
    for case in 0..2500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 8) as usize;
        let mut w = vec![vec![0i64; n]; n];
        let density = rng.range(1, 9) as u64;
        let mut edges = Vec::new();
        for u in 0..n {
            for v in u + 1..n {
                if rng.chance(density, 10) {
                    let x = rng.range(0, 6);
                    w[u][v] = x;
                    w[v][u] = x;
                    edges.push((u, v, x));
                }
            }
        }
        let t = extreme_set_tree(n, &edges);
        assert_eq!(t.nodes(), 2 * n - 1);
        assert_eq!(t.parent[t.root()], NONE);
        // leaf sets of every node
        let mut set = vec![0u32; t.nodes()];
        for v in 0..n {
            set[v] = 1 << v;
        }
        for v in n..t.nodes() {
            assert_eq!(t.children[v].len(), 2);
            set[v] = t.children[v].iter().fold(0, |m, &c| {
                assert!(m & set[c] == 0, "children overlap");
                m | set[c]
            });
        }
        assert_eq!(set[t.root()], (1 << n) - 1);
        let cut_of = |mask: u32| -> i64 {
            let mut c = 0;
            for u in 0..n {
                for v in 0..n {
                    if mask >> u & 1 == 1 && mask >> v & 1 == 0 {
                        c += w[u][v];
                    }
                }
            }
            c
        };
        for v in 0..t.nodes() {
            assert_eq!(t.cut[v], cut_of(set[v]), "case {case}: cut of node {v}");
        }
        for _ in 0..4 {
            let a: Vec<i64> = (0..n).map(|_| if rng.chance(1, 4) { 0 } else { rng.range(0, 12) }).collect();
            let value = |mask: u32| cut_of(mask) + (0..n).filter(|&v| mask >> v & 1 == 1).map(|v| a[v]).sum::<i64>();
            let brute = (1u32..1 << n).map(value).min().unwrap();
            let via_tree = (0..t.nodes()).map(|v| value(set[v])).min().unwrap();
            assert_eq!(via_tree, brute, "case {case}: n={n} {edges:?} a={a:?}");
        }
        st.case();
    }
    st.done();
}
