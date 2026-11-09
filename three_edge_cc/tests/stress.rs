//! Stress test: 3-edge-connected components against "remove any <= 2 edges and compare reachability".

use three_edge_cc::three_edge_connected_components;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn components_without(n: usize, edges: &[(usize, usize)], skip: [usize; 2]) -> Vec<usize> {
    let mut label: Vec<usize> = (0..n).collect();
    fn find(l: &mut Vec<usize>, x: usize) -> usize {
        if l[x] != x {
            let r = find(l, l[x]);
            l[x] = r;
        }
        l[x]
    }
    for (i, &(u, v)) in edges.iter().enumerate() {
        if i != skip[0] && i != skip[1] {
            let (a, b) = (find(&mut label, u), find(&mut label, v));
            label[a] = b;
        }
    }
    (0..n).map(|v| find(&mut label, v)).collect()
}

#[test]
fn matches_cut_enumeration() {
    let mut st = Stress::new("three_edge_cc");
    for case in 0..4000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 8) as usize;
        let m = rng.range(0, 14) as usize;
        let edges: Vec<(usize, usize)> = (0..m).map(|_| (rng.below(n), rng.below(n))).collect();
        // sep[u][v]: some set of <= 2 edges separates u and v.
        let mut sep = vec![vec![false; n]; n];
        let none = usize::MAX;
        let mut options = vec![[none, none]];
        for a in 0..m {
            options.push([a, none]);
            for b in a + 1..m {
                options.push([a, b]);
            }
        }
        for skip in options {
            let l = components_without(n, &edges, skip);
            for u in 0..n {
                for v in 0..n {
                    if l[u] != l[v] {
                        sep[u][v] = true;
                    }
                }
            }
        }
        let (count, comp) = three_edge_connected_components(n, &edges);
        for u in 0..n {
            for v in 0..n {
                assert_eq!(comp[u] == comp[v], !sep[u][v], "case {case}: {edges:?} u={u} v={v} comp={comp:?}");
            }
        }
        assert_eq!(count, comp.iter().max().unwrap() + 1);
        st.case();
        st.ops(n * n);
    }
    st.done();
}
