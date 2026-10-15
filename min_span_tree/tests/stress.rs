//! Stress test: Kruskal against an exhaustive search over edge subsets.

use min_span_tree::minimum_spanning_forest;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn components(n: usize, edges: &[(usize, usize, u32)], mask: u32) -> usize {
    let mut p: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        if p[x] != x {
            let r = find(p, p[x]);
            p[x] = r;
        }
        p[x]
    }
    let mut c = n;
    for (i, &(u, v, _)) in edges.iter().enumerate() {
        if mask >> i & 1 == 1 {
            let (a, b) = (find(&mut p, u), find(&mut p, v));
            if a != b {
                p[a] = b;
                c -= 1;
            }
        }
    }
    c
}

#[test]
fn matches_exhaustive_search() {
    let mut st = Stress::new("min_span_tree/naive");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 7) as usize;
        let m = rng.range(0, 10) as usize;
        let edges: Vec<(usize, usize, u32)> = (0..m).map(|_| (rng.below(n), rng.below(n), rng.below(6) as u32)).collect();
        let full = components(n, &edges, (1 << m) - 1);
        // Minimum weight over edge subsets that connect as much as everything does.
        let best = (0..1u32 << m).filter(|&mask| components(n, &edges, mask) == full).map(|mask| (0..m).filter(|&i| mask >> i & 1 == 1).map(|i| edges[i].2).sum::<u32>()).min().unwrap();
        let (total, chosen, comps) = minimum_spanning_forest(n, &edges);
        assert_eq!((total, comps), (best, full), "case {case}");
        assert_eq!(chosen.len(), n - full);
        let mask = chosen.iter().fold(0u32, |a, &e| a | 1 << e);
        assert_eq!(components(n, &edges, mask), full, "chosen edges span, case {case}");
        assert_eq!(chosen.iter().map(|&e| edges[e].2).sum::<u32>(), total);
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn large_input() {
    let mut st = Stress::new("min_span_tree/large");
    let (n, m) = (500_000, 1_000_000);
    let mut rng = Rng::new(2);
    let mut edges: Vec<(usize, usize, u64)> = (1..n).map(|i| (rng.below(i), i, rng.below(1_000_000_000) as u64)).collect();
    edges.extend((0..m - (n - 1)).map(|_| (rng.below(n), rng.below(n), rng.below(1_000_000_000) as u64)));
    st.case();
    let (_, chosen, comps) = minimum_spanning_forest(n, &edges);
    assert_eq!((comps, chosen.len()), (1, n - 1));
    st.ops(m);
    st.done();
}

#[test]
#[should_panic(expected = "outside")]
fn bad_endpoint_panics() {
    minimum_spanning_forest(2, &[(0, 2, 1u8)]);
}
