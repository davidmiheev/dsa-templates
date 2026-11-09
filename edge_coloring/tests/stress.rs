//! Stress test: the colouring is proper and uses exactly the maximum degree many colours.

use edge_coloring::bipartite_edge_coloring;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn check(l: usize, r: usize, edges: &[(usize, usize)], what: &str) {
    let (k, color) = bipartite_edge_coloring(l, r, edges);
    let (mut dl, mut dr) = (vec![0; l], vec![0; r]);
    for &(a, b) in edges {
        dl[a] += 1;
        dr[b] += 1;
    }
    let delta = dl.iter().chain(dr.iter()).copied().max().unwrap_or(0);
    assert_eq!(k, delta, "{what}");
    assert_eq!(color.len(), edges.len());
    let mut seen_l = vec![vec![false; k]; l];
    let mut seen_r = vec![vec![false; k]; r];
    for (i, &(a, b)) in edges.iter().enumerate() {
        let c = color[i];
        assert!(c < k, "{what}: colour out of range");
        assert!(!seen_l[a][c] && !seen_r[b][c], "{what}: edge {i} clashes");
        seen_l[a][c] = true;
        seen_r[b][c] = true;
    }
}

#[test]
fn proper_and_optimal() {
    let mut st = Stress::new("edge_coloring");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(case as u64);
        let (l, r) = (rng.range(1, 7) as usize, rng.range(1, 7) as usize);
        let m = rng.range(0, 30) as usize;
        let edges: Vec<(usize, usize)> = (0..m).map(|_| (rng.below(l), rng.below(r))).collect();
        check(l, r, &edges, &format!("case {case}: {edges:?}"));
        st.case();
    }
    // larger, skewed instances
    for case in 0..40 * scale() {
        let mut rng = Rng::new(case as u64 + 1000);
        let (l, r) = (rng.range(1, 60) as usize, rng.range(1, 60) as usize);
        let m = rng.range(0, 800) as usize;
        let hub = rng.below(l);
        let edges: Vec<(usize, usize)> = (0..m)
            .map(|_| (if rng.chance(1, 4) { hub } else { rng.below(l) }, rng.below(r)))
            .collect();
        check(l, r, &edges, &format!("large case {case}"));
        st.case();
    }
    st.done();
}
