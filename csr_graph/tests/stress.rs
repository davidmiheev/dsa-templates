//! Stress test: `Csr` against adjacency lists built the obvious way.

use csr_graph::Csr;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_adjacency_lists() {
    let mut st = Stress::new("csr_graph/naive");
    for case in 0..1000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 12) as usize;
        let edges: Vec<(usize, usize)> = (0..rng.range(0, 30)).map(|_| (rng.below(n), rng.below(n))).collect();
        let mut out = vec![vec![]; n];
        let mut both = vec![vec![]; n];
        let mut inn = vec![vec![]; n];
        for (e, &(u, v)) in edges.iter().enumerate() {
            out[u].push((v, e));
            inn[v].push((u, e));
            both[u].push((v, e));
            both[v].push((u, e));
        }
        let d = Csr::directed(n, &edges);
        let u = Csr::undirected(n, &edges);
        let r = d.reversed();
        for v in 0..n {
            let got: Vec<(usize, usize)> = d.arcs(v).map(|i| (d.to[i], d.id[i])).collect();
            assert_eq!(got, out[v], "directed case {case} v={v}");
            let mut got: Vec<(usize, usize)> = u.arcs(v).map(|i| (u.to[i], u.id[i])).collect();
            let mut want = both[v].clone();
            got.sort_unstable();
            want.sort_unstable();
            assert_eq!(got, want, "undirected case {case} v={v}");
            let mut got: Vec<(usize, usize)> = r.arcs(v).map(|i| (r.to[i], r.id[i])).collect();
            let mut want = inn[v].clone();
            got.sort_unstable();
            want.sort_unstable();
            assert_eq!(got, want, "reversed case {case} v={v}");
            assert_eq!(d.degree(v), out[v].len());
            assert_eq!(d.neighbors(v).len(), out[v].len());
        }
        assert_eq!((d.len(), d.num_arcs(), u.num_arcs()), (n, edges.len(), 2 * edges.len()));
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
#[should_panic(expected = "endpoint outside")]
fn out_of_range_endpoint_panics() {
    Csr::directed(2, &[(0, 2)]);
}
