//! Stress test: per-edge merge times against recomputing the SCCs after every insertion.

use incremental_scc::{incremental_scc_times, NEVER};
use scc::Scc;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_recomputation() {
    let mut st = Stress::new("incremental_scc");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 8) as usize;
        let m = rng.range(0, 20) as usize;
        let edges: Vec<(usize, usize)> = (0..m).map(|_| (rng.below(n), rng.below(n))).collect();
        let times = incremental_scc_times(n, &edges);
        let mut expect = vec![NEVER; m];
        for t in 0..m {
            let s = Scc::from_edges(n, &edges[..=t]);
            for e in 0..=t {
                if expect[e] == NEVER && s.comp[edges[e].0] == s.comp[edges[e].1] {
                    expect[e] = t;
                }
            }
        }
        assert_eq!(times, expect, "case {case}: n={n} {edges:?}");
        st.case();
        st.ops(m);
    }
    st.done();
}
