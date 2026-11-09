//! Stress test: 2-SAT against brute-force enumeration of all assignments.

use two_sat::TwoSat;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_brute_force() {
    let mut st = Stress::new("two_sat");
    for case in 0..6000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = 1 + rng.below(8);
        let m = rng.below(3 * n + 1);
        let clauses: Vec<(usize, bool, usize, bool)> = (0..m).map(|_| (rng.below(n), rng.chance(1, 2), rng.below(n), rng.chance(1, 2))).collect();
        let mut sat = TwoSat::new(n);
        for &(i, f, j, g) in &clauses {
            sat.add_clause(i, f, j, g);
        }
        let holds = |x: &[bool]| clauses.iter().all(|&(i, f, j, g)| x[i] == f || x[j] == g);
        let exists = (0u32..1 << n).any(|mask| holds(&(0..n).map(|i| mask >> i & 1 == 1).collect::<Vec<_>>()));
        let (ok, x) = sat.solve();
        assert_eq!(ok, exists, "clauses {clauses:?}");
        if ok {
            assert!(holds(&x), "the assignment must satisfy every clause: {clauses:?} {x:?}");
        } else {
            assert!(x.is_empty());
        }
        st.case();
        st.ops(n + m);
    }
    st.done();
}
