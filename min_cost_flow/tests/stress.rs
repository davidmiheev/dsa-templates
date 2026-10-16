//! Stress test: exhaustive search over all integral flows on tiny instances; every returned
//! solution is also checked against the primal and dual conditions.

use min_cost_flow::min_cost_b_flow;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn brute(n: usize, b: &[i64], edges: &[(usize, usize, i64, i64, i64)]) -> Option<i128> {
    let m = edges.len();
    let mut best: Option<i128> = None;
    let mut f = vec![0i64; m];
    fn rec(i: usize, n: usize, b: &[i64], edges: &[(usize, usize, i64, i64, i64)], f: &mut Vec<i64>, best: &mut Option<i128>) {
        if i == edges.len() {
            let mut net = vec![0i64; n];
            for (j, &(s, t, _, _, _)) in edges.iter().enumerate() {
                net[s] += f[j];
                net[t] -= f[j];
            }
            if net == b {
                let c: i128 = edges.iter().enumerate().map(|(j, e)| e.4 as i128 * f[j] as i128).sum();
                if best.map_or(true, |x| c < x) {
                    *best = Some(c);
                }
            }
            return;
        }
        for x in edges[i].2..=edges[i].3 {
            f[i] = x;
            rec(i + 1, n, b, edges, f, best);
        }
    }
    rec(0, n, b, edges, &mut f, &mut best);
    best
}

#[test]
fn matches_exhaustive_search() {
    let mut st = Stress::new("min_cost_flow");
    let (mut feasible, mut infeasible) = (0, 0);
    for case in 0..4000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 4) as usize;
        let m = rng.range(0, 6) as usize;
        let edges: Vec<(usize, usize, i64, i64, i64)> = (0..m)
            .map(|_| {
                let l = rng.range(-1, 2);
                (rng.below(n), rng.below(n), l, l + rng.range(0, 2), rng.range(-3, 3))
            })
            .collect();
        let mut b: Vec<i64> = (0..n).map(|_| rng.range(-2, 2)).collect();
        if rng.chance(2, 3) {
            // make the instance balanced
            let s: i64 = b.iter().sum();
            b[0] -= s;
        }
        let expect = brute(n, &b, &edges);
        let r = min_cost_b_flow(n, &b, &edges);
        match expect {
            None => {
                infeasible += 1;
                assert!(!r.feasible, "case {case}: should be infeasible: {b:?} {edges:?}");
            }
            Some(c) => {
                feasible += 1;
                assert!(r.feasible, "case {case}: should be feasible: {b:?} {edges:?}");
                assert_eq!(r.cost, c, "case {case}: {b:?} {edges:?}");
                // primal
                let mut net = vec![0i64; n];
                for (j, &(s, t, l, u, _)) in edges.iter().enumerate() {
                    assert!(l <= r.flow[j] && r.flow[j] <= u);
                    net[s] += r.flow[j];
                    net[t] -= r.flow[j];
                }
                assert_eq!(net, b);
                // dual
                for (j, &(s, t, l, u, cst)) in edges.iter().enumerate() {
                    let rc = cst + r.potential[s] - r.potential[t];
                    if r.flow[j] > l {
                        assert!(rc <= 0, "case {case}: slackness (down) on edge {j}");
                    }
                    if r.flow[j] < u {
                        assert!(rc >= 0, "case {case}: slackness (up) on edge {j}");
                    }
                }
            }
        }
        st.case();
    }
    assert!(feasible > 300 && infeasible > 300, "{feasible}/{infeasible}");
    st.done();
}
