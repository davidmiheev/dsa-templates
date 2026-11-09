//! Stress test: connectivity and component sums at every time against recomputation from scratch.

use offline_connectivity::{offline_dynamic_connectivity, EdgeIntervals};
use union_find::{RollbackSumUnionFind, RollbackUnionFind};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_recomputation() {
    let mut st = Stress::new("offline_connectivity");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(2, 8) as usize;
        let times = rng.range(1, 30) as usize;
        // toggle events: at each time maybe toggle one edge
        let mut ev = EdgeIntervals::new();
        let mut alive: Vec<Vec<(usize, usize)>> = Vec::new();
        let mut cur: Vec<(usize, usize)> = Vec::new();
        for t in 0..times {
            if rng.chance(2, 3) {
                let (u, v) = (rng.below(n), rng.below(n));
                if u != v {
                    ev.toggle(t, u, v);
                    let key = (u.min(v), u.max(v));
                    if let Some(p) = cur.iter().position(|&e| e == key) {
                        cur.remove(p);
                    } else {
                        cur.push(key);
                    }
                }
            }
            alive.push(cur.clone());
        }
        let lifetimes = ev.finish(times);
        let values: Vec<i64> = (0..n).map(|_| rng.range(0, 5)).collect();
        // plain connectivity
        let mut uf = RollbackUnionFind::new(n);
        let mut seen = 0;
        offline_dynamic_connectivity(times, &lifetimes, &mut uf, |t, uf| {
            seen += 1;
            let mut comp: Vec<usize> = (0..n).collect();
            for &(a, b) in &alive[t] {
                let (x, y) = (comp[a], comp[b]);
                for c in comp.iter_mut() {
                    if *c == y {
                        *c = x;
                    }
                }
            }
            for a in 0..n {
                for b in 0..n {
                    assert_eq!(uf.same(a, b), comp[a] == comp[b], "case {case} t={t}");
                }
            }
        });
        assert_eq!(seen, times);
        // component sums with point updates inside the walk
        let mut sums = RollbackSumUnionFind::new(&values);
        let mut cur_values = values.clone();
        offline_dynamic_connectivity(times, &lifetimes, &mut sums, |t, sums| {
            let (v, x) = (t % n, (t % 3) as i64);
            sums.add(v, x);
            cur_values[v] += x;
            let mut comp: Vec<usize> = (0..n).collect();
            for &(a, b) in &alive[t] {
                let (p, q) = (comp[a], comp[b]);
                for c in comp.iter_mut() {
                    if *c == q {
                        *c = p;
                    }
                }
            }
            for a in 0..n {
                let expect: i64 = (0..n).filter(|&w| comp[w] == comp[a]).map(|w| cur_values[w]).sum();
                assert_eq!(sums.sum(a), expect, "case {case} t={t}");
            }
        });
        st.case();
        st.ops(times);
    }
    st.done();
}
