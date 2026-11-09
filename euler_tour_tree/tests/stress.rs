//! Stress test: `EulerTourTree` against a naive forest (adjacency sets + BFS components) with a
//! (sum, count) monoid and additive action.

use algebra::{MapMonoid, Monoid};
use euler_tour_tree::EulerTourTree;
use std::collections::BTreeSet;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

enum SumLen {}
impl Monoid for SumLen {
    type S = (i64, i64);
    fn identity() -> (i64, i64) {
        (0, 0)
    }
    fn op(a: &(i64, i64), b: &(i64, i64)) -> (i64, i64) {
        (a.0 + b.0, a.1 + b.1)
    }
}
enum Add {}
impl MapMonoid for Add {
    type M = SumLen;
    type F = i64;
    fn id_map() -> i64 {
        0
    }
    fn apply(f: &i64, x: &(i64, i64)) -> (i64, i64) {
        (x.0 + f * x.1, x.1)
    }
    fn compose(f: &i64, g: &i64) -> i64 {
        f + g
    }
}

fn component(adj: &[BTreeSet<usize>], s: usize, banned: Option<(usize, usize)>) -> Vec<usize> {
    let mut seen = vec![false; adj.len()];
    seen[s] = true;
    let mut q = vec![s];
    let mut i = 0;
    while i < q.len() {
        let x = q[i];
        i += 1;
        for &y in &adj[x] {
            if banned == Some((x, y)) || banned == Some((y, x)) {
                continue;
            }
            if !seen[y] {
                seen[y] = true;
                q.push(y);
            }
        }
    }
    q
}

#[test]
fn matches_naive_forest() {
    let mut st = Stress::new("euler_tour_tree/naive");
    for case in 0..800 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 14) as usize;
        let mut val: Vec<i64> = (0..n).map(|_| rng.range(-50, 50)).collect();
        let init: Vec<(i64, i64)> = val.iter().map(|&x| (x, 1)).collect();
        let mut t = EulerTourTree::<Add>::new(&init);
        let mut adj = vec![BTreeSet::new(); n];
        st.case();
        for step in 0..120 {
            let ctx = format!("case {case} step {step}");
            let (u, v) = (rng.below(n), rng.below(n));
            match rng.below(8) {
                0 | 1 => {
                    let same = component(&adj, u, None).contains(&v);
                    if !same {
                        t.link(u, v);
                        adj[u].insert(v);
                        adj[v].insert(u);
                    } else {
                        assert!(catch_unwind(AssertUnwindSafe(|| t.link(u, v))).is_err(), "{ctx}");
                    }
                }
                2 => {
                    if adj[u].contains(&v) {
                        t.cut(u, v);
                        adj[u].remove(&v);
                        adj[v].remove(&u);
                    } else {
                        assert!(catch_unwind(AssertUnwindSafe(|| t.cut(u, v))).is_err(), "{ctx}");
                    }
                }
                3 => {
                    let x = rng.range(-50, 50);
                    t.set(u, (x, 1));
                    val[u] = x;
                }
                4 => {
                    let x = rng.range(-50, 50);
                    t.apply(u, x);
                    val[u] += x;
                    assert_eq!(t.get(u), (val[u], 1), "get {ctx}");
                }
                5 => {
                    let x = rng.range(-50, 50);
                    t.apply_component(u, x);
                    for w in component(&adj, u, None) {
                        val[w] += x;
                    }
                }
                6 if adj[u].contains(&v) => {
                    // subtree of u with parent v
                    let comp = component(&adj, u, Some((u, v)));
                    let want = (comp.iter().map(|&w| val[w]).sum::<i64>(), comp.len() as i64);
                    assert_eq!(t.subtree_prod(u, v), want, "subtree_prod({u}, {v}) {ctx}");
                    let x = rng.range(-50, 50);
                    t.apply_subtree(u, v, x);
                    for w in comp {
                        val[w] += x;
                    }
                }
                _ => {
                    let comp = component(&adj, u, None);
                    let want = (comp.iter().map(|&w| val[w]).sum::<i64>(), comp.len() as i64);
                    assert_eq!(t.component_prod(u), want, "component_prod({u}) {ctx}");
                    assert_eq!(t.connected(u, v), comp.contains(&v), "connected {ctx}");
                }
            }
            st.ops(1);
        }
        for v in 0..n {
            assert_eq!(t.get(v), (val[v], 1), "final value of {v}, case {case}");
        }
    }
    st.done();
}

#[test]
fn large_path_and_churn() {
    let mut st = Stress::new("euler_tour_tree/large");
    let n = 200_000;
    let mut rng = Rng::new(4);
    let init: Vec<(i64, i64)> = (0..n).map(|_| (rng.range(0, 1000), 1)).collect();
    let total: i64 = init.iter().map(|x| x.0).sum();
    let mut t = EulerTourTree::<Add>::new(&init);
    st.case();
    for i in 1..n {
        t.link(i - 1, i);
    }
    assert_eq!(t.component_prod(0), (total, n as i64));
    for _ in 0..50_000 * scale() {
        let i = 1 + rng.below(n - 1);
        // Cut a path edge, look at both halves, and glue it back.
        let (a, b) = (t.subtree_prod(i, i - 1), t.subtree_prod(i - 1, i));
        assert_eq!(a.0 + b.0, t.component_prod(0).0);
        assert_eq!(a.1 + b.1, n as i64);
        st.ops(1);
    }
    st.done();
}

#[test]
#[should_panic(expected = "no such edge")]
fn cut_missing_edge_panics() {
    EulerTourTree::<Add>::new(&[(0, 1), (0, 1)]).cut(0, 1);
}
