//! Stress test: `LinkCutTree` against a naive forest (adjacency sets + BFS), using the
//! non-commutative affine-composition monoid so that path orientation is checked.

use algebra::Monoid;
use link_cut_tree::LinkCutTree;
use std::collections::BTreeSet;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const P: u64 = 998_244_353;

/// `op(f, g)` = apply `f` first, then `g`.
enum Compose {}
impl Monoid for Compose {
    type S = (u64, u64);
    fn identity() -> (u64, u64) {
        (1, 0)
    }
    fn op(f: &(u64, u64), g: &(u64, u64)) -> (u64, u64) {
        (f.0 * g.0 % P, (g.0 * f.1 + g.1) % P)
    }
}

/// The path from `u` to `v` in the forest `adj`, if connected.
fn naive_path(adj: &[BTreeSet<usize>], u: usize, v: usize) -> Option<Vec<usize>> {
    let n = adj.len();
    let mut parent = vec![usize::MAX; n];
    let mut seen = vec![false; n];
    let mut queue = vec![u];
    seen[u] = true;
    let mut i = 0;
    while i < queue.len() {
        let x = queue[i];
        i += 1;
        for &y in &adj[x] {
            if !seen[y] {
                seen[y] = true;
                parent[y] = x;
                queue.push(y);
            }
        }
    }
    if !seen[v] {
        return None;
    }
    let mut path = vec![v];
    while *path.last().unwrap() != u {
        path.push(parent[*path.last().unwrap()]);
    }
    path.reverse();
    Some(path)
}

#[test]
fn matches_naive_forest() {
    let mut st = Stress::new("link_cut_tree/naive");
    for case in 0..800 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 14) as usize;
        let mut val: Vec<(u64, u64)> = (0..n).map(|_| (1 + rng.below(P as usize - 1) as u64, rng.below(P as usize) as u64)).collect();
        let mut t = LinkCutTree::<Compose>::new(&val);
        let mut adj = vec![BTreeSet::new(); n];
        st.case();
        for step in 0..120 {
            let ctx = format!("case {case} step {step}");
            let (u, v) = (rng.below(n), rng.below(n));
            match rng.below(6) {
                0 | 1 => {
                    let connected = naive_path(&adj, u, v).is_some();
                    if !connected {
                        t.link(u, v);
                        adj[u].insert(v);
                        adj[v].insert(u);
                    } else {
                        assert!(catch_unwind(AssertUnwindSafe(|| t.link(u, v))).is_err(), "link of connected vertices must panic {ctx}");
                    }
                }
                2 => {
                    if adj[u].contains(&v) {
                        t.cut(u, v);
                        adj[u].remove(&v);
                        adj[v].remove(&u);
                    } else {
                        assert!(catch_unwind(AssertUnwindSafe(|| t.cut(u, v))).is_err(), "cut of a missing edge must panic {ctx}");
                    }
                }
                3 => {
                    let f = (1 + rng.below(P as usize - 1) as u64, rng.below(P as usize) as u64);
                    t.set(u, f);
                    val[u] = f;
                    assert_eq!(t.get(u), f, "get {ctx}");
                }
                4 => {
                    assert_eq!(t.connected(u, v), naive_path(&adj, u, v).is_some(), "connected({u}, {v}) {ctx}");
                    if let Some(path) = naive_path(&adj, u, v) {
                        let want = path.iter().fold((1, 0), |acc, &x| Compose::op(&acc, &val[x]));
                        assert_eq!(t.path_prod(u, v), want, "path_prod({u}, {v}) {ctx}");
                    } else {
                        assert!(catch_unwind(AssertUnwindSafe(|| t.path_prod(u, v))).is_err());
                    }
                }
                _ => {
                    let r = rng.below(n);
                    if let (Some(pu), Some(pv)) = (naive_path(&adj, r, u), naive_path(&adj, r, v)) {
                        let common = pu.iter().zip(&pv).take_while(|(a, b)| a == b).count();
                        assert_eq!(t.lca(r, u, v), pu[common - 1], "lca(root {r}, {u}, {v}) {ctx}");
                    }
                }
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn large_path_and_random_trees() {
    let mut st = Stress::new("link_cut_tree/large");
    let n = 200_000;
    let mut rng = Rng::new(8);
    let val: Vec<(u64, u64)> = (0..n).map(|_| (1 + rng.below(P as usize - 1) as u64, rng.below(P as usize) as u64)).collect();
    let mut t = LinkCutTree::<Compose>::new(&val);
    st.case();
    for i in 1..n {
        t.link(i - 1, i); // a path: the worst case for naive splay depth
    }
    let want = val.iter().fold((1, 0), |acc, f| Compose::op(&acc, f));
    assert_eq!(t.path_prod(0, n - 1), want);
    for k in 0..100_000 * scale() {
        let (u, v) = (rng.below(n), rng.below(n));
        let _ = t.path_prod(u, v);
        if k % 3 == 0 {
            t.set(u, (1, 1));
        }
        st.ops(1);
    }
    st.done();
}

#[test]
#[should_panic(expected = "out of range")]
fn vertex_out_of_range_panics() {
    LinkCutTree::<Compose>::new(&[(1, 0)]).set(1, (1, 0));
}
