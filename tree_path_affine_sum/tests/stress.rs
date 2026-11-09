//! Stress test: `TreePathAffineSum` against recomputing every path sum by DFS.

use tree_path_affine_sum::TreePathAffineSum;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const P: u64 = 998_244_353;

fn brute(n: usize, a: &[u64], edges: &[(usize, usize, u64, u64)], r: usize) -> u64 {
    let mut adj = vec![vec![]; n];
    for (e, &(u, v, _, _)) in edges.iter().enumerate() {
        adj[u].push((v, e));
        adj[v].push((u, e));
    }
    // Iterative DFS computing P(r, y) = f_{e0}(P(next, y)) top-down: process in reverse BFS.
    let mut parent = vec![usize::MAX; n];
    let mut pedge = vec![usize::MAX; n];
    let mut order = vec![r];
    let mut i = 0;
    while i < order.len() {
        let v = order[i];
        i += 1;
        for &(u, e) in &adj[v] {
            if u != parent[v] {
                parent[u] = v;
                pedge[u] = e;
                order.push(u);
            }
        }
    }
    // (sum, count) of the subtree, mapped across edges.
    let mut s = vec![0u64; n];
    let mut cnt = vec![1u64; n];
    for &v in order.iter().rev() {
        s[v] = (s[v] + a[v]) % P;
        if v != r {
            let (b, c) = (edges[pedge[v]].2, edges[pedge[v]].3);
            let p = parent[v];
            s[p] = (s[p] + b * s[v] + c * (cnt[v] % P)) % P;
            cnt[p] += cnt[v];
        }
    }
    s[r]
}

#[test]
fn matches_brute_force() {
    let mut st = Stress::new("tree_path_affine_sum/naive");
    for case in 0..800 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 18) as usize;
        let shape = case % 4;
        let mut edges: Vec<(usize, usize, u64, u64)> = Vec::new();
        for i in 1..n {
            let p = match shape {
                0 => rng.below(i),
                1 => i - 1,
                2 => 0,
                _ => (i - 1) / 2,
            };
            let (u, v) = if rng.chance(1, 2) { (p, i) } else { (i, p) };
            edges.push((u, v, 1 + rng.below(P as usize - 1) as u64, rng.below(P as usize) as u64));
        }
        let mut a: Vec<u64> = (0..n).map(|_| rng.below(P as usize) as u64).collect();
        let root = rng.below(n);
        let mut t = TreePathAffineSum::<P>::new(&a, &edges, root);
        st.case();
        for step in 0..60 {
            let ctx = format!("case {case} step {step} n={n} root={root}");
            if n > 1 && rng.chance(1, 2) {
                let e = rng.below(n - 1);
                let (b, c) = (1 + rng.below(P as usize - 1) as u64, rng.below(P as usize) as u64);
                edges[e].2 = b;
                edges[e].3 = c;
                t.set_edge(e, b, c);
            } else {
                let (v, x) = (rng.below(n), rng.below(P as usize) as u64);
                a[v] = x;
                t.set_vertex(v, x);
            }
            assert_eq!(t.sum_from_root(), brute(n, &a, &edges, root), "sum_from_root {ctx}");
            let r = rng.below(n);
            assert_eq!(t.sum_from(r), brute(n, &a, &edges, r), "sum_from({r}) {ctx}");
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn large_random_and_path_trees() {
    let mut st = Stress::new("tree_path_affine_sum/large");
    for shape in 0..2 {
        let n = 200_000;
        let mut rng = Rng::new(shape + 3);
        let edges: Vec<(usize, usize, u64, u64)> = (1..n)
            .map(|i| (if shape == 0 { rng.below(i) } else { i - 1 }, i, 1 + rng.below(P as usize - 1) as u64, rng.below(P as usize) as u64))
            .collect();
        let a: Vec<u64> = (0..n).map(|_| rng.below(P as usize) as u64).collect();
        let mut t = TreePathAffineSum::<P>::new(&a, &edges, 0);
        st.case();
        let mut acc = 0u64;
        for _ in 0..50_000 * scale() {
            if rng.chance(1, 2) {
                t.set_vertex(rng.below(n), rng.below(P as usize) as u64);
            } else {
                t.set_edge(rng.below(n - 1), 1 + rng.below(P as usize - 1) as u64, rng.below(P as usize) as u64);
            }
            acc ^= t.sum_from(rng.below(n));
            st.ops(1);
        }
        // Rooting at the construction root agrees with the fixed-root query.
        assert_eq!(t.sum_from(0), t.sum_from_root());
        let _ = acc;
    }
    st.done();
}
