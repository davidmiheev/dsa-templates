//! Stress test: counts against explicit enumeration of edge subsets / edge permutations.

use modint::ModInt998244353 as Mint;
use spanning_count::{count_arborescences, count_eulerian_circuits, count_spanning_trees};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const P: u32 = 998_244_353;

#[test]
fn undirected_spanning_trees() {
    let mut st = Stress::new("spanning_count/undirected");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 6) as usize;
        let m = rng.range(0, 9) as usize;
        let edges: Vec<(usize, usize)> = (0..m).map(|_| (rng.below(n), rng.below(n))).collect();
        let mut count = 0u64;
        for mask in 0u32..1 << m {
            if mask.count_ones() as usize != n - 1 {
                continue;
            }
            let mut comp: Vec<usize> = (0..n).collect();
            let mut ok = true;
            for (i, &(u, v)) in edges.iter().enumerate() {
                if mask >> i & 1 == 1 {
                    let (a, b) = (comp[u], comp[v]);
                    if a == b {
                        ok = false;
                        break;
                    }
                    for c in comp.iter_mut() {
                        if *c == b {
                            *c = a;
                        }
                    }
                }
            }
            if ok {
                count += 1;
            }
        }
        assert_eq!(count_spanning_trees::<P>(n, &edges), Mint::new(count), "case {case}: {edges:?}");
        st.case();
    }
    st.done();
}

fn arborescences_brute(n: usize, root: usize, edges: &[(usize, usize)]) -> u64 {
    // choose one incoming edge for every non-root vertex; valid iff no cycle
    let mut by_head: Vec<Vec<usize>> = vec![vec![]; n];
    for &(u, v) in edges {
        if u != v {
            by_head[v].push(u);
        }
    }
    let mut total = 0;
    let mut parent = vec![root; n];
    fn rec(v: usize, n: usize, root: usize, by_head: &[Vec<usize>], parent: &mut Vec<usize>, total: &mut u64) {
        if v == n {
            for s in 0..n {
                let (mut x, mut steps) = (s, 0);
                while x != root {
                    x = parent[x];
                    steps += 1;
                    if steps > n {
                        return;
                    }
                }
            }
            *total += 1;
            return;
        }
        if v == root {
            return rec(v + 1, n, root, by_head, parent, total);
        }
        for &u in &by_head[v] {
            parent[v] = u;
            rec(v + 1, n, root, by_head, parent, total);
        }
    }
    rec(0, n, root, &by_head, &mut parent, &mut total);
    total
}

#[test]
fn directed_arborescences() {
    let mut st = Stress::new("spanning_count/directed");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64 + 5);
        let n = rng.range(1, 6) as usize;
        let m = rng.range(0, 10) as usize;
        let edges: Vec<(usize, usize)> = (0..m).map(|_| (rng.below(n), rng.below(n))).collect();
        let root = rng.below(n);
        assert_eq!(count_arborescences::<P>(n, root, &edges), Mint::new(arborescences_brute(n, root, &edges)), "case {case}");
        st.case();
    }
    st.done();
}

fn circuits_brute(n: usize, edges: &[(usize, usize)]) -> u64 {
    // count edge orderings starting with edge 0 that chain head-to-tail and close up
    let m = edges.len();
    let mut used = vec![false; m];
    used[0] = true;
    fn rec(cur: usize, count: usize, m: usize, first_tail: usize, edges: &[(usize, usize)], used: &mut Vec<bool>) -> u64 {
        if count == m {
            return (cur == first_tail) as u64;
        }
        let mut total = 0;
        for e in 0..m {
            if !used[e] && edges[e].0 == cur {
                used[e] = true;
                total += rec(edges[e].1, count + 1, m, first_tail, edges, used);
                used[e] = false;
            }
        }
        total
    }
    let _ = n;
    rec(edges[0].1, 1, m, edges[0].0, edges, &mut used)
}

#[test]
fn eulerian_circuits() {
    let mut st = Stress::new("spanning_count/euler");
    let mut nonzero = 0;
    for case in 0..2000 * scale() {
        let mut rng = Rng::new(case as u64 + 11);
        let n = rng.range(1, 5) as usize;
        // random closed walks glued together make balanced graphs
        let mut edges = Vec::new();
        for _ in 0..rng.range(1, 3) {
            let len = rng.range(1, 4) as usize;
            let vs: Vec<usize> = (0..len).map(|_| rng.below(n)).collect();
            for i in 0..len {
                edges.push((vs[i], vs[(i + 1) % len]));
            }
        }
        if rng.chance(1, 5) {
            edges.push((rng.below(n), rng.below(n))); // often unbalanced
        }
        if edges.len() > 9 {
            continue;
        }
        let expect = circuits_brute(n, &edges);
        nonzero += (expect > 0) as usize;
        assert_eq!(count_eulerian_circuits::<P>(n, &edges), Mint::new(expect), "case {case}: {edges:?}");
        st.case();
    }
    assert!(nonzero > 100);
    st.done();
}
