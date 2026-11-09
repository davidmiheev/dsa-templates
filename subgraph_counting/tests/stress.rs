//! Stress test: triangles, cliques and per-edge 4-cycles against brute-force enumeration.

use subgraph_counting::{for_each_clique, for_each_triangle, four_cycles_per_edge};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn random_simple(rng: &mut Rng, n: usize, density: u64) -> Vec<(usize, usize)> {
    let mut edges = Vec::new();
    for u in 0..n {
        for v in u + 1..n {
            if rng.chance(density, 10) {
                edges.push(if rng.chance(1, 2) { (u, v) } else { (v, u) });
            }
        }
    }
    edges
}

#[test]
fn triangles_and_cliques() {
    let mut st = Stress::new("subgraph/cliques");
    for case in 0..2000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 9) as usize;
        let density = rng.range(1, 9) as u64;
        let edges = random_simple(&mut rng, n, density);
        let mut adj = vec![vec![false; n]; n];
        for &(u, v) in &edges {
            adj[u][v] = true;
            adj[v][u] = true;
        }
        let mut tri = Vec::new();
        for_each_triangle(n, &edges, |a, b, c| {
            let mut t = [a, b, c];
            t.sort();
            tri.push(t);
        });
        tri.sort();
        let mut expect = Vec::new();
        let mut cliques_expect = std::collections::BTreeSet::new();
        for mask in 1u32..1 << n {
            let vs: Vec<usize> = (0..n).filter(|&v| mask >> v & 1 == 1).collect();
            if vs.iter().all(|&a| vs.iter().all(|&b| a == b || adj[a][b])) {
                cliques_expect.insert(mask);
                if vs.len() == 3 {
                    expect.push([vs[0], vs[1], vs[2]]);
                }
            }
        }
        expect.sort();
        assert_eq!(tri, expect, "case {case}: {edges:?}");
        let mut got = std::collections::BTreeSet::new();
        let mut total = 0;
        for_each_clique(n, &edges, |c| {
            let mask = c.iter().fold(0u32, |m, &v| m | 1 << v);
            assert!(got.insert(mask), "clique reported twice");
            total += 1;
        });
        assert_eq!(got, cliques_expect, "case {case}: {edges:?}");
        st.case();
        st.ops(total);
    }
    st.done();
}

#[test]
fn four_cycles() {
    let mut st = Stress::new("subgraph/c4");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(case as u64 + 77);
        let n = rng.range(4, 8) as usize;
        let m = rng.range(0, 14) as usize;
        let edges: Vec<(usize, usize)> = (0..m)
            .map(|_| {
                let u = rng.below(n);
                (u, (u + 1 + rng.below(n - 1)) % n)
            })
            .collect();
        let got = four_cycles_per_edge(n, &edges);
        // Brute force over 4-subsets of edges.
        let mut expect = vec![0u64; m];
        for a in 0..m {
            for b in a + 1..m {
                for c in b + 1..m {
                    for d in c + 1..m {
                        let es = [edges[a], edges[b], edges[c], edges[d]];
                        let mut deg = std::collections::HashMap::new();
                        for &(u, v) in &es {
                            *deg.entry(u).or_insert(0) += 1;
                            *deg.entry(v).or_insert(0) += 1;
                        }
                        if deg.len() == 4 && deg.values().all(|&x| x == 2) {
                            // 2-regular on 4 vertices: a 4-cycle, or two double edges (not a cycle)
                            let mut seen = vec![es[0].0];
                            let mut cur = es[0].0;
                            let mut used = [false; 4];
                            let mut ok = true;
                            for _ in 0..4 {
                                let mut found = false;
                                for (i, &(u, v)) in es.iter().enumerate() {
                                    if !used[i] && (u == cur || v == cur) {
                                        used[i] = true;
                                        cur = if u == cur { v } else { u };
                                        found = true;
                                        break;
                                    }
                                }
                                if !found {
                                    ok = false;
                                    break;
                                }
                                seen.push(cur);
                            }
                            // a closed walk visiting 4 distinct vertices
                            let closed = ok && seen[4] == seen[0] && {
                                let mut distinct = seen[..4].to_vec();
                                distinct.sort();
                                distinct.dedup();
                                distinct.len() == 4
                            };
                            if closed {
                                for &i in &[a, b, c, d] {
                                    expect[i] += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(got, expect, "case {case}: n={n} {edges:?}");
        st.case();
        st.ops(m);
    }
    st.done();
}
