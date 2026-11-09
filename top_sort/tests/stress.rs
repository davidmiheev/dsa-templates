//! Stress test: `topological_sort` returns a valid order exactly when the graph is acyclic
//! (acyclicity decided independently by a colouring DFS).

use top_sort::topological_sort;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

/// Iterative three-colour DFS cycle detection.
fn has_cycle(g: &[Vec<usize>]) -> bool {
    let mut color = vec![0u8; g.len()]; // 0 new, 1 on stack, 2 done
    for s in 0..g.len() {
        if color[s] != 0 {
            continue;
        }
        let mut stack = vec![(s, 0usize)];
        color[s] = 1;
        while let Some(&mut (v, ref mut i)) = stack.last_mut() {
            if *i < g[v].len() {
                let u = g[v][*i];
                *i += 1;
                match color[u] {
                    0 => {
                        color[u] = 1;
                        stack.push((u, 0));
                    }
                    1 => return true,
                    _ => {}
                }
            } else {
                color[v] = 2;
                stack.pop();
            }
        }
    }
    false
}

fn is_valid_order(g: &[Vec<usize>], order: &[usize]) -> bool {
    let n = g.len();
    let mut pos = vec![usize::MAX; n];
    for (i, &v) in order.iter().enumerate() {
        if v >= n || pos[v] != usize::MAX {
            return false; // out of range or repeated
        }
        pos[v] = i;
    }
    order.len() == n && (0..n).all(|u| g[u].iter().all(|&v| pos[u] < pos[v]))
}

#[test]
fn valid_order_iff_acyclic() {
    let mut st = Stress::new("top_sort/naive");
    let (mut dags, mut cyclic) = (0, 0);
    for case in 0..2000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 25) as usize;
        let mut g = vec![Vec::new(); n];
        // Hidden random ranking: edges going "forward" in it form a DAG.
        let mut rank: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            rank.swap(i, rng.below(i + 1));
        }
        for _ in 0..rng.below(3 * n + 1) {
            let (a, b) = (rng.below(n), rng.below(n));
            match case % 3 {
                0 => {
                    // acyclic by construction (parallel edges allowed, no self-loops)
                    if rank[a] < rank[b] {
                        g[a].push(b)
                    } else if rank[b] < rank[a] {
                        g[b].push(a)
                    }
                }
                _ => g[a].push(b), // arbitrary: loops, cycles, parallel edges
            }
        }
        let expect_cycle = has_cycle(&g);
        if expect_cycle { cyclic += 1 } else { dags += 1 }
        match topological_sort(&g) {
            Some(order) => {
                assert!(!expect_cycle, "returned an order for a cyclic graph, case {case}: {g:?}");
                assert!(is_valid_order(&g, &order), "invalid order {order:?}, case {case}: {g:?}");
            }
            None => assert!(expect_cycle, "returned None for a DAG, case {case}: {g:?}"),
        }
        st.case();
        st.ops(1);
    }
    assert!(dags > 300 && cyclic > 300, "generator must produce both outcomes ({dags} DAGs, {cyclic} cyclic)");
    st.done();
}

#[test]
fn large_dag_and_single_back_edge() {
    let mut st = Stress::new("top_sort/large");
    let n = 200_000;
    let mut rng = Rng::new(4);
    let mut g = vec![Vec::new(); n];
    for _ in 0..600_000 {
        let (a, b) = (rng.below(n), rng.below(n));
        if a < b {
            g[a].push(b);
        } else if b < a {
            g[b].push(a);
        }
    }
    st.case();
    let order = topological_sort(&g).expect("edges only go from smaller to larger ids");
    assert!(is_valid_order(&g, &order));
    // Closing a cycle with one back edge must be detected.
    g[n - 1].push(0);
    g[0].push(1);
    g[1].push(n - 1);
    assert!(topological_sort(&g).is_none());
    st.ops(2 * n);
    st.done();
}
