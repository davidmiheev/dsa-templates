//! Minimum-cost **b-flow** with lower and upper bounds, negative costs, and optimal dual
//! potentials, by capacity-scaling successive shortest paths (Ahuja–Magnanti–Orlin).
//!
//! A b-flow is a flow `f` with `l_e <= f_e <= u_e` whose net outflow at each vertex `v` equals
//! `b_v` (supplies positive, demands negative). The returned potentials `p` certify optimality:
//! `f_e > l_e => c_e + p_s - p_t <= 0` and `f_e < u_e => c_e + p_s - p_t >= 0`.
//!
//! Method: shift out the lower bounds, saturate the negative-cost arcs so that all reduced costs
//! of residual arcs are non-negative, then for `delta = 2^k, ..., 2, 1` first saturate the arcs
//! of negative reduced cost with residual `>= delta` and push `delta` units along Dijkstra
//! shortest paths (w.r.t. reduced costs, potentials updated with the capped distances) from a
//! vertex of excess `>= delta` to one of deficit `>= delta` in the `delta`-residual graph. This
//! needs `O(m log U)` augmentations, each one `O(m log n)`; it is polynomial even on inputs that
//! defeat plain successive shortest paths. At `delta = 1` an excess that cannot reach a deficit
//! proves infeasibility.
//!
//! # Example
//!
//! ```
//! use min_cost_flow::min_cost_b_flow;
//!
//! // Send 2 units from vertex 0 to vertex 2: through 1 costs 1+1 per unit (capacity 1),
//! // the direct edge costs 5 per unit.
//! let r = min_cost_b_flow(3, &[2, 0, -2], &[(0, 1, 0, 1, 1), (1, 2, 0, 1, 1), (0, 2, 0, 5, 5)]);
//! assert!(r.feasible);
//! assert_eq!(r.cost, 2 + 5);
//! assert_eq!(r.flow, vec![1, 1, 1]);
//! ```

use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// Result of [`min_cost_b_flow`]. If `feasible` is false the vectors are empty.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BFlow {
    pub feasible: bool,
    /// Total cost `sum c_e f_e` (can exceed 64 bits).
    pub cost: i128,
    /// Flow on every edge, in input order.
    pub flow: Vec<i64>,
    /// Dual potential of every vertex.
    pub potential: Vec<i64>,
}

/// Solve the b-flow problem on `n` vertices. `b[v]` is the supply of `v`; every edge is
/// `(from, to, lower, upper, cost)` with `lower <= upper`. Self-loops are allowed.
///
/// # Panics
/// Panics on out-of-range endpoints, `b.len() != n` or `lower > upper`.
pub fn min_cost_b_flow(n: usize, b: &[i64], edges: &[(usize, usize, i64, i64, i64)]) -> BFlow {
    assert_eq!(b.len(), n, "need one supply per vertex");
    let infeasible = BFlow { feasible: false, cost: 0, flow: Vec::new(), potential: Vec::new() };
    if b.iter().map(|&x| x as i128).sum::<i128>() != 0 {
        return infeasible;
    }
    let m = edges.len();
    // Arc 2i: s -> t with capacity u - l, cost c; arc 2i+1: its reverse.
    let mut to = vec![0usize; 2 * m];
    let mut res = vec![0i64; 2 * m];
    let mut cost = vec![0i64; 2 * m];
    let mut excess: Vec<i64> = b.to_vec();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut loop_flow = vec![0i64; m];
    for (i, &(s, t, l, u, c)) in edges.iter().enumerate() {
        assert!(s < n && t < n, "edge {i} has an endpoint outside 0..{n}");
        assert!(l <= u, "edge {i}: lower bound exceeds upper bound");
        // shift out the lower bound: the l units leave s and arrive at t
        excess[s] -= l;
        excess[t] += l;
        if s == t {
            loop_flow[i] = if c < 0 { u - l } else { 0 };
            continue;
        }
        to[2 * i] = t;
        to[2 * i + 1] = s;
        res[2 * i] = u - l;
        cost[2 * i] = c;
        cost[2 * i + 1] = -c;
        adj[s].push(2 * i);
        adj[t].push(2 * i + 1);
    }
    // from[a] = tail of arc a
    let from = |a: usize| to[a ^ 1];
    let mut potential = vec![0i64; n];
    // Saturate negative-cost arcs so that all residual arcs have non-negative reduced cost.
    let push = |a: usize, amount: i64, res: &mut Vec<i64>, excess: &mut Vec<i64>| {
        res[a] -= amount;
        res[a ^ 1] += amount;
        excess[to[a ^ 1]] -= amount;
        excess[to[a]] += amount;
    };
    for i in 0..m {
        if edges[i].0 != edges[i].1 && cost[2 * i] < 0 {
            let r = res[2 * i];
            push(2 * i, r, &mut res, &mut excess);
        }
    }
    let max_excess = excess.iter().map(|x| x.abs()).max().unwrap_or(0).max(1);
    let mut delta = 1i64 << (63 - (max_excess as u64).leading_zeros());
    let mut dist = vec![i64::MAX; n];
    let mut prev_arc = vec![usize::MAX; n];
    loop {
        // Re-establish non-negative reduced costs on the delta-residual graph.
        for a in 0..2 * m {
            if res[a] >= delta && cost[a] + potential[from(a)] - potential[to[a]] < 0 {
                let r = res[a];
                push(a, r, &mut res, &mut excess);
            }
        }
        let mut blocked = vec![false; n];
        loop {
            let Some(k) = (0..n).find(|&v| excess[v] >= delta && !blocked[v]) else { break };
            // Dijkstra from k over arcs with residual >= delta, until a deficit node is settled.
            dist.iter_mut().for_each(|d| *d = i64::MAX);
            let mut heap = BinaryHeap::new();
            dist[k] = 0;
            heap.push(Reverse((0i64, k)));
            let mut settled: Vec<usize> = Vec::new();
            let mut target = usize::MAX;
            while let Some(Reverse((d, v))) = heap.pop() {
                if d > dist[v] {
                    continue;
                }
                settled.push(v);
                if excess[v] <= -delta {
                    target = v;
                    break;
                }
                for &a in &adj[v] {
                    if res[a] < delta {
                        continue;
                    }
                    let w = to[a];
                    let nd = d + cost[a] + potential[v] - potential[w];
                    if nd < dist[w] {
                        dist[w] = nd;
                        prev_arc[w] = a;
                        heap.push(Reverse((nd, w)));
                    }
                }
            }
            if target == usize::MAX {
                if delta == 1 {
                    return infeasible;
                }
                blocked[k] = true;
                continue;
            }
            let cap = dist[target];
            for v in 0..n {
                potential[v] += dist[v].min(cap);
            }
            let mut v = target;
            while v != k {
                let a = prev_arc[v];
                res[a] -= delta;
                res[a ^ 1] += delta;
                v = from(a);
            }
            excess[k] -= delta;
            excess[target] += delta;
        }
        if delta == 1 {
            break;
        }
        delta >>= 1;
    }
    let mut flow = Vec::with_capacity(m);
    let mut total: i128 = 0;
    for (i, &(s, t, l, _, c)) in edges.iter().enumerate() {
        let f = if s == t { l + loop_flow[i] } else { l + res[2 * i + 1] };
        total += c as i128 * f as i128;
        flow.push(f);
    }
    BFlow { feasible: true, cost: total, flow, potential }
}
