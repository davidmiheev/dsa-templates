//! Offline incremental strongly connected components in `O(m log m)` SCC work.
//!
//! For a sequence of edge insertions this computes, for every edge, the **time** — the index of
//! the insertion after which — its two endpoints lie in one SCC ([`NEVER`] if they never do).
//! Merging the endpoints of the edges in order of that time with a union-find then reproduces
//! the SCC partition after every insertion.
//!
//! Divide and conquer over time: for the time range `[lo, hi]` and the edges whose answer lies in
//! it (endpoints already contracted by the SCCs of time `lo - 1`), compute the SCCs of the edges
//! inserted by `mid`. Edges inside an SCC have an answer `<= mid` and go left; the others go right
//! with their endpoints contracted into the new SCCs.
//!
//! # Example
//!
//! ```
//! use incremental_scc::{incremental_scc_times, NEVER};
//!
//! // 0 -> 1, 1 -> 2, 2 -> 0 closes the cycle at time 2; 2 -> 3 is never on a cycle.
//! let t = incremental_scc_times(4, &[(0, 1), (1, 2), (2, 0), (2, 3)]);
//! assert_eq!(t, vec![2, 2, 2, NEVER]);
//! ```

use scc::Scc;

/// "This edge is never inside an SCC."
pub const NEVER: usize = usize::MAX;

struct Ctx {
    ans: Vec<usize>,
    stamp: Vec<usize>,
    local: Vec<usize>,
    m: usize,
    round: usize,
}

/// For edge `i` of the directed multigraph `edges` on `0..n`, the smallest `t >= i` such that
/// the endpoints are in one SCC of the graph formed by edges `0..=t`, or [`NEVER`].
pub fn incremental_scc_times(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let m = edges.len();
    let mut ctx = Ctx { ans: vec![NEVER; m], stamp: vec![usize::MAX; n], local: vec![0; n], m, round: 0 };
    let mut list = Vec::with_capacity(m);
    for (i, &(u, v)) in edges.iter().enumerate() {
        assert!(u < n && v < n, "edge {i} has an endpoint outside 0..{n}");
        if u == v {
            ctx.ans[i] = i;
        } else {
            list.push((u, v, i));
        }
    }
    solve(0, m, list, &mut ctx);
    ctx.ans
}

fn solve(lo: usize, hi: usize, list: Vec<(usize, usize, usize)>, ctx: &mut Ctx) {
    if list.is_empty() {
        return;
    }
    if lo == hi {
        for &(_, _, id) in &list {
            ctx.ans[id] = if lo == ctx.m { NEVER } else { lo };
        }
        return;
    }
    let mid = (lo + hi) / 2;
    // Relabel the vertices touched by the edges inserted by `mid`.
    ctx.round += 1;
    let round = ctx.round;
    let mut verts: Vec<usize> = Vec::new();
    let mut local_edges: Vec<(usize, usize)> = Vec::new();
    for &(u, v, id) in &list {
        if id > mid {
            continue;
        }
        for x in [u, v] {
            if ctx.stamp[x] != round {
                ctx.stamp[x] = round;
                ctx.local[x] = verts.len();
                verts.push(x);
            }
        }
        local_edges.push((ctx.local[u], ctx.local[v]));
    }
    let comps = Scc::from_edges(verts.len(), &local_edges);
    // representative vertex of every component
    let mut rep = vec![usize::MAX; comps.count];
    for (i, &x) in verts.iter().enumerate() {
        if rep[comps.comp[i]] == usize::MAX {
            rep[comps.comp[i]] = x;
        }
    }
    let (mut left, mut right) = (Vec::new(), Vec::new());
    for &(u, v, id) in &list {
        if id <= mid {
            let (cu, cv) = (comps.comp[ctx.local[u]], comps.comp[ctx.local[v]]);
            if cu == cv {
                left.push((u, v, id));
                continue;
            }
        }
        let map = |x: usize, ctx: &Ctx| if ctx.stamp[x] == round { rep[comps.comp[ctx.local[x]]] } else { x };
        right.push((map(u, ctx), map(v, ctx), id));
    }
    drop(list);
    // Edges whose endpoints merged into one vertex without being inside a local SCC cannot
    // exist: both endpoints of an edge inserted by `mid` are in `verts`, and an edge between two
    // vertices of one SCC is a `left` edge. An edge inserted after `mid` may however join two
    // vertices of one SCC; it is a self-loop now and becomes part of the SCC at its own time.
    let mut loops = Vec::new();
    right.retain(|&(u, v, id)| {
        if u == v {
            loops.push(id);
            false
        } else {
            true
        }
    });
    for id in loops {
        ctx.ans[id] = id;
    }
    solve(lo, mid, left, ctx);
    solve(mid + 1, hi, right, ctx);
}
