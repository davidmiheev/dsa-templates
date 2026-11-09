//! Edge colouring of a bipartite multigraph with exactly `Δ` (maximum degree) colours — the
//! optimum by König's theorem.
//!
//! 1. **Regularise.** Vertices of each side are packed greedily into groups whose total degree is
//!    at most `Δ` (merging vertices never invalidates a colouring), the two sides are padded to the
//!    same number of groups and dummy edges are added until the multigraph is `Δ`-regular. There
//!    are `O(m / Δ)` groups, so `O(m)` edges.
//! 2. **Split.** An even-degree regular graph is cut by Euler tours into two regular graphs of half
//!    the degree (alternate the edges along every closed trail). An odd-degree one has a perfect
//!    matching (Hall); that matching is one colour class
//!    ([`matching::hopcroft_karp`]), after which the degree is even.
//! 3. Recurse until the degree is 1: one more colour class.
//!
//! The work is `O(m log Δ)` for the splits plus `O(m sqrt n)` per odd level.
//!
//! # Example
//!
//! ```
//! use edge_coloring::bipartite_edge_coloring;
//!
//! // K_{2,2} plus a parallel edge: maximum degree 3.
//! let (k, color) = bipartite_edge_coloring(2, 2, &[(0, 0), (0, 1), (1, 0), (1, 1), (0, 0)]);
//! assert_eq!(k, 3);
//! assert_eq!(color.len(), 5);
//! assert_ne!(color[0], color[4]); // parallel edges differ
//! ```

use matching::hopcroft_karp;

/// Colours the edges `(left, right)` of the bipartite multigraph with `l` left and `r` right
/// vertices. Returns `(k, color)` with `k` the maximum degree and `color[i] < k` for edge `i`;
/// two edges sharing an endpoint always differ. An edgeless graph gets `k = 0`.
///
/// # Panics
/// Panics if an endpoint is out of range.
pub fn bipartite_edge_coloring(l: usize, r: usize, edges: &[(usize, usize)]) -> (usize, Vec<usize>) {
    let m = edges.len();
    let (mut dl, mut dr) = (vec![0usize; l], vec![0usize; r]);
    for &(a, b) in edges {
        assert!(a < l && b < r, "edge ({a}, {b}) has an endpoint outside the {l} x {r} sides");
        dl[a] += 1;
        dr[b] += 1;
    }
    let delta = dl.iter().chain(dr.iter()).copied().max().unwrap_or(0);
    if delta == 0 {
        return (0, Vec::new());
    }
    // Pack vertices into groups of total degree <= delta.
    let pack = |deg: &[usize]| -> (Vec<usize>, Vec<usize>) {
        let (mut group, mut load) = (vec![0usize; deg.len()], vec![0usize]);
        for (v, &d) in deg.iter().enumerate() {
            if *load.last().unwrap() + d > delta {
                load.push(0);
            }
            group[v] = load.len() - 1;
            *load.last_mut().unwrap() += d;
        }
        (group, load)
    };
    let (gl, mut load_l) = pack(&dl);
    let (gr, mut load_r) = pack(&dr);
    let groups = load_l.len().max(load_r.len());
    load_l.resize(groups, 0);
    load_r.resize(groups, 0);
    // All edges of the regular multigraph: real ones first, then the dummy padding.
    let mut es: Vec<(usize, usize)> = edges.iter().map(|&(a, b)| (gl[a], gr[b])).collect();
    let (mut i, mut j) = (0, 0);
    loop {
        while i < groups && load_l[i] == delta {
            i += 1;
        }
        while j < groups && load_r[j] == delta {
            j += 1;
        }
        if i == groups || j == groups {
            break;
        }
        let add = (delta - load_l[i]).min(delta - load_r[j]);
        for _ in 0..add {
            es.push((i, j));
        }
        load_l[i] += add;
        load_r[j] += add;
    }
    let mut color = vec![usize::MAX; es.len()];
    colour_regular(groups, delta, (0..es.len()).collect(), 0, &es, &mut color);
    color.truncate(m);
    (delta, color)
}

/// Colours the `delta`-regular bipartite multigraph formed by `ids` (indices into `es`) on
/// `groups + groups` vertices with the colours `offset..offset + delta`.
fn colour_regular(groups: usize, delta: usize, ids: Vec<usize>, offset: usize, es: &[(usize, usize)], color: &mut [usize]) {
    let mut work = vec![(ids, delta, offset)];
    while let Some((mut ids, mut delta, mut offset)) = work.pop() {
        if delta == 0 {
            continue;
        }
        if delta % 2 == 1 {
            // Peel a perfect matching as one colour class.
            let local: Vec<(usize, usize)> = ids.iter().map(|&e| es[e]).collect();
            let mt = hopcroft_karp(groups, groups, &local);
            assert_eq!(mt.size, groups, "a regular bipartite multigraph has a perfect matching");
            let mut taken = vec![false; groups];
            let mut rest = Vec::with_capacity(ids.len() - groups);
            for (idx, &e) in ids.iter().enumerate() {
                let (a, b) = local[idx];
                if !taken[a] && mt.left_to_right[a] == b {
                    taken[a] = true;
                    color[e] = offset;
                } else {
                    rest.push(e);
                }
            }
            ids = rest;
            offset += 1;
            delta -= 1;
            if delta == 0 {
                continue;
            }
        }
        // Even degree: split into two halves with Euler tours.
        let (half_a, half_b) = euler_split(groups, &ids, es);
        work.push((half_a, delta / 2, offset));
        work.push((half_b, delta / 2, offset + delta / 2));
    }
    // Degree-1 graphs are perfect matchings: handled by the odd branch above (delta = 1).
}

/// Splits an even-regular bipartite multigraph into two halves, alternating edges along trails.
fn euler_split(groups: usize, ids: &[usize], es: &[(usize, usize)]) -> (Vec<usize>, Vec<usize>) {
    // vertices 0..groups are left, groups..2*groups are right
    let nv = 2 * groups;
    let mut start = vec![0usize; nv + 1];
    for &e in ids {
        let (a, b) = es[e];
        start[a + 1] += 1;
        start[groups + b + 1] += 1;
    }
    for v in 0..nv {
        start[v + 1] += start[v];
    }
    let mut slot = start.clone();
    let mut arc = vec![(0usize, 0usize); 2 * ids.len()]; // (edge position, other end)
    for (pos, &e) in ids.iter().enumerate() {
        let (a, b) = es[e];
        arc[slot[a]] = (pos, groups + b);
        slot[a] += 1;
        arc[slot[groups + b]] = (pos, a);
        slot[groups + b] += 1;
    }
    let mut next = start[..nv].to_vec();
    let mut used = vec![false; ids.len()];
    let (mut a_half, mut b_half) = (Vec::with_capacity(ids.len() / 2), Vec::with_capacity(ids.len() / 2));
    for s in 0..nv {
        loop {
            // advance to an unused edge at s
            while next[s] < start[s + 1] && used[arc[next[s]].0] {
                next[s] += 1;
            }
            if next[s] == start[s + 1] {
                break;
            }
            // walk a closed trail from s
            let (mut v, mut parity) = (s, 0usize);
            loop {
                while next[v] < start[v + 1] && used[arc[next[v]].0] {
                    next[v] += 1;
                }
                if next[v] == start[v + 1] {
                    break;
                }
                let (pos, w) = arc[next[v]];
                used[pos] = true;
                if parity == 0 {
                    a_half.push(ids[pos]);
                } else {
                    b_half.push(ids[pos]);
                }
                parity ^= 1;
                v = w;
            }
        }
    }
    (a_half, b_half)
}
