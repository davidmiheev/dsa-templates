//! Stress test: `Hld` against naive parent/depth/size computations, the structural
//! properties of heavy-light layouts, and path decomposition into position segments.

use heavy_light_decomposition::Hld;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn naive_tree(adj: &[Vec<usize>], root: usize) -> (Vec<Option<usize>>, Vec<usize>, Vec<usize>) {
    let n = adj.len();
    let (mut parent, mut depth, mut order) = (vec![None; n], vec![0; n], vec![root]);
    let mut i = 0;
    while i < order.len() {
        let v = order[i];
        i += 1;
        for &u in &adj[v] {
            if Some(u) != parent[v] {
                parent[u] = Some(v);
                depth[u] = depth[v] + 1;
                order.push(u);
            }
        }
    }
    let mut size = vec![1; n];
    for &v in order.iter().rev() {
        if let Some(p) = parent[v] {
            size[p] += size[v];
        }
    }
    (parent, depth, size)
}

fn path_vertices(parent: &[Option<usize>], depth: &[usize], mut u: usize, mut v: usize) -> Vec<usize> {
    let mut out = vec![];
    while u != v {
        if depth[u] >= depth[v] {
            out.push(u);
            u = parent[u].unwrap();
        } else {
            out.push(v);
            v = parent[v].unwrap();
        }
    }
    out.push(u);
    out
}

fn check_tree(adj: &[Vec<usize>], root: usize, rng: &mut Rng, ctx: &str, queries: usize) -> usize {
    let n = adj.len();
    let h = Hld::new(adj, root);
    let (parent, depth, size) = naive_tree(adj, root);
    assert_eq!((&h.parent, &h.depth, &h.size), (&parent, &depth, &size), "parent/depth/size {ctx}");

    // pos is a permutation of 0..n and every subtree occupies a contiguous block starting at pos[v].
    let mut seen = vec![false; n];
    for &p in &h.pos {
        assert!(p < n && !seen[p], "pos must be a permutation {ctx}");
        seen[p] = true;
    }
    let mut at = vec![0; n];
    for v in 0..n {
        at[h.pos[v]] = v;
    }
    for v in 0..n {
        let mut members: Vec<usize> = (h.pos[v]..h.pos[v] + size[v]).map(|p| at[p]).collect();
        members.sort_unstable();
        let mut want = vec![v];
        let mut stack = vec![v];
        while let Some(x) = stack.pop() {
            for &u in &adj[x] {
                if Some(u) != parent[x] {
                    want.push(u);
                    stack.push(u);
                }
            }
        }
        want.sort_unstable();
        assert_eq!(members, want, "subtree of {v} must be contiguous in pos order {ctx}");

        // chain structure: a vertex is either a chain head or continues its parent's chain
        match parent[v] {
            None => assert_eq!(h.head[v], v, "root heads its chain {ctx}"),
            Some(p) => {
                if h.head[v] == v {
                    assert_ne!(h.head[p], h.head[v], "a head starts a new chain {ctx}");
                } else {
                    assert_eq!(h.head[v], h.head[p], "chain continues through the parent {ctx}");
                    assert_eq!(h.pos[v], h.pos[p] + 1, "heavy child directly follows its parent {ctx}");
                    let biggest = adj[p].iter().filter(|&&u| Some(u) != parent[p]).map(|&u| size[u]).max().unwrap();
                    assert_eq!(size[v], biggest, "chain continues into a largest child {ctx}");
                }
            }
        }
    }

    // path_segments: disjoint segments covering exactly the path's positions, O(log n) of them.
    let limit = 2 * (usize::BITS - n.leading_zeros()) as usize + 2;
    for _ in 0..queries {
        let (u, v) = (rng.below(n), rng.below(n));
        let segs = h.path_segments(u, v);
        let mut got: Vec<usize> = segs.iter().flat_map(|&(lo, hi)| {
            assert!(lo < hi && hi <= n, "segment [{lo},{hi}) must be non-empty and in range {ctx}");
            lo..hi
        }).collect();
        got.sort_unstable();
        let mut want: Vec<usize> = path_vertices(&parent, &depth, u, v).into_iter().map(|x| h.pos[x]).collect();
        want.sort_unstable();
        assert_eq!(got, want, "segments of path {u}-{v} {ctx}");
        assert!(segs.len() <= limit, "{} segments for path {u}-{v} exceeds O(log n) bound {limit} {ctx}", segs.len());
    }
    n + queries
}

#[test]
fn matches_naive_on_small_trees() {
    let mut st = Stress::new("hld/naive");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 60) as usize;
        let adj = rng.tree(n, case % 4);
        let root = rng.below(n);
        let ops = check_tree(&adj, root, &mut rng, &format!("case {case} n={n} shape={}", case % 4), 30);
        st.case();
        st.ops(ops);
    }
    st.done();
}

#[test]
fn large_trees_on_a_big_stack() {
    // The crate's decomposition recurses along heavy paths, so deep trees need more than the
    // default 2 MiB test-thread stack; run on a 512 MiB stack.
    std::thread::Builder::new()
        .stack_size(512 << 20)
        .spawn(|| {
            let mut st = Stress::new("hld/large");
            for shape in 0..4 {
                let mut rng = Rng::new(100 + shape as u64);
                let n = if cfg!(debug_assertions) { 50_000 } else { 200_000 }; // smaller in debug builds
                let adj = rng.tree(n, shape);
                let ops = check_tree_light(&adj, &mut rng, &format!("shape {shape}"), 2_000);
                st.case();
                st.ops(ops);
            }
            st.done();
        })
        .unwrap()
        .join()
        .unwrap();
}

/// Cheaper variant for n = 200k: property checks that are O(n) or O(log n) per query.
fn check_tree_light(adj: &[Vec<usize>], rng: &mut Rng, ctx: &str, queries: usize) -> usize {
    let n = adj.len();
    let h = Hld::new(adj, 0);
    let (parent, depth, size) = naive_tree(adj, 0);
    assert_eq!((&h.parent, &h.depth, &h.size), (&parent, &depth, &size), "parent/depth/size {ctx}");
    let mut seen = vec![false; n];
    for &p in &h.pos {
        assert!(!seen[p], "pos must be a permutation {ctx}");
        seen[p] = true;
    }
    for v in 1..n {
        let p = parent[v].unwrap();
        if h.head[v] != v {
            assert_eq!((h.head[v], h.pos[v]), (h.head[p], h.pos[p] + 1), "chain step at {v} {ctx}");
        }
        // subtree block: the last position in [pos[v], pos[v] + size[v]) belongs to v's subtree
        assert!(h.pos[v] >= h.pos[p] + 1 && h.pos[v] + size[v] <= h.pos[p] + size[p], "nesting at {v} {ctx}");
    }
    let limit = 2 * (usize::BITS - n.leading_zeros()) as usize + 2;
    for _ in 0..queries {
        let (u, v) = (rng.below(n), rng.below(n));
        let segs = h.path_segments(u, v);
        let covered: usize = segs.iter().map(|&(lo, hi)| hi - lo).sum();
        let path = path_vertices(&parent, &depth, u, v);
        assert_eq!(covered, path.len(), "path {u}-{v} length {ctx}");
        assert!(path.iter().all(|&x| segs.iter().any(|&(lo, hi)| lo <= h.pos[x] && h.pos[x] < hi)), "cover {ctx}");
        assert!(segs.len() <= limit, "{} segments exceeds bound {limit} {ctx}", segs.len());
    }
    n + queries
}
