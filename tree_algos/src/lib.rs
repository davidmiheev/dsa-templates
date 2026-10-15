//! Small, self-contained tree algorithms.
//!
//! * [`cartesian_tree`]: the min-rooted Cartesian tree of a sequence in `O(n)`.
//! * [`tree_diameter`]: weighted diameter and one longest path in `O(n)`.
//! * [`rooted_isomorphism_classes`]: classes of isomorphic rooted subtrees in `O(n log n)`.
//! * [`min_inversion_order`]: tree-ordered schedule minimizing `sum c_i d_j` over pairs `j` before `i`.
//! * [`rerooting`]: the answer of a tree DP for *every* choice of root in `O(n)` merges.
//!
//! Parent arrays use [`NONE`] as the parent of a root.
//!
//! # Example
//!
//! ```
//! use tree_algos::{cartesian_tree, tree_diameter, NONE};
//!
//! assert_eq!(cartesian_tree(&[3, 1, 2]), vec![1, NONE, 1]);
//! let (len, path) = tree_diameter(4, &[(0, 1, 5u64), (1, 2, 1), (1, 3, 7)]);
//! assert_eq!((len, path.len()), (12, 3)); // 0 - 1 - 3
//! ```

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::ops::Add;

/// The parent of a root.
pub const NONE: usize = usize::MAX;

/// The Cartesian tree of `a` as a parent array: the root holds the smallest element (the
/// leftmost one among equals), and the in-order traversal is `0..n`.
pub fn cartesian_tree<T: Ord>(a: &[T]) -> Vec<usize> {
    let n = a.len();
    let mut parent = vec![NONE; n];
    let mut stack: Vec<usize> = Vec::new();
    for i in 0..n {
        let mut last = NONE;
        while let Some(&top) = stack.last() {
            if a[top] > a[i] {
                last = top;
                stack.pop();
            } else {
                break;
            }
        }
        if last != NONE {
            parent[last] = i;
        }
        if let Some(&top) = stack.last() {
            parent[i] = top;
        }
        stack.push(i);
    }
    parent
}

/// The diameter of a weighted tree on `n` vertices (non-negative weights): its length and the
/// vertices of one longest path, in order. Returns `(0, [0])` for a single vertex.
///
/// # Panics
/// Panics if `n == 0`, an endpoint is out of range, or the edges do not form a tree.
pub fn tree_diameter<W: Copy + Ord + Default + Add<Output = W>>(n: usize, edges: &[(usize, usize, W)]) -> (W, Vec<usize>) {
    assert!(n > 0, "the tree must have at least one vertex");
    assert_eq!(edges.len() + 1, n, "a tree on {n} vertices has {} edges", n - 1);
    let mut offsets = vec![0usize; n + 2];
    for &(u, v, _) in edges {
        assert!(u < n && v < n, "edge endpoint out of range");
        offsets[u + 2] += 1;
        offsets[v + 2] += 1;
    }
    for i in 2..n + 2 {
        offsets[i] += offsets[i - 1];
    }
    let mut adj = vec![(0usize, W::default()); 2 * edges.len()];
    for &(u, v, w) in edges {
        adj[offsets[u + 1]] = (v, w);
        offsets[u + 1] += 1;
        adj[offsets[v + 1]] = (u, w);
        offsets[v + 1] += 1;
    }
    // Farthest vertex from `s`, with parents for path recovery.
    let sweep = |s: usize| -> (usize, W, Vec<usize>) {
        let mut parent = vec![NONE; n];
        let mut dist = vec![W::default(); n];
        let mut order = vec![s];
        let mut i = 0;
        while i < order.len() {
            let v = order[i];
            i += 1;
            for &(u, w) in &adj[offsets[v]..offsets[v + 1]] {
                if u != parent[v] {
                    parent[u] = v;
                    dist[u] = dist[v] + w;
                    order.push(u);
                }
            }
        }
        assert_eq!(order.len(), n, "the graph is not connected");
        let far = *order.iter().max_by_key(|&&v| dist[v]).expect("non-empty");
        (far, dist[far], parent)
    };
    let (a, _, _) = sweep(0);
    let (b, len, parent) = sweep(a);
    let mut path = vec![b];
    while *path.last().unwrap() != a {
        path.push(parent[*path.last().unwrap()]);
    }
    path.reverse();
    (len, path)
}

/// Classify the subtrees of a rooted tree by isomorphism. Returns the number of classes `k` and
/// a class id in `0..k` for the subtree rooted at every vertex; two vertices get the same id
/// exactly when their subtrees are isomorphic as rooted trees.
///
/// # Panics
/// Panics if the parent array does not describe a tree with the single root `root`.
pub fn rooted_isomorphism_classes(parent: &[usize], root: usize) -> (usize, Vec<usize>) {
    let n = parent.len();
    assert!(root < n, "root out of range");
    let mut offsets = vec![0usize; n + 2];
    for (v, &p) in parent.iter().enumerate() {
        if v != root {
            assert!(p < n, "parent out of range");
            offsets[p + 2] += 1;
        }
    }
    for i in 2..n + 2 {
        offsets[i] += offsets[i - 1];
    }
    let mut children = vec![0usize; n - 1];
    for (v, &p) in parent.iter().enumerate() {
        if v != root {
            children[offsets[p + 1]] = v;
            offsets[p + 1] += 1;
        }
    }
    offsets.pop();
    let mut order = vec![root];
    let mut i = 0;
    while i < order.len() {
        let v = order[i];
        i += 1;
        order.extend_from_slice(&children[offsets[v]..offsets[v + 1]]);
    }
    assert_eq!(order.len(), n, "the parent array does not describe a tree rooted at {root}");
    let mut class = vec![0usize; n];
    let mut ids: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut key: Vec<u32> = Vec::new();
    for &v in order.iter().rev() {
        key.clear();
        key.extend(children[offsets[v]..offsets[v + 1]].iter().map(|&c| class[c] as u32));
        key.sort_unstable();
        let next = ids.len();
        class[v] = match ids.get(&key) {
            Some(&id) => id,
            None => {
                ids.insert(key.clone(), next);
                next
            }
        };
    }
    (ids.len(), class)
}

/// Tree DP for every root at once ("rerooting").
///
/// For a root `r` the DP value of a vertex `v` is
/// `vertex(v, fold over neighbours u (away from r) of edge(dp(u), e, u, v))`, where `fold` is
/// the commutative-or-not monoid `(identity, merge)` (the neighbours of `v` are merged in the
/// order of the edge list, which only matters for non-commutative `merge`), `edge(value, e, u, v)`
/// carries the DP value `value` of the subtree at `u` across the edge `e = edges[e] = {u, v}` to
/// `v`, and `vertex(acc, v)` finishes the value of `v` from the merged contributions.
///
/// Returns `ans[r]` = the DP value of `r` when the tree is rooted at `r`, for every `r`.
/// No inverse of `merge` is needed (prefix/suffix products are used), `O(n)` closure calls.
///
/// ```
/// use tree_algos::rerooting;
///
/// // Farthest vertex distance (eccentricity) of every vertex of the path 0 - 1 - 2 - 3.
/// let edges = [(0, 1), (1, 2), (2, 3)];
/// let ecc = rerooting(4, &edges, 0u32, |a, b| *a.max(b), |&d, _, _, _| d + 1, |&acc, _| acc);
/// assert_eq!(ecc, vec![3, 2, 2, 3]);
/// ```
///
/// # Panics
/// Panics if the edges do not form a tree on `n >= 1` vertices.
pub fn rerooting<V, Merge, Edge, Vertex>(n: usize, edges: &[(usize, usize)], identity: V, merge: Merge, edge: Edge, vertex: Vertex) -> Vec<V>
where
    V: Clone,
    Merge: Fn(&V, &V) -> V,
    Edge: Fn(&V, usize, usize, usize) -> V,
    Vertex: Fn(&V, usize) -> V,
{
    assert!(n > 0, "the tree must have at least one vertex");
    assert_eq!(edges.len() + 1, n, "a tree on {n} vertices has {} edges", n - 1);
    // Adjacency in CSR form: (neighbour, edge index).
    let mut offsets = vec![0usize; n + 2];
    for &(u, v) in edges {
        assert!(u < n && v < n, "edge endpoint out of range");
        offsets[u + 2] += 1;
        offsets[v + 2] += 1;
    }
    for i in 2..n + 2 {
        offsets[i] += offsets[i - 1];
    }
    let mut adj = vec![(0usize, 0usize); 2 * edges.len()];
    for (e, &(u, v)) in edges.iter().enumerate() {
        adj[offsets[u + 1]] = (v, e);
        offsets[u + 1] += 1;
        adj[offsets[v + 1]] = (u, e);
        offsets[v + 1] += 1;
    }
    offsets.pop();
    // Root at 0.
    let mut parent = vec![NONE; n];
    let mut parent_edge = vec![NONE; n];
    let mut order = vec![0usize];
    let mut i = 0;
    while i < order.len() {
        let v = order[i];
        i += 1;
        for &(u, e) in &adj[offsets[v]..offsets[v + 1]] {
            if u != parent[v] {
                parent[u] = v;
                parent_edge[u] = e;
                order.push(u);
            }
        }
    }
    assert_eq!(order.len(), n, "the graph is not connected");
    // Bottom-up: `up_value[c]` = contribution of the subtree of `c` to its parent.
    let mut up_value: Vec<Option<V>> = vec![None; n];
    for &v in order.iter().rev() {
        let mut acc = identity.clone();
        for &(u, _) in &adj[offsets[v]..offsets[v + 1]] {
            if u != parent[v] {
                acc = merge(&acc, up_value[u].as_ref().expect("children come first"));
            }
        }
        let value = vertex(&acc, v);
        if parent[v] != NONE {
            up_value[v] = Some(edge(&value, parent_edge[v], v, parent[v]));
        }
    }
    // Top-down: `from_parent[c]` = contribution of everything outside the subtree of `c`.
    let mut from_parent: Vec<Option<V>> = vec![None; n];
    let mut ans: Vec<Option<V>> = vec![None; n];
    for &v in &order {
        let nbrs = &adj[offsets[v]..offsets[v + 1]];
        // Contributions in adjacency order (the parent's slot holds `from_parent[v]`).
        let contrib: Vec<&V> = nbrs.iter().map(|&(u, _)| if u == parent[v] { from_parent[v].as_ref().expect("parent first") } else { up_value[u].as_ref().expect("computed") }).collect();
        let mut suffix = vec![identity.clone(); contrib.len() + 1];
        for k in (0..contrib.len()).rev() {
            suffix[k] = merge(contrib[k], &suffix[k + 1]);
        }
        let mut prefix = identity.clone();
        let mut sent: Vec<(usize, V)> = Vec::new();
        for (k, &(u, e)) in nbrs.iter().enumerate() {
            if u != parent[v] {
                let rest = vertex(&merge(&prefix, &suffix[k + 1]), v);
                sent.push((u, edge(&rest, e, v, u)));
            }
            prefix = merge(&prefix, contrib[k]);
        }
        ans[v] = Some(vertex(&prefix, v));
        drop(contrib);
        for (u, value) in sent {
            from_parent[u] = Some(value);
        }
    }
    ans.into_iter().map(|a| a.expect("every vertex is reached")).collect()
}

/// An order of the vertices of a rooted tree in which every vertex comes after its parent, that
/// minimizes `sum over pairs (j before i) of c[order[i]] * d[order[j]]` (a precedence-constrained
/// single-machine scheduling problem). Returns the order and the minimum cost.
///
/// `parent[v]` is the parent of `v` (`NONE` for the root `root`); `c`, `d` are non-negative and
/// their products must fit in `u128` sums (the cost itself is returned as `u128`).
///
/// Greedy (Horn / Adolphson-Hu): repeatedly take the non-root group with the smallest
/// `D / C` ratio and append it to its parent's group, where `C`, `D` are sums over the group.
/// `O(n log n)`.
///
/// ```
/// use tree_algos::{min_inversion_order, NONE};
/// // 0 is the root of 1 and 2. Vertex 2 has a small d / c ratio, so it goes first:
/// // cost = c[2] * d[0] + c[1] * (d[0] + d[2]) = 0 + 1 * 1.
/// let (order, cost) = min_inversion_order(&[NONE, 0, 0], 0, &[1, 1, 5], &[0, 10, 1]);
/// assert_eq!(order, vec![0, 2, 1]);
/// assert_eq!(cost, 1);
/// ```
///
/// # Panics
/// Panics if the parent array does not describe a tree rooted at `root`, or the slices differ in
/// length.
pub fn min_inversion_order(parent: &[usize], root: usize, c: &[u64], d: &[u64]) -> (Vec<usize>, u128) {
    let n = parent.len();
    assert!(c.len() == n && d.len() == n, "c and d need one entry per vertex");
    assert!(root < n, "root out of range");
    // A group is identified by its top vertex. `ds`/`cs` hold group sums; the sequences are
    // singly linked lists (`next`) with head = the top vertex and a tail pointer.
    let mut cs: Vec<u128> = c.iter().map(|&x| x as u128).collect();
    let mut ds: Vec<u128> = d.iter().map(|&x| x as u128).collect();
    let mut tail: Vec<usize> = (0..n).collect();
    let mut next = vec![NONE; n];
    let mut dsu: Vec<usize> = (0..n).collect(); // group top reachable by find()
    fn find(dsu: &mut [usize], mut x: usize) -> usize {
        let mut r = x;
        while dsu[r] != r {
            r = dsu[r];
        }
        while dsu[x] != r {
            let nx = dsu[x];
            dsu[x] = r;
            x = nx;
        }
        r
    }
    // Min-heap on D / C (cross-multiplied), stale entries skipped by comparing sums.
    #[derive(Clone, Copy)]
    struct Entry {
        d: u128,
        c: u128,
        v: usize,
    }
    impl PartialEq for Entry {
        fn eq(&self, o: &Self) -> bool {
            self.cmp(o) == Ordering::Equal
        }
    }
    impl Eq for Entry {}
    impl PartialOrd for Entry {
        fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
            Some(self.cmp(o))
        }
    }
    impl Ord for Entry {
        // BinaryHeap is a max-heap: the "greatest" entry has the smallest d / c. A group with
        // c == 0 and d == 0 is neutral and goes first; c == 0 with d > 0 has ratio infinity.
        fn cmp(&self, o: &Self) -> Ordering {
            let neutral = |e: &Entry| e.c == 0 && e.d == 0;
            match (neutral(self), neutral(o)) {
                (true, true) => return o.v.cmp(&self.v),
                (true, false) => return Ordering::Greater,
                (false, true) => return Ordering::Less,
                _ => {}
            }
            // self is greater when self.d / self.c < o.d / o.c, i.e. self.d * o.c < o.d * self.c.
            (o.d * self.c).cmp(&(self.d * o.c)).then(o.v.cmp(&self.v))
        }
    }
    let mut heap: BinaryHeap<Entry> = BinaryHeap::new();
    let mut seen_root = false;
    for v in 0..n {
        if v == root {
            assert!(parent[v] == NONE, "the root must have parent NONE");
            seen_root = true;
        } else {
            assert!(parent[v] < n, "parent out of range");
            heap.push(Entry { d: ds[v], c: cs[v], v });
        }
    }
    assert!(seen_root);
    let mut merged = 0;
    while let Some(e) = heap.pop() {
        if dsu[e.v] != e.v || ds[e.v] != e.d || cs[e.v] != e.c {
            continue; // stale
        }
        let p = find(&mut dsu, parent[e.v]);
        assert!(p != e.v, "the parent array contains a cycle");
        // Append the group of e.v behind the group of p.
        next[tail[p]] = e.v;
        tail[p] = tail[e.v];
        ds[p] += ds[e.v];
        cs[p] += cs[e.v];
        dsu[e.v] = p;
        merged += 1;
        if p != root {
            heap.push(Entry { d: ds[p], c: cs[p], v: p });
        }
    }
    assert_eq!(merged, n - 1, "the parent array does not describe a tree rooted at {root}");
    let mut order = Vec::with_capacity(n);
    let mut v = root;
    while v != NONE {
        order.push(v);
        v = next[v];
    }
    let (mut sum_d, mut cost) = (0u128, 0u128);
    for &v in &order {
        cost += c[v] as u128 * sum_d;
        sum_d += d[v] as u128;
    }
    (order, cost)
}
