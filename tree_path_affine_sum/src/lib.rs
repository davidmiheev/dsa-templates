//! Sum over all vertices of a function composed along the path from a chosen root.
//!
//! Every vertex `v` carries a number `a_v` and every edge `e` an affine function
//! `f_e(x) = b_e * x + c_e` (modulo the const-generic prime `MOD`). For a root `r` and a vertex
//! `y`, let `P(r, y) = f_{e_0}(f_{e_1}(... f_{e_k}(a_y)))` where `e_0, ..., e_k` are the edges
//! on the path from `r` to `y`, starting at `r`. [`TreePathAffineSum`] maintains
//! `sum_y P(r, y)`:
//!
//! * [`set_vertex`](TreePathAffineSum::set_vertex) / [`set_edge`](TreePathAffineSum::set_edge),
//! * [`sum_from_root`](TreePathAffineSum::sum_from_root) for the root fixed at construction,
//! * [`sum_from`](TreePathAffineSum::sum_from) for any root.
//!
//! All operations are `O(log^2 n)`: a "dynamic DP" over heavy-light chains. Along a heavy chain
//! the value of a vertex is an affine function of the value of its heavy child, so chain
//! suffixes are products of affine maps in a segment tree; light children are summed into the
//! parent. For roots other than the construction root, a second segment tree holds the maps
//! that carry "everything outside the subtree" down the chains.
//!
//! # Example
//!
//! ```
//! use tree_path_affine_sum::TreePathAffineSum;
//!
//! // Path 0 - 1 - 2, a = [1, 2, 3]; edge 0 = (0, 1): f(x) = 2x + 1, edge 1 = (1, 2): f(x) = x + 5.
//! let mut t = TreePathAffineSum::<998_244_353>::new(
//!     &[1, 2, 3],
//!     &[(0, 1, 2, 1), (1, 2, 1, 5)],
//!     0,
//! );
//! // From root 0: P(0,0)=1, P(0,1)=2*2+1=5, P(0,2)=2*(3+5)+1=17.
//! assert_eq!(t.sum_from_root(), 1 + 5 + 17);
//! // From root 2: P(2,2)=3, P(2,1)=f_1(2)=7, P(2,0)=f_1(f_0(1))=1*(2*1+1)+5=8.
//! assert_eq!(t.sum_from(2), 3 + 7 + 8);
//! t.set_vertex(1, 10);
//! assert_eq!(t.sum_from_root(), 1 + 21 + 17);
//! ```

use heavy_light_decomposition::Hld;

type Aff = (u64, u64);

const NONE: usize = usize::MAX;

/// See the [module documentation](self).
pub struct TreePathAffineSum<const MOD: u64> {
    n: usize,
    root: usize,
    hld: Hld,
    a: Vec<u64>,
    /// `(b, c)` of the edge from each non-root vertex to its parent.
    edge_fn: Vec<Aff>,
    /// Edge index -> the child endpoint.
    edge_child: Vec<usize>,
    heavy: Vec<usize>,
    /// By head: position of the last vertex of its chain.
    chain_end: Vec<usize>,
    /// For every vertex: sum of the contributions of its light children.
    light: Vec<u64>,
    /// For every chain head: its current contribution `b * S + c * size` to its parent.
    contribution: Vec<u64>,
    down: AffSeg<MOD>,
    up: AffSeg<MOD>,
}

/// Iterative segment tree over affine maps. `reversed = false`: the fold of positions
/// `l..r` is `f_l ∘ f_{l+1} ∘ ... ∘ f_{r-1}`; `reversed = true`: `f_{r-1} ∘ ... ∘ f_l`.
struct AffSeg<const MOD: u64> {
    size: usize,
    reversed: bool,
    tree: Vec<Aff>,
}

fn compose<const MOD: u64>(f: Aff, g: Aff) -> Aff {
    // f ∘ g
    (f.0 * g.0 % MOD, (f.0 * g.1 + f.1) % MOD)
}

impl<const MOD: u64> AffSeg<MOD> {
    fn new(items: &[Aff], reversed: bool) -> Self {
        let size = items.len().max(1).next_power_of_two();
        let mut tree = vec![(1, 0); 2 * size];
        tree[size..size + items.len()].copy_from_slice(items);
        let mut s = AffSeg { size, reversed, tree };
        for i in (1..size).rev() {
            s.tree[i] = s.combine(s.tree[2 * i], s.tree[2 * i + 1]);
        }
        s
    }

    fn combine(&self, left: Aff, right: Aff) -> Aff {
        if self.reversed {
            compose::<MOD>(right, left)
        } else {
            compose::<MOD>(left, right)
        }
    }

    fn set(&mut self, i: usize, v: Aff) {
        let mut i = i + self.size;
        self.tree[i] = v;
        while i > 1 {
            i >>= 1;
            self.tree[i] = self.combine(self.tree[2 * i], self.tree[2 * i + 1]);
        }
    }

    fn prod(&self, l: usize, r: usize) -> Aff {
        let (mut l, mut r) = (l + self.size, r + self.size);
        let (mut left, mut right) = ((1, 0), (1, 0));
        while l < r {
            if l & 1 == 1 {
                left = self.combine(left, self.tree[l]);
                l += 1;
            }
            if r & 1 == 1 {
                r -= 1;
                right = self.combine(self.tree[r], right);
            }
            l >>= 1;
            r >>= 1;
        }
        self.combine(left, right)
    }
}

impl<const MOD: u64> TreePathAffineSum<MOD> {
    /// Build over a tree on `a.len()` vertices with `edges[e] = (u, v, b_e, c_e)`; the root used
    /// by [`sum_from_root`](Self::sum_from_root) is `root`. All numbers must be `< MOD`.
    ///
    /// # Panics
    /// Panics if the edges do not form a tree or `root` is out of range.
    pub fn new(a: &[u64], edges: &[(usize, usize, u64, u64)], root: usize) -> Self {
        let n = a.len();
        assert!(MOD < 1 << 31, "the modulus must be below 2^31");
        let plain: Vec<(usize, usize)> = edges.iter().map(|&(u, v, _, _)| (u, v)).collect();
        let hld = Hld::from_edges(n, &plain, root);
        let mut edge_fn = vec![(0, 0); n];
        let mut edge_child = Vec::with_capacity(edges.len());
        for &(u, v, b, c) in edges {
            let child = if hld.parent[v] == u { v } else { u };
            edge_fn[child] = (b % MOD, c % MOD);
            edge_child.push(child);
        }
        // Heavy child: the child that follows its parent in the preorder, if it is on its chain.
        let mut heavy = vec![NONE; n];
        let mut chain_end = vec![0usize; n];
        for (v, h) in heavy.iter_mut().enumerate() {
            let p = hld.pos[v];
            if p + 1 < n {
                let w = hld.order[p + 1];
                if hld.parent[w] == v && hld.head[w] == hld.head[v] {
                    *h = w;
                }
            }
        }
        for p in 0..n {
            let v = hld.order[p];
            chain_end[hld.head[v]] = p; // positions ascend, the last assignment wins
        }
        let mut t = TreePathAffineSum {
            n,
            root,
            a: a.iter().map(|&x| x % MOD).collect(),
            edge_fn,
            edge_child,
            heavy,
            chain_end,
            light: vec![0; n],
            contribution: vec![0; n],
            down: AffSeg::new(&[], false),
            up: AffSeg::new(&[], true),
            hld,
        };
        t.build();
        t
    }

    fn size_of(&self, v: usize) -> u64 {
        self.hld.size[v] as u64
    }

    /// The contribution `b * s + c * count` of a subtree with sum `s` and `count` vertices
    /// across the edge above `x`.
    fn across(&self, x: usize, s: u64, count: u64) -> u64 {
        let (b, c) = self.edge_fn[x];
        (b * s + c * (count % MOD)) % MOD
    }

    fn build(&mut self) {
        let n = self.n;
        // Subtree sums bottom-up (children have larger positions than their parent).
        let mut s = vec![0u64; n];
        for p in (0..n).rev() {
            let v = self.hld.order[p];
            let mut total = self.a[v] + s[v];
            total %= MOD;
            s[v] = total;
            if v != self.root {
                let par = self.hld.parent[v];
                let g = self.across(v, total, self.size_of(v));
                s[par] = (s[par] + g) % MOD;
                if self.heavy[par] != v {
                    self.light[par] = (self.light[par] + g) % MOD;
                    self.contribution[v] = g;
                }
            }
        }
        let mut down_items = vec![(1u64, 0u64); n];
        let mut up_items = vec![(1u64, 0u64); n];
        for v in 0..n {
            down_items[self.hld.pos[v]] = self.down_map(v);
            if self.hld.head[v] != v {
                up_items[self.hld.pos[v]] = self.up_map_heavy(self.hld.parent[v], v);
            }
        }
        self.down = AffSeg::new(&down_items, false);
        self.up = AffSeg::new(&up_items, true);
    }

    /// `S_v = alpha * S_heavy + beta`.
    fn down_map(&self, v: usize) -> Aff {
        let base = (self.a[v] + self.light[v]) % MOD;
        match self.heavy[v] {
            NONE => (0, base),
            h => {
                let (b, c) = self.edge_fn[h];
                (b, (base + c * (self.size_of(h) % MOD)) % MOD)
            }
        }
    }

    /// The map `U(h) = A * U(p) + B` for the heavy child `h` of `p`.
    fn up_map_heavy(&self, p: usize, h: usize) -> Aff {
        let (b, c) = self.edge_fn[h];
        let outside = (self.n as u64 - self.size_of(h)) % MOD;
        (b, (b * ((self.a[p] + self.light[p]) % MOD) + c * outside) % MOD)
    }

    /// Sum over the subtree of `v` of the values seen from `v` (rooted at `self.root`).
    fn subtree_sum(&self, v: usize) -> u64 {
        let h = self.hld.head[v];
        self.down.prod(self.hld.pos[v], self.chain_end[h] + 1).1
    }

    /// Refresh the segment-tree entries of `v` and propagate the changed contribution of its
    /// chain to the parent chains.
    fn refresh_from(&mut self, mut v: usize) {
        loop {
            self.down.set(self.hld.pos[v], self.down_map(v));
            let h = self.heavy[v];
            if h != NONE {
                self.up.set(self.hld.pos[h], self.up_map_heavy(v, h));
            }
            let head = self.hld.head[v];
            if head == self.root {
                return;
            }
            let g = self.across(head, self.subtree_sum(head), self.size_of(head));
            let parent = self.hld.parent[head];
            self.light[parent] = (self.light[parent] + g + MOD - self.contribution[head]) % MOD;
            self.contribution[head] = g;
            v = parent;
        }
    }

    /// Set `a_v = x`.
    pub fn set_vertex(&mut self, v: usize, x: u64) {
        assert!(v < self.n, "vertex {v} out of range for {} vertices", self.n);
        self.a[v] = x % MOD;
        self.refresh_from(v);
    }

    /// Set `f_e(x) = b x + c` for the edge with index `e` (in the order given to `new`).
    pub fn set_edge(&mut self, e: usize, b: u64, c: u64) {
        assert!(e < self.edge_child.len(), "edge {e} out of range");
        let child = self.edge_child[e];
        self.edge_fn[child] = (b % MOD, c % MOD);
        let parent = self.hld.parent[child];
        if self.heavy[parent] == child {
            // The map of the parent and the outside-map of the child both depend on the edge.
            self.refresh_from(parent);
            self.up.set(self.hld.pos[child], self.up_map_heavy(parent, child));
        } else {
            // A light child: its contribution to the parent changes.
            let g = self.across(child, self.subtree_sum(child), self.size_of(child));
            self.light[parent] = (self.light[parent] + g + MOD - self.contribution[child]) % MOD;
            self.contribution[child] = g;
            self.refresh_from(parent);
        }
    }

    /// `sum_y P(root, y)` for the root given to [`new`](Self::new).
    pub fn sum_from_root(&self) -> u64 {
        self.subtree_sum(self.root)
    }

    /// `sum_y P(r, y)` for an arbitrary root `r`.
    ///
    /// # Panics
    /// Panics if `r` is out of range.
    pub fn sum_from(&self, r: usize) -> u64 {
        assert!(r < self.n, "vertex {r} out of range for {} vertices", self.n);
        let inside = self.subtree_sum(r);
        // U(r) = h_r(h_{parent}(... h_{c}(0))): collect the maps from r upwards.
        let mut acc: Aff = (1, 0);
        let mut x = r;
        loop {
            let head = self.hld.head[x];
            if x != head {
                // Heavy-child maps of the chain pieces head+1 ..= x, deepest outermost.
                acc = compose::<MOD>(acc, self.up.prod(self.hld.pos[head] + 1, self.hld.pos[x] + 1));
            }
            if head == self.root {
                break;
            }
            let p = self.hld.parent[head];
            // The map of the light child `head` of `p`: everything around `p` except `head`'s
            // subtree, which includes the heavy child's contribution.
            let heavy = self.heavy[p];
            let heavy_part = if heavy == NONE { 0 } else { self.across(heavy, self.subtree_sum(heavy), self.size_of(heavy)) };
            let others = (self.a[p] + self.light[p] + MOD - self.contribution[head] + heavy_part) % MOD;
            let (b, c) = self.edge_fn[head];
            let outside = (self.n as u64 - self.size_of(head)) % MOD;
            let h = (b, (b * others + c * outside) % MOD);
            acc = compose::<MOD>(acc, h);
            x = p;
        }
        (inside + acc.1) % MOD
    }
}
