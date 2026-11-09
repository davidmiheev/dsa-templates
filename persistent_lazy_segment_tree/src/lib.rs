//! Fully persistent lazy segment tree.
//!
//! Every update returns the id of a new *version*; all earlier versions stay valid and can be
//! queried or updated again. Updates copy the `O(log n)` nodes they touch (path copying) and
//! share everything else, so `k` updates cost `O(k log n)` time and memory.
//!
//! * [`PersistentLazySegmentTree::apply_range`]: apply a [`MapMonoid`] action to `[l, r)`,
//! * [`PersistentLazySegmentTree::set`]: replace one element,
//! * [`PersistentLazySegmentTree::copy_range`]: make `[l, r)` of a version equal to the same
//!   range of *another* version, by sharing that version's subtrees,
//! * [`PersistentLazySegmentTree::prod`] / [`get`](PersistentLazySegmentTree::get): read-only.
//!
//! Pending actions are pushed down while descending on an update (copying the children they
//! land on), so the actions need not commute; queries do not modify anything and carry the
//! accumulated pending action down instead.
//!
//! # Example
//!
//! ```
//! use algebra::{MapMonoid, Monoid};
//! use persistent_lazy_segment_tree::PersistentLazySegmentTree;
//!
//! enum Sum {}
//! impl Monoid for Sum {
//!     type S = (i64, i64); // (sum, length)
//!     fn identity() -> (i64, i64) { (0, 0) }
//!     fn op(a: &(i64, i64), b: &(i64, i64)) -> (i64, i64) { (a.0 + b.0, a.1 + b.1) }
//! }
//! enum Add {}
//! impl MapMonoid for Add {
//!     type M = Sum;
//!     type F = i64;
//!     fn id_map() -> i64 { 0 }
//!     fn apply(f: &i64, x: &(i64, i64)) -> (i64, i64) { (x.0 + f * x.1, x.1) }
//!     fn compose(f: &i64, g: &i64) -> i64 { f + g }
//! }
//!
//! let mut t = PersistentLazySegmentTree::<Add>::from_slice(&[(1, 1), (2, 1), (3, 1), (4, 1)]);
//! let v1 = t.apply_range(0, 0, 3, 10); // version 1: [11, 12, 13, 4]
//! let v2 = t.copy_range(v1, 0, 1, 3);  // version 2: [11, 2, 3, 4] (range [1, 3) from version 0)
//! assert_eq!(t.prod(v1, 0, 4).0, 40);
//! assert_eq!(t.prod(v2, 0, 4).0, 20);
//! assert_eq!(t.prod(0, 0, 4).0, 10); // the original is untouched
//! ```

use algebra::{MapMonoid, Monoid};

type S<A> = <<A as MapMonoid>::M as Monoid>::S;

struct Node<A: MapMonoid> {
    left: u32,
    right: u32,
    /// Aggregate of the subtree, including this node's own pending action.
    agg: S<A>,
    /// Pending action for the children (already reflected in `agg`).
    lazy: A::F,
    has_lazy: bool,
}

impl<A: MapMonoid> Clone for Node<A> {
    fn clone(&self) -> Self {
        Node { left: self.left, right: self.right, agg: self.agg.clone(), lazy: self.lazy.clone(), has_lazy: self.has_lazy }
    }
}

/// Persistent lazy segment tree over the monoid `A::M` with range actions `A::F`.
pub struct PersistentLazySegmentTree<A: MapMonoid> {
    n: usize,
    nodes: Vec<Node<A>>,
    /// `roots[v]` is the root of version `v`.
    roots: Vec<u32>,
}

impl<A: MapMonoid> PersistentLazySegmentTree<A> {
    /// Version 0 holds `v`.
    ///
    /// # Panics
    /// Panics if `v` is empty.
    pub fn from_slice(v: &[S<A>]) -> Self {
        assert!(!v.is_empty(), "a persistent segment tree needs at least one element");
        let mut t = PersistentLazySegmentTree { n: v.len(), nodes: Vec::with_capacity(2 * v.len()), roots: Vec::new() };
        let root = t.build(v, 0, v.len());
        t.roots.push(root);
        t
    }

    fn build(&mut self, v: &[S<A>], lo: usize, hi: usize) -> u32 {
        if hi - lo == 1 {
            return self.alloc(Node { left: 0, right: 0, agg: v[lo].clone(), lazy: A::id_map(), has_lazy: false });
        }
        let mid = (lo + hi) / 2;
        let (l, r) = (self.build(v, lo, mid), self.build(v, mid, hi));
        let agg = <A::M as Monoid>::op(&self.nodes[l as usize].agg, &self.nodes[r as usize].agg);
        self.alloc(Node { left: l, right: r, agg, lazy: A::id_map(), has_lazy: false })
    }

    fn alloc(&mut self, node: Node<A>) -> u32 {
        self.nodes.push(node);
        (self.nodes.len() - 1) as u32
    }

    /// Number of versions created so far (version 0 is the initial one).
    pub fn num_versions(&self) -> usize {
        self.roots.len()
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    /// Total number of tree nodes allocated across all versions (a memory indicator).
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    fn check_version(&self, v: usize) -> u32 {
        assert!(v < self.roots.len(), "unknown version {v}");
        self.roots[v]
    }

    fn check_range(&self, l: usize, r: usize) {
        assert!(l <= r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
    }

    fn apply_to(&mut self, u: u32, f: &A::F) {
        let n = &mut self.nodes[u as usize];
        n.agg = A::apply(f, &n.agg);
        n.lazy = if n.has_lazy { A::compose(f, &n.lazy) } else { f.clone() };
        n.has_lazy = true;
    }

    /// `u` itself if it was created by the current operation (id >= `start`), else a copy.
    fn own(&mut self, u: u32, start: u32) -> u32 {
        if u >= start {
            u
        } else {
            let copy = self.nodes[u as usize].clone();
            self.alloc(copy)
        }
    }

    /// Push the pending action of the owned node `u` into (owned copies of) its children.
    fn push_down(&mut self, u: u32, start: u32) {
        if !self.nodes[u as usize].has_lazy {
            return;
        }
        let f = self.nodes[u as usize].lazy.clone();
        let (l, r) = (self.nodes[u as usize].left, self.nodes[u as usize].right);
        let (l, r) = (self.own(l, start), self.own(r, start));
        self.apply_to(l, &f);
        self.apply_to(r, &f);
        let n = &mut self.nodes[u as usize];
        n.left = l;
        n.right = r;
        n.lazy = A::id_map();
        n.has_lazy = false;
    }

    fn pull(&mut self, u: u32) {
        let (l, r) = (self.nodes[u as usize].left as usize, self.nodes[u as usize].right as usize);
        self.nodes[u as usize].agg = <A::M as Monoid>::op(&self.nodes[l].agg, &self.nodes[r].agg);
    }

    /// Apply `f` to `[l, r)` of version `version`; returns the new version's id.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()` and the version exists.
    pub fn apply_range(&mut self, version: usize, l: usize, r: usize, f: A::F) -> usize {
        let root = self.check_version(version);
        self.check_range(l, r);
        let new_root = if l == r {
            root
        } else {
            let start = self.nodes.len() as u32;
            self.update(root, start, 0, self.n, l, r, &f)
        };
        self.roots.push(new_root);
        self.roots.len() - 1
    }

    #[allow(clippy::too_many_arguments)]
    fn update(&mut self, u: u32, start: u32, lo: usize, hi: usize, l: usize, r: usize, f: &A::F) -> u32 {
        let u = self.own(u, start);
        if l <= lo && hi <= r {
            self.apply_to(u, f);
            return u;
        }
        self.push_down(u, start);
        let mid = (lo + hi) / 2;
        if l < mid {
            let c = self.update(self.nodes[u as usize].left, start, lo, mid, l, r, f);
            self.nodes[u as usize].left = c;
        }
        if r > mid {
            let c = self.update(self.nodes[u as usize].right, start, mid, hi, l, r, f);
            self.nodes[u as usize].right = c;
        }
        self.pull(u);
        u
    }

    /// Replace element `i` of `version` by `x`; returns the new version's id.
    ///
    /// # Panics
    /// Panics if `i >= len()` or the version does not exist.
    pub fn set(&mut self, version: usize, i: usize, x: S<A>) -> usize {
        let root = self.check_version(version);
        assert!(i < self.n, "index {i} out of range for length {}", self.n);
        let start = self.nodes.len() as u32;
        let new_root = self.set_rec(root, start, 0, self.n, i, x);
        self.roots.push(new_root);
        self.roots.len() - 1
    }

    fn set_rec(&mut self, u: u32, start: u32, lo: usize, hi: usize, i: usize, x: S<A>) -> u32 {
        let u = self.own(u, start);
        if hi - lo == 1 {
            let n = &mut self.nodes[u as usize];
            n.agg = x;
            n.has_lazy = false;
            n.lazy = A::id_map();
            return u;
        }
        self.push_down(u, start);
        let mid = (lo + hi) / 2;
        if i < mid {
            let c = self.set_rec(self.nodes[u as usize].left, start, lo, mid, i, x);
            self.nodes[u as usize].left = c;
        } else {
            let c = self.set_rec(self.nodes[u as usize].right, start, mid, hi, i, x);
            self.nodes[u as usize].right = c;
        }
        self.pull(u);
        u
    }

    /// The version `dst` with `[l, r)` replaced by the same range of version `src`; returns the
    /// new version's id. The copied subtrees are shared, not duplicated.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()` and both versions exist.
    pub fn copy_range(&mut self, dst: usize, src: usize, l: usize, r: usize) -> usize {
        let (dst_root, src_root) = (self.check_version(dst), self.check_version(src));
        self.check_range(l, r);
        let new_root = if l == r {
            dst_root
        } else {
            let start = self.nodes.len() as u32;
            self.graft(dst_root, src_root, start, 0, self.n, l, r)
        };
        self.roots.push(new_root);
        self.roots.len() - 1
    }

    /// Walk `u` (destination) and `v` (source) together; `v` must be a consistent view, i.e. all
    /// pending actions of its ancestors are already pushed into it.
    #[allow(clippy::too_many_arguments)]
    fn graft(&mut self, u: u32, v: u32, start: u32, lo: usize, hi: usize, l: usize, r: usize) -> u32 {
        if l <= lo && hi <= r {
            return v;
        }
        let u = self.own(u, start);
        self.push_down(u, start);
        // Views of the source's children with the source's pending action applied.
        let (mut vl, mut vr) = (self.nodes[v as usize].left, self.nodes[v as usize].right);
        if self.nodes[v as usize].has_lazy {
            let f = self.nodes[v as usize].lazy.clone();
            vl = self.own(vl, start);
            vr = self.own(vr, start);
            self.apply_to(vl, &f);
            self.apply_to(vr, &f);
        }
        let mid = (lo + hi) / 2;
        if l < mid {
            let c = self.graft(self.nodes[u as usize].left, vl, start, lo, mid, l, r);
            self.nodes[u as usize].left = c;
        }
        if r > mid {
            let c = self.graft(self.nodes[u as usize].right, vr, start, mid, hi, l, r);
            self.nodes[u as usize].right = c;
        }
        self.pull(u);
        u
    }

    /// The fold of `[l, r)` in `version` (the identity for an empty range). Read-only.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()` and the version exists.
    pub fn prod(&self, version: usize, l: usize, r: usize) -> S<A> {
        let root = self.check_version(version);
        self.check_range(l, r);
        self.prod_rec(root, 0, self.n, l, r, None)
    }

    /// `pending` is the composed action of the ancestors that has not been pushed into `u`.
    fn prod_rec(&self, u: u32, lo: usize, hi: usize, l: usize, r: usize, pending: Option<&A::F>) -> S<A> {
        if r <= lo || hi <= l {
            return <A::M as Monoid>::identity();
        }
        let node = &self.nodes[u as usize];
        if l <= lo && hi <= r {
            return match pending {
                Some(f) => A::apply(f, &node.agg),
                None => node.agg.clone(),
            };
        }
        // The children have not seen `node.lazy` (nearer, so applied first) nor `pending`.
        let own;
        let below = if node.has_lazy {
            own = match pending {
                Some(f) => A::compose(f, &node.lazy),
                None => node.lazy.clone(),
            };
            Some(&own)
        } else {
            pending
        };
        let mid = (lo + hi) / 2;
        let left = self.prod_rec(node.left, lo, mid, l, r, below);
        let right = self.prod_rec(node.right, mid, hi, l, r, below);
        <A::M as Monoid>::op(&left, &right)
    }

    /// Element `i` of `version`.
    pub fn get(&self, version: usize, i: usize) -> S<A> {
        assert!(i < self.n, "index {i} out of range for length {}", self.n);
        self.prod(version, i, i + 1)
    }
}
