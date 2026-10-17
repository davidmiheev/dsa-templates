//! Euler tour tree: a dynamic forest whose vertices carry values of a [`MapMonoid`], supporting
//!
//! * [`link`](EulerTourTree::link) / [`cut`](EulerTourTree::cut) of edges,
//! * [`component_prod`](EulerTourTree::component_prod): the fold of all values in a tree,
//! * [`apply_component`](EulerTourTree::apply_component): apply an action to all values in a tree,
//! * [`get`](EulerTourTree::get) / [`set`](EulerTourTree::set) / [`apply`](EulerTourTree::apply)
//!   for a single vertex,
//! * the "subtree of `v` when `p` is its parent" versions
//!   [`subtree_prod`](EulerTourTree::subtree_prod) and
//!   [`apply_subtree`](EulerTourTree::apply_subtree), which cut the edge `(v, p)`, work on the
//!   component of `v` and link it back.
//!
//! Every tree is stored as the cyclic Euler tour of its vertices in a treap (expected
//! `O(log n)` per operation). The fold order inside a tree is unspecified, so the monoid must
//! be **commutative**, and the action must map the identity to the identity (use a monoid that
//! carries sizes, e.g. `(sum, count)` with "add" as the action).
//!
//! # Example
//!
//! ```
//! use algebra::{MapMonoid, Monoid};
//! use euler_tour_tree::EulerTourTree;
//!
//! enum SumLen {}
//! impl Monoid for SumLen {
//!     type S = (i64, i64); // (sum, number of vertices)
//!     fn identity() -> (i64, i64) { (0, 0) }
//!     fn op(a: &(i64, i64), b: &(i64, i64)) -> (i64, i64) { (a.0 + b.0, a.1 + b.1) }
//! }
//! enum Add {}
//! impl MapMonoid for Add {
//!     type M = SumLen;
//!     type F = i64;
//!     fn id_map() -> i64 { 0 }
//!     fn apply(f: &i64, x: &(i64, i64)) -> (i64, i64) { (x.0 + f * x.1, x.1) }
//!     fn compose(f: &i64, g: &i64) -> i64 { f + g }
//! }
//!
//! // Path 0 - 1 - 2 with values 1, 10, 100.
//! let mut t = EulerTourTree::<Add>::new(&[(1, 1), (10, 1), (100, 1)]);
//! t.link(0, 1);
//! t.link(1, 2);
//! assert_eq!(t.subtree_prod(1, 0).0, 110); // vertex 1 with parent 0: {1, 2}
//! t.apply_subtree(1, 0, 5);               // +5 on vertices 1 and 2
//! assert_eq!(t.component_prod(0).0, 1 + 15 + 105);
//! ```

use algebra::{MapMonoid, Monoid};
use std::collections::HashMap;

type S<A> = <<A as MapMonoid>::M as Monoid>::S;

const NIL: u32 = 0;

struct Node<A: MapMonoid> {
    left: u32,
    right: u32,
    parent: u32,
    size: u32,
    priority: u32,
    value: S<A>,
    agg: S<A>,
    lazy: Option<A::F>,
}

/// See the [module documentation](self).
pub struct EulerTourTree<A: MapMonoid> {
    n: usize,
    /// Node `v + 1` is the token of vertex `v`; node 0 is the empty sequence; later nodes are
    /// edge tokens.
    nodes: Vec<Node<A>>,
    free: Vec<u32>,
    /// Edge token of the directed edge `(u, v)`.
    edge: HashMap<(u32, u32), u32>,
    rng: u64,
}

impl<A: MapMonoid> EulerTourTree<A> {
    /// A forest of `values.len()` isolated vertices carrying `values`.
    pub fn new(values: &[S<A>]) -> Self {
        let id = <A::M as Monoid>::identity();
        let mut t = EulerTourTree { n: values.len(), nodes: Vec::with_capacity(3 * values.len() + 1), free: Vec::new(), edge: HashMap::new(), rng: 0x9E37_79B9_7F4A_7C15 };
        t.nodes.push(Node { left: NIL, right: NIL, parent: NIL, size: 0, priority: 0, value: id.clone(), agg: id, lazy: None });
        for v in values {
            let priority = t.next_priority();
            t.nodes.push(Node { left: NIL, right: NIL, parent: NIL, size: 1, priority, value: v.clone(), agg: v.clone(), lazy: None });
        }
        t
    }

    /// Number of vertices.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    fn next_priority(&mut self) -> u32 {
        self.rng ^= self.rng << 13;
        self.rng ^= self.rng >> 7;
        self.rng ^= self.rng << 17;
        (self.rng >> 32) as u32
    }

    fn check(&self, v: usize) -> u32 {
        assert!(v < self.n, "vertex {v} out of range for {} vertices", self.n);
        v as u32 + 1
    }

    fn is_vertex(&self, x: u32) -> bool {
        (x as usize) <= self.n
    }

    fn alloc_edge_token(&mut self) -> u32 {
        let id = <A::M as Monoid>::identity();
        let priority = self.next_priority();
        let node = Node { left: NIL, right: NIL, parent: NIL, size: 1, priority, value: id.clone(), agg: id, lazy: None };
        if let Some(x) = self.free.pop() {
            self.nodes[x as usize] = node;
            x
        } else {
            self.nodes.push(node);
            (self.nodes.len() - 1) as u32
        }
    }

    fn apply_to(&mut self, x: u32, f: &A::F) {
        if x == NIL {
            return;
        }
        let is_vertex = self.is_vertex(x);
        let node = &mut self.nodes[x as usize];
        node.agg = A::apply(f, &node.agg);
        if is_vertex {
            node.value = A::apply(f, &node.value);
        }
        node.lazy = Some(match node.lazy.take() {
            Some(g) => A::compose(f, &g),
            None => f.clone(),
        });
    }

    fn push(&mut self, x: u32) {
        if let Some(f) = self.nodes[x as usize].lazy.take() {
            let (l, r) = (self.nodes[x as usize].left, self.nodes[x as usize].right);
            self.apply_to(l, &f);
            self.apply_to(r, &f);
        }
    }

    fn pull(&mut self, x: u32) {
        let (l, r) = (self.nodes[x as usize].left, self.nodes[x as usize].right);
        let size = self.nodes[l as usize].size + 1 + self.nodes[r as usize].size;
        let agg = <A::M as Monoid>::op(&<A::M as Monoid>::op(&self.nodes[l as usize].agg, &self.nodes[x as usize].value), &self.nodes[r as usize].agg);
        let node = &mut self.nodes[x as usize];
        node.size = size;
        node.agg = agg;
        if l != NIL {
            self.nodes[l as usize].parent = x;
        }
        if r != NIL {
            self.nodes[r as usize].parent = x;
        }
    }

    fn merge(&mut self, a: u32, b: u32) -> u32 {
        if a == NIL {
            return b;
        }
        if b == NIL {
            return a;
        }
        if self.nodes[a as usize].priority > self.nodes[b as usize].priority {
            self.push(a);
            let r = self.nodes[a as usize].right;
            let m = self.merge(r, b);
            self.nodes[a as usize].right = m;
            self.pull(a);
            a
        } else {
            self.push(b);
            let l = self.nodes[b as usize].left;
            let m = self.merge(a, l);
            self.nodes[b as usize].left = m;
            self.pull(b);
            b
        }
    }

    /// Split into the first `k` elements and the rest; both results get parent `NIL`.
    fn split(&mut self, t: u32, k: u32) -> (u32, u32) {
        if t == NIL {
            return (NIL, NIL);
        }
        self.push(t);
        let l = self.nodes[t as usize].left;
        let lsize = self.nodes[l as usize].size;
        let (a, b) = if k <= lsize {
            let (a, b) = self.split(l, k);
            self.nodes[t as usize].left = b;
            self.pull(t);
            (a, t)
        } else {
            let r = self.nodes[t as usize].right;
            let (a, b) = self.split(r, k - lsize - 1);
            self.nodes[t as usize].right = a;
            self.pull(t);
            (t, b)
        };
        self.nodes[a as usize].parent = NIL;
        self.nodes[b as usize].parent = NIL;
        (a, b)
    }

    fn root_of(&self, mut x: u32) -> u32 {
        while self.nodes[x as usize].parent != NIL {
            x = self.nodes[x as usize].parent;
        }
        x
    }

    /// Number of tokens before `x` in its sequence.
    fn rank_of(&self, mut x: u32) -> u32 {
        let mut rank = self.nodes[self.nodes[x as usize].left as usize].size;
        while self.nodes[x as usize].parent != NIL {
            let p = self.nodes[x as usize].parent;
            if self.nodes[p as usize].right == x {
                rank += self.nodes[self.nodes[p as usize].left as usize].size + 1;
            }
            x = p;
        }
        rank
    }

    /// Push the lazy tags on the path from the root of `x`'s treap down to `x`.
    fn push_path(&mut self, x: u32) {
        let mut path = vec![x];
        let mut y = x;
        while self.nodes[y as usize].parent != NIL {
            y = self.nodes[y as usize].parent;
            path.push(y);
        }
        for &y in path.iter().rev() {
            self.push(y);
        }
    }

    /// Rotate the tour of `v`'s tree so that it starts at `v`; returns the new treap root.
    fn reroot(&mut self, x: u32) -> u32 {
        let root = self.root_of(x);
        let k = self.rank_of(x);
        let (a, b) = self.split(root, k);
        self.merge(b, a)
    }

    /// Whether `u` and `v` are in the same tree.
    pub fn connected(&self, u: usize, v: usize) -> bool {
        let (x, y) = (self.check(u), self.check(v));
        self.root_of(x) == self.root_of(y)
    }

    /// Add the edge `(u, v)`.
    ///
    /// # Panics
    /// Panics if `u` and `v` are already connected.
    pub fn link(&mut self, u: usize, v: usize) {
        assert!(!self.connected(u, v), "link({u}, {v}): the vertices are already connected");
        let (x, y) = (self.check(u), self.check(v));
        let tu = self.reroot(x);
        let tv = self.reroot(y);
        let (e_uv, e_vu) = (self.alloc_edge_token(), self.alloc_edge_token());
        self.edge.insert((u as u32, v as u32), e_uv);
        self.edge.insert((v as u32, u as u32), e_vu);
        let t = self.merge(tu, e_uv);
        let t = self.merge(t, tv);
        self.merge(t, e_vu);
    }

    /// Remove the edge `(u, v)`.
    ///
    /// # Panics
    /// Panics if there is no such edge.
    pub fn cut(&mut self, u: usize, v: usize) {
        self.check(u);
        self.check(v);
        let (&e1, &e2) = (
            self.edge.get(&(u as u32, v as u32)).unwrap_or_else(|| panic!("cut({u}, {v}): no such edge")),
            self.edge.get(&(v as u32, u as u32)).unwrap_or_else(|| panic!("cut({u}, {v}): no such edge")),
        );
        self.edge.remove(&(u as u32, v as u32));
        self.edge.remove(&(v as u32, u as u32));
        let (mut a, mut b) = (e1, e2);
        let (mut pa, mut pb) = (self.rank_of(a), self.rank_of(b));
        if pa > pb {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut pa, &mut pb);
        }
        let root = self.root_of(a);
        // root = X a Y b Z
        let (x, rest) = self.split(root, pa);
        let (_, rest) = self.split(rest, 1);
        let (y, rest) = self.split(rest, pb - pa - 1);
        let (_, z) = self.split(rest, 1);
        self.merge(x, z);
        let _ = y; // `y` is now a treap of its own (parent NIL)
        self.free.push(e1);
        self.free.push(e2);
    }

    /// The fold of the values of all vertices in the tree of `v`.
    pub fn component_prod(&self, v: usize) -> S<A> {
        let x = self.check(v);
        self.nodes[self.root_of(x) as usize].agg.clone()
    }

    /// Apply `f` to the value of every vertex in the tree of `v`.
    pub fn apply_component(&mut self, v: usize, f: A::F) {
        let x = self.check(v);
        let root = self.root_of(x);
        self.apply_to(root, &f);
    }

    /// The value of `v`.
    pub fn get(&mut self, v: usize) -> S<A> {
        let x = self.check(v);
        self.push_path(x);
        self.nodes[x as usize].value.clone()
    }

    /// Replace the value of `v`.
    pub fn set(&mut self, v: usize, value: S<A>) {
        let x = self.check(v);
        self.push_path(x);
        self.nodes[x as usize].value = value;
        let mut y = x;
        while y != NIL {
            self.pull(y);
            y = self.nodes[y as usize].parent;
        }
    }

    /// Apply `f` to the value of `v` only.
    pub fn apply(&mut self, v: usize, f: A::F) {
        let cur = self.get(v);
        self.set(v, A::apply(&f, &cur));
    }

    /// The fold over the subtree of `v` when `parent` is regarded as its parent, i.e. the
    /// component of `v` after removing the edge `(v, parent)`.
    ///
    /// # Panics
    /// Panics if `(v, parent)` is not an edge.
    pub fn subtree_prod(&mut self, v: usize, parent: usize) -> S<A> {
        self.cut(v, parent);
        let r = self.component_prod(v);
        self.link(v, parent);
        r
    }

    /// Apply `f` to every vertex of the subtree of `v` when `parent` is regarded as its parent.
    ///
    /// # Panics
    /// Panics if `(v, parent)` is not an edge.
    pub fn apply_subtree(&mut self, v: usize, parent: usize, f: A::F) {
        self.cut(v, parent);
        self.apply_component(v, f);
        self.link(v, parent);
    }
}
