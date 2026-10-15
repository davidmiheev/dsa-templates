//! Segment trees over a huge index range that allocate only the nodes they touch.
//!
//! The array has `n` elements (`n` up to `2^32`, e.g. `10^9`), all equal to a given `default`
//! value at the start. Nodes are created on demand, so `q` operations use `O(q log n)` memory.
//! Both trees are non-persistent and work with non-commutative monoids.
//!
//! * [`DynamicSegmentTree`]: point `set` / `get` and range `prod` over a [`Monoid`].
//! * [`DynamicLazySegmentTree`]: additionally `apply_range` with a [`MapMonoid`] action.
//!
//! The aggregate of an untouched block of `2^k` default elements is precomputed, so the
//! monoid's `op` is never asked for anything but real values.
//!
//! # Example
//!
//! ```
//! use algebra::Monoid;
//! use dynamic_segment_tree::DynamicSegmentTree;
//!
//! enum Sum {}
//! impl Monoid for Sum {
//!     type S = i64;
//!     fn identity() -> i64 { 0 }
//!     fn op(a: &i64, b: &i64) -> i64 { a + b }
//! }
//!
//! let mut t = DynamicSegmentTree::<Sum>::new(1_000_000_000, 0);
//! t.set(5, 7);
//! t.set(999_999_999, 3);
//! assert_eq!(t.prod(0, 1_000_000_000), 10);
//! assert_eq!(t.prod(6, 999_999_999), 0);
//! assert_eq!(t.get(5), 7);
//! ```

use algebra::{MapMonoid, Monoid};

type S<A> = <<A as MapMonoid>::M as Monoid>::S;

const NIL: u32 = 0;

/// `def[k]` = the aggregate of `2^k` default elements.
fn default_aggregates<M: Monoid>(bits: u32, default: &M::S) -> Vec<M::S> {
    let mut def = vec![default.clone()];
    for k in 0..bits as usize {
        let next = M::op(&def[k], &def[k]);
        def.push(next);
    }
    def
}

fn bits_for(n: usize) -> u32 {
    assert!(n >= 1 && (n as u64) <= (1 << 32), "length must be in 1..=2^32, got {n}");
    n.next_power_of_two().trailing_zeros()
}

struct Node<S> {
    left: u32,
    right: u32,
    agg: S,
}

/// Point-update, range-fold segment tree over `n <= 2^32` elements, all initially `default`.
pub struct DynamicSegmentTree<M: Monoid> {
    n: usize,
    bits: u32,
    def: Vec<M::S>,
    /// `nodes[0]` is unused (`NIL`).
    nodes: Vec<Node<M::S>>,
    root: u32,
}

impl<M: Monoid> DynamicSegmentTree<M> {
    /// A tree of `n` copies of `default`.
    ///
    /// # Panics
    /// Panics unless `1 <= n <= 2^32`.
    pub fn new(n: usize, default: M::S) -> Self {
        let bits = bits_for(n);
        let def = default_aggregates::<M>(bits, &default);
        let nodes = vec![Node { left: NIL, right: NIL, agg: default }];
        DynamicSegmentTree { n, bits, def, nodes, root: NIL }
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    /// Number of allocated nodes (a memory indicator).
    pub fn node_count(&self) -> usize {
        self.nodes.len() - 1
    }

    fn agg_of(&self, u: u32, level: u32) -> &M::S {
        if u == NIL {
            &self.def[level as usize]
        } else {
            &self.nodes[u as usize].agg
        }
    }

    /// Set element `i` to `x`.
    ///
    /// # Panics
    /// Panics if `i >= len()`.
    pub fn set(&mut self, i: usize, x: M::S) {
        assert!(i < self.n, "index {i} out of range for length {}", self.n);
        // Walk down, creating missing nodes, remembering the path.
        let mut path = [NIL; 33];
        let mut depth = 0;
        let mut u = self.root;
        let mut slot_is_root = true;
        let mut parent = NIL;
        let mut went_right = false;
        for level in (0..=self.bits).rev() {
            if u == NIL {
                let agg = self.def[level as usize].clone();
                self.nodes.push(Node { left: NIL, right: NIL, agg });
                u = (self.nodes.len() - 1) as u32;
                if slot_is_root {
                    self.root = u;
                } else if went_right {
                    self.nodes[parent as usize].right = u;
                } else {
                    self.nodes[parent as usize].left = u;
                }
            }
            path[depth] = u;
            depth += 1;
            if level == 0 {
                break;
            }
            slot_is_root = false;
            parent = u;
            went_right = (i >> (level - 1)) & 1 == 1;
            u = if went_right { self.nodes[u as usize].right } else { self.nodes[u as usize].left };
        }
        let leaf = path[depth - 1] as usize;
        self.nodes[leaf].agg = x;
        for level in 1..=self.bits as usize {
            let u = path[depth - 1 - level] as usize;
            let (l, r) = (self.nodes[u].left, self.nodes[u].right);
            let agg = M::op(self.agg_of(l, level as u32 - 1), self.agg_of(r, level as u32 - 1));
            self.nodes[u].agg = agg;
        }
    }

    /// Element `i`.
    pub fn get(&self, i: usize) -> M::S {
        assert!(i < self.n, "index {i} out of range for length {}", self.n);
        self.prod(i, i + 1)
    }

    /// The fold of `[l, r)` (the identity for an empty range).
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn prod(&self, l: usize, r: usize) -> M::S {
        assert!(l <= r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
        self.prod_rec(self.root, self.bits, 0, l, r)
    }

    fn prod_rec(&self, u: u32, level: u32, lo: usize, l: usize, r: usize) -> M::S {
        let hi = lo + (1usize << level);
        if r <= lo || hi <= l {
            return M::identity();
        }
        if l <= lo && hi <= r {
            return self.agg_of(u, level).clone();
        }
        let (left, right) = if u == NIL { (NIL, NIL) } else { (self.nodes[u as usize].left, self.nodes[u as usize].right) };
        let mid = lo + (1usize << (level - 1));
        M::op(&self.prod_rec(left, level - 1, lo, l, r), &self.prod_rec(right, level - 1, mid, l, r))
    }
}

struct LazyNode<A: MapMonoid> {
    left: u32,
    right: u32,
    agg: S<A>,
    /// Pending action for the children (already reflected in `agg`); `None` = identity.
    lazy: Option<A::F>,
}

/// Lazy segment tree (range action, range fold, point set) over `n <= 2^32` elements, all
/// initially `default`.
pub struct DynamicLazySegmentTree<A: MapMonoid> {
    n: usize,
    bits: u32,
    def: Vec<S<A>>,
    nodes: Vec<LazyNode<A>>,
    root: u32,
}

impl<A: MapMonoid> DynamicLazySegmentTree<A> {
    /// A tree of `n` copies of `default`.
    ///
    /// # Panics
    /// Panics unless `1 <= n <= 2^32`.
    pub fn new(n: usize, default: S<A>) -> Self {
        let bits = bits_for(n);
        let def = default_aggregates::<A::M>(bits, &default);
        let nodes = vec![LazyNode { left: NIL, right: NIL, agg: default, lazy: None }];
        DynamicLazySegmentTree { n, bits, def, nodes, root: NIL }
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    /// Number of allocated nodes (a memory indicator).
    pub fn node_count(&self) -> usize {
        self.nodes.len() - 1
    }

    fn agg_of(&self, u: u32, level: u32) -> &S<A> {
        if u == NIL {
            &self.def[level as usize]
        } else {
            &self.nodes[u as usize].agg
        }
    }

    /// A fresh node for an untouched block of `2^level` default elements.
    fn fresh(&mut self, level: u32) -> u32 {
        let agg = self.def[level as usize].clone();
        self.nodes.push(LazyNode { left: NIL, right: NIL, agg, lazy: None });
        (self.nodes.len() - 1) as u32
    }

    fn apply_to(&mut self, u: u32, f: &A::F) {
        let n = &mut self.nodes[u as usize];
        n.agg = A::apply(f, &n.agg);
        n.lazy = Some(match &n.lazy {
            Some(g) => A::compose(f, g),
            None => f.clone(),
        });
    }

    /// Make sure `u` (at `level >= 1`) has both children and has pushed its action into them.
    fn push_down(&mut self, u: u32, level: u32) {
        let (mut l, mut r) = (self.nodes[u as usize].left, self.nodes[u as usize].right);
        if l == NIL {
            l = self.fresh(level - 1);
            self.nodes[u as usize].left = l;
        }
        if r == NIL {
            r = self.fresh(level - 1);
            self.nodes[u as usize].right = r;
        }
        if let Some(f) = self.nodes[u as usize].lazy.take() {
            self.apply_to(l, &f);
            self.apply_to(r, &f);
        }
    }

    fn pull(&mut self, u: u32, level: u32) {
        let (l, r) = (self.nodes[u as usize].left, self.nodes[u as usize].right);
        let agg = <A::M as Monoid>::op(self.agg_of(l, level - 1), self.agg_of(r, level - 1));
        self.nodes[u as usize].agg = agg;
    }

    fn ensure_root(&mut self) -> u32 {
        if self.root == NIL {
            self.root = self.fresh(self.bits);
        }
        self.root
    }

    /// Apply `f` to every element of `[l, r)`.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn apply_range(&mut self, l: usize, r: usize, f: A::F) {
        assert!(l <= r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
        if l == r {
            return;
        }
        let root = self.ensure_root();
        self.apply_rec(root, self.bits, 0, l, r, &f);
    }

    fn apply_rec(&mut self, u: u32, level: u32, lo: usize, l: usize, r: usize, f: &A::F) {
        let hi = lo + (1usize << level);
        if r <= lo || hi <= l {
            return;
        }
        if l <= lo && hi <= r {
            self.apply_to(u, f);
            return;
        }
        self.push_down(u, level);
        let mid = lo + (1usize << (level - 1));
        let (left, right) = (self.nodes[u as usize].left, self.nodes[u as usize].right);
        self.apply_rec(left, level - 1, lo, l, r, f);
        self.apply_rec(right, level - 1, mid, l, r, f);
        self.pull(u, level);
    }

    /// Set element `i` to `x`.
    ///
    /// # Panics
    /// Panics if `i >= len()`.
    pub fn set(&mut self, i: usize, x: S<A>) {
        assert!(i < self.n, "index {i} out of range for length {}", self.n);
        let mut u = self.ensure_root();
        let mut path = Vec::with_capacity(self.bits as usize);
        for level in (1..=self.bits).rev() {
            self.push_down(u, level);
            path.push((u, level));
            let right = (i >> (level - 1)) & 1 == 1;
            u = if right { self.nodes[u as usize].right } else { self.nodes[u as usize].left };
        }
        self.nodes[u as usize].agg = x;
        for &(p, level) in path.iter().rev() {
            self.pull(p, level);
        }
    }

    /// Element `i`.
    pub fn get(&self, i: usize) -> S<A> {
        assert!(i < self.n, "index {i} out of range for length {}", self.n);
        self.prod(i, i + 1)
    }

    /// The fold of `[l, r)` (the identity for an empty range). Does not allocate.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn prod(&self, l: usize, r: usize) -> S<A> {
        assert!(l <= r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
        self.prod_rec(self.root, self.bits, 0, l, r, None)
    }

    /// `pending` is the composed action of the ancestors that has not been pushed into `u`.
    fn prod_rec(&self, u: u32, level: u32, lo: usize, l: usize, r: usize, pending: Option<&A::F>) -> S<A> {
        let hi = lo + (1usize << level);
        if r <= lo || hi <= l {
            return <A::M as Monoid>::identity();
        }
        if l <= lo && hi <= r {
            let agg = self.agg_of(u, level);
            return match pending {
                Some(f) => A::apply(f, agg),
                None => agg.clone(),
            };
        }
        // `u` is partially covered, so it is a real node unless it is NIL (then all below is NIL).
        let (left, right, own) = if u == NIL {
            (NIL, NIL, None)
        } else {
            let node = &self.nodes[u as usize];
            (node.left, node.right, node.lazy.as_ref())
        };
        let composed;
        let below = match (own, pending) {
            (Some(g), Some(f)) => {
                composed = A::compose(f, g);
                Some(&composed)
            }
            (Some(g), None) => Some(g),
            (None, p) => p,
        };
        let mid = lo + (1usize << (level - 1));
        let a = self.prod_rec(left, level - 1, lo, l, r, below);
        let b = self.prod_rec(right, level - 1, mid, l, r, below);
        <A::M as Monoid>::op(&a, &b)
    }
}
