//! Implicit treap: a dynamic sequence with `O(log n)` expected time per operation.
//!
//! Elements are keyed by their *position*, not by a stored key: every node knows the size of its
//! subtree, and `split` / `merge` cut and join the sequence at positions. On top of that:
//!
//! * [`ImplicitTreap::insert`] / [`ImplicitTreap::remove`] at any position,
//! * [`ImplicitTreap::reverse`] of a range (a lazy "swap children" flag),
//! * [`ImplicitTreap::apply`] of a lazy action from a [`MapMonoid`] to a range, and
//! * [`ImplicitTreap::prod`], the fold of a monoid over a range.
//!
//! **Reversal and non-commutative monoids.** `reverse` flips the order of the elements but keeps
//! each node's stored aggregate, which is only correct when the aggregate does not depend on the
//! order, i.e. when the monoid is commutative (sums, min, max, xor, counts). Do not combine
//! `reverse` with order-sensitive aggregates such as function composition.
//!
//! The nodes live in one `Vec` (index 0 is the null node), so there is no per-node allocation
//! and no pointer chasing through `Box`es.
//!
//! # Example: range add, range sum, reversal
//!
//! ```
//! use algebra::{MapMonoid, Monoid};
//! use implicit_treap::ImplicitTreap;
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
//! let mut t = ImplicitTreap::<Add>::from_slice(&[(1, 1), (2, 1), (3, 1)]);
//! t.insert(1, (10, 1)); // [1, 10, 2, 3]
//! t.reverse(0, 3); // [2, 10, 1, 3]
//! t.apply(0, 2, 5); // [7, 15, 1, 3]
//! assert_eq!(t.prod(0, 4).0, 26);
//! assert_eq!(t.get(1).0, 15);
//! assert_eq!(t.remove(0).0, 7);
//! assert_eq!(t.len(), 3);
//! ```

use algebra::{MapMonoid, Monoid};

type S<A> = <<A as MapMonoid>::M as Monoid>::S;

struct Node<A: MapMonoid> {
    left: u32,
    right: u32,
    prio: u32,
    size: u32,
    val: S<A>,
    /// Fold of the whole subtree (in sequence order).
    agg: S<A>,
    /// Pending action for the *children* of this node (already applied to `val` and `agg`).
    lazy: A::F,
    has_lazy: bool,
    /// Pending reversal of the subtree below this node.
    rev: bool,
}

/// A sequence of monoid elements with lazy range actions and reversal.
pub struct ImplicitTreap<A: MapMonoid> {
    /// `nodes[0]` is the null node (size 0, aggregate = identity); it is never modified.
    nodes: Vec<Node<A>>,
    root: u32,
    free: Vec<u32>,
    seed: u64,
}

impl<A: MapMonoid> Default for ImplicitTreap<A> {
    fn default() -> Self {
        Self::new()
    }
}

impl<A: MapMonoid> ImplicitTreap<A> {
    /// An empty sequence.
    pub fn new() -> Self {
        let null = Node {
            left: 0,
            right: 0,
            prio: 0,
            size: 0,
            val: <A::M as Monoid>::identity(),
            agg: <A::M as Monoid>::identity(),
            lazy: A::id_map(),
            has_lazy: false,
            rev: false,
        };
        ImplicitTreap { nodes: vec![null], root: 0, free: Vec::new(), seed: 0x9E37_79B9_7F4A_7C15 }
    }

    /// Build from a slice in `O(n)` (a Cartesian tree over random priorities).
    pub fn from_slice(v: &[S<A>]) -> Self {
        let mut t = Self::new();
        t.nodes.reserve(v.len());
        let mut stack: Vec<u32> = Vec::new();
        for x in v {
            let id = t.new_node(x.clone());
            let mut last = 0;
            while let Some(&top) = stack.last() {
                if t.nodes[top as usize].prio < t.nodes[id as usize].prio {
                    last = stack.pop().unwrap();
                } else {
                    break;
                }
            }
            t.nodes[id as usize].left = last;
            if let Some(&top) = stack.last() {
                t.nodes[top as usize].right = id;
            }
            stack.push(id);
        }
        t.root = stack.first().copied().unwrap_or(0);
        t.pull_all(t.root);
        t
    }

    fn random(&mut self) -> u32 {
        self.seed ^= self.seed << 13;
        self.seed ^= self.seed >> 7;
        self.seed ^= self.seed << 17;
        (self.seed >> 32) as u32
    }

    fn new_node(&mut self, val: S<A>) -> u32 {
        let prio = self.random();
        let node = Node {
            left: 0,
            right: 0,
            prio,
            size: 1,
            agg: val.clone(),
            val,
            lazy: A::id_map(),
            has_lazy: false,
            rev: false,
        };
        if let Some(id) = self.free.pop() {
            self.nodes[id as usize] = node;
            id
        } else {
            self.nodes.push(node);
            (self.nodes.len() - 1) as u32
        }
    }

    /// Recompute sizes and aggregates of a freshly built subtree (children first).
    fn pull_all(&mut self, t: u32) {
        if t == 0 {
            return;
        }
        let (l, r) = (self.nodes[t as usize].left, self.nodes[t as usize].right);
        self.pull_all(l);
        self.pull_all(r);
        self.pull(t);
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.nodes[self.root as usize].size as usize
    }

    pub fn is_empty(&self) -> bool {
        self.root == 0
    }

    fn size(&self, t: u32) -> usize {
        self.nodes[t as usize].size as usize
    }

    fn pull(&mut self, t: u32) {
        let (l, r) = (self.nodes[t as usize].left as usize, self.nodes[t as usize].right as usize);
        let size = 1 + self.nodes[l].size + self.nodes[r].size;
        let left_and_self = <A::M as Monoid>::op(&self.nodes[l].agg, &self.nodes[t as usize].val);
        let agg = <A::M as Monoid>::op(&left_and_self, &self.nodes[r].agg);
        let n = &mut self.nodes[t as usize];
        n.size = size;
        n.agg = agg;
    }

    fn apply_to(&mut self, t: u32, f: &A::F) {
        if t == 0 {
            return;
        }
        let n = &mut self.nodes[t as usize];
        n.val = A::apply(f, &n.val);
        n.agg = A::apply(f, &n.agg);
        n.lazy = if n.has_lazy { A::compose(f, &n.lazy) } else { f.clone() };
        n.has_lazy = true;
    }

    fn push(&mut self, t: u32) {
        let (rev, has_lazy) = (self.nodes[t as usize].rev, self.nodes[t as usize].has_lazy);
        if !rev && !has_lazy {
            return;
        }
        if rev {
            let n = &mut self.nodes[t as usize];
            std::mem::swap(&mut n.left, &mut n.right);
            n.rev = false;
            let (l, r) = (n.left as usize, n.right as usize);
            // The null node is never flagged (flags on it would be shared by every empty child).
            if l != 0 {
                self.nodes[l].rev ^= true;
            }
            if r != 0 {
                self.nodes[r].rev ^= true;
            }
        }
        if has_lazy {
            let f = self.nodes[t as usize].lazy.clone();
            let n = &mut self.nodes[t as usize];
            n.has_lazy = false;
            n.lazy = A::id_map();
            let (l, r) = (n.left, n.right);
            self.apply_to(l, &f);
            self.apply_to(r, &f);
        }
    }

    /// Split into (first `k` elements, the rest).
    fn split(&mut self, t: u32, k: usize) -> (u32, u32) {
        if t == 0 {
            return (0, 0);
        }
        self.push(t);
        let left = self.nodes[t as usize].left;
        if self.size(left) >= k {
            let (a, b) = self.split(left, k);
            self.nodes[t as usize].left = b;
            self.pull(t);
            (a, t)
        } else {
            let right = self.nodes[t as usize].right;
            let (a, b) = self.split(right, k - self.size(left) - 1);
            self.nodes[t as usize].right = a;
            self.pull(t);
            (t, b)
        }
    }

    fn merge(&mut self, a: u32, b: u32) -> u32 {
        if a == 0 {
            return b;
        }
        if b == 0 {
            return a;
        }
        if self.nodes[a as usize].prio > self.nodes[b as usize].prio {
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

    fn check_range(&self, l: usize, r: usize) {
        assert!(l <= r && r <= self.len(), "invalid range [{l}, {r}) for length {}", self.len());
    }

    /// Run `f` on the node of the range `[l, r)` (0 for an empty range), then stitch back.
    fn with_range<R>(&mut self, l: usize, r: usize, f: impl FnOnce(&mut Self, u32) -> R) -> R {
        self.check_range(l, r);
        let (a, rest) = self.split(self.root, l);
        let (m, b) = self.split(rest, r - l);
        let out = f(self, m);
        let ab = self.merge(a, m);
        self.root = self.merge(ab, b);
        out
    }

    /// Insert `x` so that it becomes the element at position `i`.
    ///
    /// # Panics
    /// Panics if `i > len()`.
    pub fn insert(&mut self, i: usize, x: S<A>) {
        assert!(i <= self.len(), "insert position {i} out of range for length {}", self.len());
        let node = self.new_node(x);
        let (a, b) = self.split(self.root, i);
        let an = self.merge(a, node);
        self.root = self.merge(an, b);
    }

    /// Remove and return the element at position `i`.
    ///
    /// # Panics
    /// Panics if `i >= len()`.
    pub fn remove(&mut self, i: usize) -> S<A> {
        assert!(i < self.len(), "index {i} out of range for length {}", self.len());
        let (a, rest) = self.split(self.root, i);
        let (m, b) = self.split(rest, 1);
        self.root = self.merge(a, b);
        let val = self.nodes[m as usize].val.clone();
        self.free.push(m);
        val
    }

    /// The element at position `i`.
    ///
    /// # Panics
    /// Panics if `i >= len()`.
    pub fn get(&mut self, mut i: usize) -> S<A> {
        assert!(i < self.len(), "index {i} out of range for length {}", self.len());
        let mut t = self.root;
        loop {
            self.push(t);
            let left = self.nodes[t as usize].left;
            let ls = self.size(left);
            if i < ls {
                t = left;
            } else if i == ls {
                return self.nodes[t as usize].val.clone();
            } else {
                i -= ls + 1;
                t = self.nodes[t as usize].right;
            }
        }
    }

    /// Reverse the elements of `[l, r)`. See the module docs: the monoid must be commutative.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn reverse(&mut self, l: usize, r: usize) {
        self.with_range(l, r, |t, m| {
            if m != 0 {
                t.nodes[m as usize].rev ^= true;
            }
        });
    }

    /// Apply the action `f` to every element of `[l, r)`.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn apply(&mut self, l: usize, r: usize, f: A::F) {
        self.with_range(l, r, |t, m| t.apply_to(m, &f));
    }

    /// The fold of `[l, r)` (the identity for an empty range).
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn prod(&mut self, l: usize, r: usize) -> S<A> {
        self.with_range(l, r, |t, m| t.nodes[m as usize].agg.clone())
    }

    /// All elements in order.
    pub fn to_vec(&mut self) -> Vec<S<A>> {
        let mut out = Vec::with_capacity(self.len());
        let mut stack: Vec<u32> = Vec::new();
        let mut t = self.root;
        loop {
            while t != 0 {
                self.push(t);
                stack.push(t);
                t = self.nodes[t as usize].left;
            }
            match stack.pop() {
                None => return out,
                Some(n) => {
                    out.push(self.nodes[n as usize].val.clone());
                    t = self.nodes[n as usize].right;
                }
            }
        }
    }
}
