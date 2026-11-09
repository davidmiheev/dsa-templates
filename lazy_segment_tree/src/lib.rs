//! Lazy segment tree: range update, range query.
//!
//! Elements form a [`Monoid`](algebra::Monoid) `M`; updates are actions from a
//! [`MapMonoid`](algebra::MapMonoid) `A` (range add, range assign, range affine, ...). Each
//! operation takes `O(log n)`:
//!
//! * [`LazySegmentTree::apply_range`]: apply an action to every element of `[l, r)`,
//! * [`LazySegmentTree::prod`]: combine the elements of `[l, r)`,
//! * [`LazySegmentTree::max_right`] / [`LazySegmentTree::min_left`]: binary search on the
//!   tree for the longest prefix / suffix whose aggregate satisfies a monotone predicate.
//!
//! Ranges are half-open `[l, r)`. The layout is the usual power-of-two bottom-up one, with the
//! pending action of each inner node stored in `lazy[node]`.
//!
//! # Example: range add, range minimum
//!
//! ```
//! use algebra::{MapMonoid, Monoid};
//! use lazy_segment_tree::LazySegmentTree;
//!
//! enum Min {}
//! impl Monoid for Min {
//!     type S = i64;
//!     fn identity() -> i64 { i64::MAX }
//!     fn op(a: &i64, b: &i64) -> i64 { *a.min(b) }
//! }
//! enum AddMin {}
//! impl MapMonoid for AddMin {
//!     type M = Min;
//!     type F = i64;
//!     fn id_map() -> i64 { 0 }
//!     fn apply(f: &i64, x: &i64) -> i64 { if *x == i64::MAX { *x } else { x + f } }
//!     fn compose(f: &i64, g: &i64) -> i64 { f + g }
//! }
//!
//! let mut t = LazySegmentTree::<AddMin>::from_slice(&[5, 1, 4, 2, 8]);
//! assert_eq!(t.prod(0, 5), 1);
//! t.apply_range(1, 3, 10); // [5, 11, 14, 2, 8]
//! assert_eq!(t.prod(0, 3), 5);
//! assert_eq!(t.prod(1, 3), 11);
//! ```

use algebra::{MapMonoid, Monoid};

type S<A> = <<A as MapMonoid>::M as Monoid>::S;

/// Lazy segment tree over the monoid `A::M` with range actions `A::F`.
pub struct LazySegmentTree<A: MapMonoid> {
    n: usize,
    size: usize,
    log: usize,
    d: Vec<S<A>>,
    lazy: Vec<A::F>,
}

impl<A: MapMonoid> LazySegmentTree<A> {
    /// A tree of `n` identity elements.
    pub fn new(n: usize) -> Self {
        Self::from_slice(&vec![<A::M as Monoid>::identity(); n])
    }

    /// Build from a slice in `O(n)`.
    pub fn from_slice(v: &[S<A>]) -> Self {
        let n = v.len();
        let size = n.next_power_of_two();
        let log = size.trailing_zeros() as usize;
        let mut d = vec![<A::M as Monoid>::identity(); 2 * size];
        d[size..size + n].clone_from_slice(v);
        let mut t = LazySegmentTree { n, size, log, d, lazy: vec![A::id_map(); size] };
        for k in (1..size).rev() {
            t.pull(k);
        }
        t
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    #[inline]
    fn pull(&mut self, k: usize) {
        self.d[k] = <A::M as Monoid>::op(&self.d[2 * k], &self.d[2 * k + 1]);
    }

    #[inline]
    fn all_apply(&mut self, k: usize, f: &A::F) {
        self.d[k] = A::apply(f, &self.d[k]);
        if k < self.size {
            self.lazy[k] = A::compose(f, &self.lazy[k]);
        }
    }

    #[inline]
    fn push(&mut self, k: usize) {
        let f = self.lazy[k].clone();
        self.all_apply(2 * k, &f);
        self.all_apply(2 * k + 1, &f);
        self.lazy[k] = A::id_map();
    }

    /// Push pending actions down along the boundary paths of leaves `l` and `r` (heap indices).
    #[inline]
    fn push_boundaries(&mut self, l: usize, r: usize) {
        for i in (1..=self.log).rev() {
            if ((l >> i) << i) != l {
                self.push(l >> i);
            }
            if ((r >> i) << i) != r {
                self.push((r - 1) >> i);
            }
        }
    }

    /// Replace element `p` by `x`.
    ///
    /// # Panics
    /// Panics if `p >= len()`.
    pub fn set(&mut self, p: usize, x: S<A>) {
        assert!(p < self.n, "index {p} out of range for length {}", self.n);
        let p = p + self.size;
        for i in (1..=self.log).rev() {
            self.push(p >> i);
        }
        self.d[p] = x;
        for i in 1..=self.log {
            self.pull(p >> i);
        }
    }

    /// The current value of element `p`.
    ///
    /// # Panics
    /// Panics if `p >= len()`.
    pub fn get(&mut self, p: usize) -> S<A> {
        assert!(p < self.n, "index {p} out of range for length {}", self.n);
        let p = p + self.size;
        for i in (1..=self.log).rev() {
            self.push(p >> i);
        }
        self.d[p].clone()
    }

    /// `op` over the elements of `[l, r)` (the identity if the range is empty).
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn prod(&mut self, l: usize, r: usize) -> S<A> {
        assert!(l <= r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
        if l == r {
            return <A::M as Monoid>::identity();
        }
        let (mut l, mut r) = (l + self.size, r + self.size);
        self.push_boundaries(l, r);
        let (mut left, mut right) = (<A::M as Monoid>::identity(), <A::M as Monoid>::identity());
        while l < r {
            if l & 1 == 1 {
                left = <A::M as Monoid>::op(&left, &self.d[l]);
                l += 1;
            }
            if r & 1 == 1 {
                r -= 1;
                right = <A::M as Monoid>::op(&self.d[r], &right);
            }
            l >>= 1;
            r >>= 1;
        }
        <A::M as Monoid>::op(&left, &right)
    }

    /// `op` over all elements.
    pub fn all_prod(&self) -> S<A> {
        self.d[1].clone()
    }

    /// Apply the action `f` to element `p`.
    ///
    /// # Panics
    /// Panics if `p >= len()`.
    pub fn apply(&mut self, p: usize, f: A::F) {
        assert!(p < self.n, "index {p} out of range for length {}", self.n);
        let p = p + self.size;
        for i in (1..=self.log).rev() {
            self.push(p >> i);
        }
        self.d[p] = A::apply(&f, &self.d[p]);
        for i in 1..=self.log {
            self.pull(p >> i);
        }
    }

    /// Apply the action `f` to every element of `[l, r)`.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn apply_range(&mut self, l: usize, r: usize, f: A::F) {
        assert!(l <= r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
        if l == r {
            return;
        }
        let (l, r) = (l + self.size, r + self.size);
        self.push_boundaries(l, r);
        let (mut a, mut b) = (l, r);
        while a < b {
            if a & 1 == 1 {
                self.all_apply(a, &f);
                a += 1;
            }
            if b & 1 == 1 {
                b -= 1;
                self.all_apply(b, &f);
            }
            a >>= 1;
            b >>= 1;
        }
        for i in 1..=self.log {
            if ((l >> i) << i) != l {
                self.pull(l >> i);
            }
            if ((r >> i) << i) != r {
                self.pull((r - 1) >> i);
            }
        }
    }

    /// The largest `r >= l` such that `pred(op(a[l..r]))` holds, for a predicate that is
    /// monotone (true on a prefix of growing ranges) with `pred(identity())` true.
    ///
    /// # Panics
    /// Panics if `l > len()`, or if `pred(identity())` is false.
    pub fn max_right<G: Fn(&S<A>) -> bool>(&mut self, l: usize, pred: G) -> usize {
        assert!(l <= self.n, "index {l} out of range for length {}", self.n);
        assert!(pred(&<A::M as Monoid>::identity()), "pred(identity) must hold");
        if l == self.n {
            return self.n;
        }
        let mut l = l + self.size;
        for i in (1..=self.log).rev() {
            self.push(l >> i);
        }
        let mut acc = <A::M as Monoid>::identity();
        loop {
            while l % 2 == 0 {
                l >>= 1;
            }
            if !pred(&<A::M as Monoid>::op(&acc, &self.d[l])) {
                while l < self.size {
                    self.push(l);
                    l *= 2;
                    let next = <A::M as Monoid>::op(&acc, &self.d[l]);
                    if pred(&next) {
                        acc = next;
                        l += 1;
                    }
                }
                return l - self.size;
            }
            acc = <A::M as Monoid>::op(&acc, &self.d[l]);
            l += 1;
            if l & l.wrapping_neg() == l {
                return self.n;
            }
        }
    }

    /// The smallest `l <= r` such that `pred(op(a[l..r]))` holds (mirror image of
    /// [`Self::max_right`]).
    ///
    /// # Panics
    /// Panics if `r > len()`, or if `pred(identity())` is false.
    pub fn min_left<G: Fn(&S<A>) -> bool>(&mut self, r: usize, pred: G) -> usize {
        assert!(r <= self.n, "index {r} out of range for length {}", self.n);
        assert!(pred(&<A::M as Monoid>::identity()), "pred(identity) must hold");
        if r == 0 {
            return 0;
        }
        let mut r = r + self.size;
        for i in (1..=self.log).rev() {
            self.push((r - 1) >> i);
        }
        let mut acc = <A::M as Monoid>::identity();
        loop {
            r -= 1;
            while r > 1 && r % 2 == 1 {
                r >>= 1;
            }
            if !pred(&<A::M as Monoid>::op(&self.d[r], &acc)) {
                while r < self.size {
                    self.push(r);
                    r = 2 * r + 1;
                    let next = <A::M as Monoid>::op(&self.d[r], &acc);
                    if pred(&next) {
                        acc = next;
                        r -= 1;
                    }
                }
                return r + 1 - self.size;
            }
            acc = <A::M as Monoid>::op(&self.d[r], &acc);
            if r & r.wrapping_neg() == r {
                return 0;
            }
        }
    }
}
