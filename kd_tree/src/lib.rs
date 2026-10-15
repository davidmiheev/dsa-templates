//! A 2D kd-tree over a fixed set of points with a value per point, supporting
//!
//! * [`KdTree::prod_rect`]: fold the values of the points inside a rectangle,
//! * [`KdTree::apply_rect`]: apply a [`MapMonoid`] action to the points inside a rectangle,
//! * [`KdTree::set`] / [`KdTree::get`]: replace / read the value of one point.
//!
//! The *shape* is fixed at construction (points cannot be added), but a point can start
//! "inactive" by carrying the monoid's identity — pair that with a count in the monoid and a
//! point appears later with `set`. Rectangle operations cost `O(sqrt n)`.
//!
//! The fold order inside a rectangle is unspecified, so the monoid must be **commutative**.
//! Rectangles are half-open: `x_lo <= x < x_hi` and `y_lo <= y < y_hi`; coordinates are any
//! `Ord + Copy` type.
//!
//! # Example
//!
//! ```
//! use algebra::{MapMonoid, Monoid};
//! use kd_tree::KdTree;
//!
//! enum Sum {}
//! impl Monoid for Sum {
//!     type S = (i64, i64); // (sum of weights, number of points)
//!     fn identity() -> (i64, i64) { (0, 0) }
//!     fn op(a: &(i64, i64), b: &(i64, i64)) -> (i64, i64) { (a.0 + b.0, a.1 + b.1) }
//! }
//! enum Affine {}
//! impl MapMonoid for Affine {
//!     type M = Sum;
//!     type F = (i64, i64); // w -> a*w + b
//!     fn id_map() -> (i64, i64) { (1, 0) }
//!     fn apply(f: &(i64, i64), x: &(i64, i64)) -> (i64, i64) { (f.0 * x.0 + f.1 * x.1, x.1) }
//!     fn compose(f: &(i64, i64), g: &(i64, i64)) -> (i64, i64) { (f.0 * g.0, f.0 * g.1 + f.1) }
//! }
//!
//! let points = [(0, 0), (5, 5), (9, 1)];
//! let mut t = KdTree::<Affine, i32, i32>::new(&points, &[(1, 1), (10, 1), (100, 1)]);
//! assert_eq!(t.prod_rect(0, 6, 0, 6).0, 11);
//! t.apply_rect(4, 10, 0, 10, (2, 1)); // the points (5, 5) and (9, 1): w -> 2w + 1
//! assert_eq!(t.prod_rect(0, 10, 0, 10).0, 1 + 21 + 201);
//! t.set(0, (7, 1));
//! assert_eq!(t.get(0), (7, 1));
//! ```

use algebra::{MapMonoid, Monoid};

type S<A> = <<A as MapMonoid>::M as Monoid>::S;

/// See the [module documentation](self).
pub struct KdTree<A: MapMonoid, X, Y> {
    n: usize,
    /// Indexed by *tree position*; the subtree of the range `[lo, hi)` is rooted at
    /// `(lo + hi) / 2`. A node keeps all its fields together: a visit is one cache miss.
    nodes: Vec<Node<A, X, Y>>,
    /// `pos[i]` = tree position of input point `i`.
    pos: Vec<u32>,
}

struct Node<A: MapMonoid, X, Y> {
    px: X,
    py: Y,
    /// Inclusive bounding box `(x_min, x_max, y_min, y_max)` of the subtree.
    bbox: (X, X, Y, Y),
    val: S<A>,
    agg: S<A>,
    /// Pending action for the children (already reflected in `val` and `agg` of the node).
    lazy: Option<A::F>,
}

impl<A: MapMonoid, X: Clone, Y: Clone> Clone for Node<A, X, Y> {
    fn clone(&self) -> Self {
        Node { px: self.px.clone(), py: self.py.clone(), bbox: self.bbox.clone(), val: self.val.clone(), agg: self.agg.clone(), lazy: self.lazy.clone() }
    }
}

impl<A: MapMonoid, X: Ord + Copy, Y: Ord + Copy> KdTree<A, X, Y> {
    /// Build over `points` with initial `values` (`values[i]` belongs to `points[i]`).
    ///
    /// # Panics
    /// Panics if the slices differ in length or there are `2^32` or more points.
    pub fn new(points: &[(X, Y)], values: &[S<A>]) -> Self {
        assert_eq!(points.len(), values.len(), "one value per point");
        let n = points.len();
        assert!(n < u32::MAX as usize, "too many points");
        let mut order: Vec<u32> = (0..n as u32).collect();
        let mut t = KdTree { n, nodes: Vec::new(), pos: vec![0; n] };
        if n > 0 {
            // Placeholders that `build` overwrites position by position.
            let (x0, y0) = points[0];
            let blank = Node { px: x0, py: y0, bbox: (x0, x0, y0, y0), val: values[0].clone(), agg: values[0].clone(), lazy: None };
            t.nodes = vec![blank; n];
            t.build(points, values, &mut order, 0, n);
        }
        t
    }

    fn build(&mut self, points: &[(X, Y)], values: &[S<A>], order: &mut [u32], lo: usize, hi: usize) {
        if lo >= hi {
            return;
        }
        let slice = &mut order[lo..hi];
        let (mut xl, mut xh) = (points[slice[0] as usize].0, points[slice[0] as usize].0);
        let (mut yl, mut yh) = (points[slice[0] as usize].1, points[slice[0] as usize].1);
        for &i in slice.iter() {
            let (x, y) = points[i as usize];
            xl = xl.min(x);
            xh = xh.max(x);
            yl = yl.min(y);
            yh = yh.max(y);
        }
        let mid = (lo + hi) / 2;
        // Alternate the split axis by level (subtree sizes halve at every level, so the parity
        // of `log2(size)` flips). Splitting on the larger *spread* instead is tempting, but
        // skewed inputs (one axis with few distinct values) then produce boxes that are thin
        // in one dimension and long in the other, and rectangle queries degrade to `O(n)`.
        let split_x = (hi - lo).ilog2() % 2 == 0;
        if split_x {
            slice.select_nth_unstable_by_key(mid - lo, |&i| points[i as usize].0);
        } else {
            slice.select_nth_unstable_by_key(mid - lo, |&i| points[i as usize].1);
        }
        let i = order[mid] as usize;
        self.pos[i] = mid as u32;
        let node = &mut self.nodes[mid];
        node.px = points[i].0;
        node.py = points[i].1;
        node.val = values[i].clone();
        node.bbox = (xl, xh, yl, yh);
        self.build(points, values, order, lo, mid);
        self.build(points, values, order, mid + 1, hi);
        self.pull(lo, hi);
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Recompute the aggregate of the subtree of `[lo, hi)` from its children.
    fn pull(&mut self, lo: usize, hi: usize) {
        let mid = (lo + hi) / 2;
        let mut agg = self.nodes[mid].val.clone();
        if lo < mid {
            agg = <A::M as Monoid>::op(&self.nodes[(lo + mid) / 2].agg, &agg);
        }
        if mid + 1 < hi {
            agg = <A::M as Monoid>::op(&agg, &self.nodes[(mid + 1 + hi) / 2].agg);
        }
        self.nodes[mid].agg = agg;
    }

    fn apply_to(&mut self, node: usize, f: &A::F) {
        let node = &mut self.nodes[node];
        node.val = A::apply(f, &node.val);
        node.agg = A::apply(f, &node.agg);
        node.lazy = Some(match node.lazy.take() {
            Some(g) => A::compose(f, &g),
            None => f.clone(),
        });
    }

    fn push_down(&mut self, lo: usize, hi: usize) {
        let mid = (lo + hi) / 2;
        if let Some(f) = self.nodes[mid].lazy.take() {
            if lo < mid {
                self.apply_to((lo + mid) / 2, &f);
            }
            if mid + 1 < hi {
                self.apply_to((mid + 1 + hi) / 2, &f);
            }
        }
    }

    fn inside(&self, node: usize, r: &(X, X, Y, Y)) -> bool {
        let n = &self.nodes[node];
        r.0 <= n.px && n.px < r.1 && r.2 <= n.py && n.py < r.3
    }

    fn disjoint(bb: &(X, X, Y, Y), r: &(X, X, Y, Y)) -> bool {
        bb.1 < r.0 || r.1 <= bb.0 || bb.3 < r.2 || r.3 <= bb.2
    }

    fn contained(bb: &(X, X, Y, Y), r: &(X, X, Y, Y)) -> bool {
        r.0 <= bb.0 && bb.1 < r.1 && r.2 <= bb.2 && bb.3 < r.3
    }

    /// Apply `f` to the values of all points with `x_lo <= x < x_hi` and `y_lo <= y < y_hi`.
    pub fn apply_rect(&mut self, x_lo: X, x_hi: X, y_lo: Y, y_hi: Y, f: A::F) {
        self.apply_rec(0, self.n, &(x_lo, x_hi, y_lo, y_hi), &f);
    }

    fn apply_rec(&mut self, lo: usize, hi: usize, r: &(X, X, Y, Y), f: &A::F) {
        if lo >= hi {
            return;
        }
        let mid = (lo + hi) / 2;
        if Self::disjoint(&self.nodes[mid].bbox, r) {
            return;
        }
        if Self::contained(&self.nodes[mid].bbox, r) {
            self.apply_to(mid, f);
            return;
        }
        self.push_down(lo, hi);
        if self.inside(mid, r) {
            self.nodes[mid].val = A::apply(f, &self.nodes[mid].val);
        }
        self.apply_rec(lo, mid, r, f);
        self.apply_rec(mid + 1, hi, r, f);
        self.pull(lo, hi);
    }

    /// The fold of the values of all points with `x_lo <= x < x_hi` and `y_lo <= y < y_hi`
    /// (the identity if there are none). The monoid must be commutative.
    pub fn prod_rect(&mut self, x_lo: X, x_hi: X, y_lo: Y, y_hi: Y) -> S<A> {
        self.prod_rec(0, self.n, &(x_lo, x_hi, y_lo, y_hi))
    }

    fn prod_rec(&mut self, lo: usize, hi: usize, r: &(X, X, Y, Y)) -> S<A> {
        if lo >= hi {
            return <A::M as Monoid>::identity();
        }
        let mid = (lo + hi) / 2;
        if Self::disjoint(&self.nodes[mid].bbox, r) {
            return <A::M as Monoid>::identity();
        }
        if Self::contained(&self.nodes[mid].bbox, r) {
            return self.nodes[mid].agg.clone();
        }
        self.push_down(lo, hi);
        let mut res = if self.inside(mid, r) { self.nodes[mid].val.clone() } else { <A::M as Monoid>::identity() };
        let left = self.prod_rec(lo, mid, r);
        let right = self.prod_rec(mid + 1, hi, r);
        res = <A::M as Monoid>::op(&res, &left);
        <A::M as Monoid>::op(&res, &right)
    }

    /// Replace the value of input point `i` by `x`.
    ///
    /// # Panics
    /// Panics if `i >= len()`.
    pub fn set(&mut self, i: usize, x: S<A>) {
        assert!(i < self.n, "point {i} out of range for {} points", self.n);
        self.set_rec(0, self.n, self.pos[i] as usize, x);
    }

    fn set_rec(&mut self, lo: usize, hi: usize, target: usize, x: S<A>) {
        let mid = (lo + hi) / 2;
        self.push_down(lo, hi);
        if target == mid {
            self.nodes[mid].val = x;
        } else if target < mid {
            self.set_rec(lo, mid, target, x);
        } else {
            self.set_rec(mid + 1, hi, target, x);
        }
        self.pull(lo, hi);
    }

    /// The value of input point `i`.
    ///
    /// # Panics
    /// Panics if `i >= len()`.
    pub fn get(&mut self, i: usize) -> S<A> {
        assert!(i < self.n, "point {i} out of range for {} points", self.n);
        let target = self.pos[i] as usize;
        let (mut lo, mut hi) = (0, self.n);
        loop {
            let mid = (lo + hi) / 2;
            if target == mid {
                return self.nodes[mid].val.clone();
            }
            self.push_down(lo, hi);
            if target < mid {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
    }
}
