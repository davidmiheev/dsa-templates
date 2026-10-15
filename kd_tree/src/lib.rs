//! A 2D kd-tree over a fixed set of points with a value per point, supporting
//!
//! * [`KdTree::prod_rect`]: fold the values of the points inside a rectangle,
//! * [`KdTree::apply_rect`]: apply a [`MapMonoid`] action to the points inside a rectangle,
//! * [`KdTree::set`] / [`KdTree::get`]: replace / read the value of one point.
//!
//! Small subtrees (at most 32 points) are stored as flat buckets that are scanned linearly.
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

/// Ranges with at most this many points are leaf buckets.
const LEAF: usize = 32;

/// See the [module documentation](self).
pub struct KdTree<A: MapMonoid, X, Y> {
    n: usize,
    /// Points and their values by *tree position*; the subtree of the range `[lo, hi)` is rooted
    /// at `(lo + hi) / 2`. Ranges of at most [`LEAF`] points are buckets.
    pts: Vec<(X, Y)>,
    vals: Vec<S<A>>,
    /// Per-subtree data, stored at the root position of the subtree.
    nodes: Vec<Node<A, X, Y>>,
    /// `pos[i]` = tree position of input point `i`.
    pos: Vec<u32>,
}

struct Node<A: MapMonoid, X, Y> {
    /// Inclusive bounding box `(x_min, x_max, y_min, y_max)` of the subtree.
    bbox: (X, X, Y, Y),
    agg: S<A>,
    /// Pending action for the two children (internal nodes only; the node's own point and `agg`
    /// already reflect it).
    lazy: Option<A::F>,
}

impl<A: MapMonoid, X: Clone, Y: Clone> Clone for Node<A, X, Y> {
    fn clone(&self) -> Self {
        Node { bbox: self.bbox.clone(), agg: self.agg.clone(), lazy: self.lazy.clone() }
    }
}

type Rect<X, Y> = (X, X, Y, Y);

fn inside<X: Ord, Y: Ord>(p: &(X, Y), r: &Rect<X, Y>) -> bool {
    r.0 <= p.0 && p.0 < r.1 && r.2 <= p.1 && p.1 < r.3
}

fn disjoint<X: Ord, Y: Ord>(bb: &Rect<X, Y>, r: &Rect<X, Y>) -> bool {
    bb.1 < r.0 || r.1 <= bb.0 || bb.3 < r.2 || r.3 <= bb.2
}

fn contained<X: Ord, Y: Ord>(bb: &Rect<X, Y>, r: &Rect<X, Y>) -> bool {
    r.0 <= bb.0 && bb.1 < r.1 && r.2 <= bb.2 && bb.3 < r.3
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
        let mut t = KdTree { n, pts: Vec::new(), vals: Vec::new(), nodes: Vec::new(), pos: vec![0; n] };
        if n > 0 {
            // Placeholders that `build` overwrites position by position.
            t.pts = vec![points[0]; n];
            t.vals = vec![values[0].clone(); n];
            let (x0, y0) = points[0];
            t.nodes = vec![Node { bbox: (x0, x0, y0, y0), agg: values[0].clone(), lazy: None }; n];
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
        self.nodes[mid].bbox = (xl, xh, yl, yh);
        if hi - lo <= LEAF {
            for (p, &i) in (lo..hi).zip(&order[lo..hi]) {
                let i = i as usize;
                self.pos[i] = p as u32;
                self.pts[p] = points[i];
                self.vals[p] = values[i].clone();
            }
            self.nodes[mid].agg = self.fold_bucket(lo, hi);
            return;
        }
        // Alternate the split axis by level (subtree sizes halve at every level, so the parity
        // of `log2(size)` flips). Splitting on the larger *spread* instead is tempting, but
        // skewed inputs (one axis with few distinct values) then produce boxes that are thin
        // in one dimension and long in the other, and rectangle queries degrade to `O(n)`.
        if (hi - lo).ilog2() & 1 == 0 {
            slice.select_nth_unstable_by_key(mid - lo, |&i| points[i as usize].0);
        } else {
            slice.select_nth_unstable_by_key(mid - lo, |&i| points[i as usize].1);
        }
        let i = order[mid] as usize;
        self.pos[i] = mid as u32;
        self.pts[mid] = points[i];
        self.vals[mid] = values[i].clone();
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

    /// The fold of the values at positions `lo..hi` (a bucket).
    fn fold_bucket(&self, lo: usize, hi: usize) -> S<A> {
        let mut agg = <A::M as Monoid>::identity();
        for v in &self.vals[lo..hi] {
            agg = <A::M as Monoid>::op(&agg, v);
        }
        agg
    }

    /// Recompute the aggregate of the internal node over `[lo, hi)` from its children.
    fn pull(&mut self, lo: usize, hi: usize) {
        let mid = (lo + hi) / 2;
        let mut agg = self.vals[mid].clone();
        if lo < mid {
            agg = <A::M as Monoid>::op(&self.nodes[(lo + mid) / 2].agg, &agg);
        }
        if mid + 1 < hi {
            agg = <A::M as Monoid>::op(&agg, &self.nodes[(mid + 1 + hi) / 2].agg);
        }
        self.nodes[mid].agg = agg;
    }

    /// Apply `f` to *every* point of the subtree over `[lo, hi)`.
    fn apply_all(&mut self, lo: usize, hi: usize, f: &A::F) {
        let mid = (lo + hi) / 2;
        let agg = A::apply(f, &self.nodes[mid].agg);
        self.nodes[mid].agg = agg;
        if hi - lo <= LEAF {
            for v in &mut self.vals[lo..hi] {
                *v = A::apply(f, v);
            }
        } else {
            self.vals[mid] = A::apply(f, &self.vals[mid]);
            self.nodes[mid].lazy = Some(match self.nodes[mid].lazy.take() {
                Some(g) => A::compose(f, &g),
                None => f.clone(),
            });
        }
    }

    /// Push the pending action of the internal node over `[lo, hi)` into its children.
    fn push_down(&mut self, lo: usize, hi: usize) {
        let mid = (lo + hi) / 2;
        if let Some(f) = self.nodes[mid].lazy.take() {
            if lo < mid {
                self.apply_all(lo, mid, &f);
            }
            if mid + 1 < hi {
                self.apply_all(mid + 1, hi, &f);
            }
        }
    }

    /// Apply `f` to the values of all points with `x_lo <= x < x_hi` and `y_lo <= y < y_hi`.
    pub fn apply_rect(&mut self, x_lo: X, x_hi: X, y_lo: Y, y_hi: Y, f: A::F) {
        self.apply_rec(0, self.n, &(x_lo, x_hi, y_lo, y_hi), &f);
    }

    fn apply_rec(&mut self, lo: usize, hi: usize, r: &Rect<X, Y>, f: &A::F) {
        if lo >= hi {
            return;
        }
        let mid = (lo + hi) / 2;
        if disjoint(&self.nodes[mid].bbox, r) {
            return;
        }
        if contained(&self.nodes[mid].bbox, r) {
            self.apply_all(lo, hi, f);
            return;
        }
        if hi - lo <= LEAF {
            for p in lo..hi {
                if inside(&self.pts[p], r) {
                    self.vals[p] = A::apply(f, &self.vals[p]);
                }
            }
            self.nodes[mid].agg = self.fold_bucket(lo, hi);
            return;
        }
        self.push_down(lo, hi);
        if inside(&self.pts[mid], r) {
            self.vals[mid] = A::apply(f, &self.vals[mid]);
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

    fn prod_rec(&mut self, lo: usize, hi: usize, r: &Rect<X, Y>) -> S<A> {
        let identity = <A::M as Monoid>::identity;
        if lo >= hi {
            return identity();
        }
        let mid = (lo + hi) / 2;
        if disjoint(&self.nodes[mid].bbox, r) {
            return identity();
        }
        if contained(&self.nodes[mid].bbox, r) {
            return self.nodes[mid].agg.clone();
        }
        if hi - lo <= LEAF {
            let mut res = identity();
            for p in lo..hi {
                if inside(&self.pts[p], r) {
                    res = <A::M as Monoid>::op(&res, &self.vals[p]);
                }
            }
            return res;
        }
        self.push_down(lo, hi);
        let mut res = if inside(&self.pts[mid], r) { self.vals[mid].clone() } else { identity() };
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
        if hi - lo <= LEAF {
            self.vals[target] = x;
            self.nodes[mid].agg = self.fold_bucket(lo, hi);
            return;
        }
        self.push_down(lo, hi);
        if target == mid {
            self.vals[mid] = x;
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
        while hi - lo > LEAF {
            let mid = (lo + hi) / 2;
            if target == mid {
                break;
            }
            self.push_down(lo, hi);
            if target < mid {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        self.vals[target].clone()
    }
}
