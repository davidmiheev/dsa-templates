//! Li Chao tree: the minimum of a set of lines (or line segments) at integer query points.
//!
//! The domain is a fixed inclusive integer interval `[lo, hi]` (it may be huge: nodes are created
//! lazily, so memory is proportional to the number of insertions, not to the domain size). Each
//! node keeps the one line that is lowest at the middle of its interval; an inserted line either
//! replaces it and pushes the loser down, or is pushed down itself, towards the single side on
//! which it can still win.
//!
//! * [`LiChaoTree::add_line`]: `O(log C)` for `C = hi - lo + 1`,
//! * [`LiChaoTree::add_segment`]: a line restricted to `[l, r)`, `O(log^2 C)`,
//! * [`LiChaoTree::query`]: `O(log C)`.
//!
//! Queries return the value passed as `empty` to [`LiChaoTree::new`] when no line covers the
//! point, so there is no `Option` to unwrap: use `i64::MAX` (or any value larger than every
//! real answer) as the "infinity" of your problem. For maxima, negate slopes and intercepts.
//!
//! The value type `V` only needs `Copy + Ord + Add + Mul + From<i64>`; use `i128` if
//! `|a * x + b|` can overflow `i64`.
//!
//! # Example
//!
//! ```
//! use li_chao_tree::LiChaoTree;
//!
//! let mut t = LiChaoTree::<i64>::new(-10, 10, i64::MAX);
//! t.add_line(1, 0); //  y = x
//! t.add_line(-1, 4); // y = -x + 4
//! assert_eq!(t.query(0), 0);
//! assert_eq!(t.query(5), -1);
//! t.add_segment(0, 3, 0, -100); // y = -100 only for x in [0, 3)
//! assert_eq!(t.query(2), -100);
//! assert_eq!(t.query(3), 1); // the segment ends before x = 3
//! ```

use std::ops::{Add, Mul};

type Line<V> = (V, V);

struct Node<V> {
    line: Option<Line<V>>,
    left: u32,
    right: u32,
}

/// Li Chao tree for minimum queries over the inclusive integer domain `[lo, hi]`.
pub struct LiChaoTree<V> {
    lo: i64,
    hi: i64,
    empty: V,
    nodes: Vec<Node<V>>,
}

impl<V> LiChaoTree<V>
where
    V: Copy + Ord + Add<Output = V> + Mul<Output = V> + From<i64>,
{
    /// An empty tree over `[lo, hi]`; `empty` is what [`Self::query`] returns when no line
    /// covers the queried point.
    ///
    /// # Panics
    /// Panics if `lo > hi`.
    pub fn new(lo: i64, hi: i64, empty: V) -> Self {
        assert!(lo <= hi, "empty domain [{lo}, {hi}]");
        LiChaoTree { lo, hi, empty, nodes: vec![Node { line: None, left: 0, right: 0 }] }
    }

    #[inline]
    fn eval(line: Line<V>, x: i64) -> V {
        line.0 * V::from(x) + line.1
    }

    fn child(&mut self, node: usize, right: bool) -> usize {
        let c = if right { self.nodes[node].right } else { self.nodes[node].left };
        if c != 0 {
            return c as usize;
        }
        self.nodes.push(Node { line: None, left: 0, right: 0 });
        let id = self.nodes.len() - 1;
        if right {
            self.nodes[node].right = id as u32;
        } else {
            self.nodes[node].left = id as u32;
        }
        id
    }

    /// Insert `line` into the subtree of `node`, which covers `[lo, hi]`.
    fn insert(&mut self, mut node: usize, mut lo: i64, mut hi: i64, mut line: Line<V>) {
        loop {
            let Some(mut cur) = self.nodes[node].line else {
                self.nodes[node].line = Some(line);
                return;
            };
            let mid = lo + (hi - lo) / 2;
            if Self::eval(line, mid) < Self::eval(cur, mid) {
                std::mem::swap(&mut cur, &mut line);
            }
            self.nodes[node].line = Some(cur); // `cur` is now the lower line at `mid`
            if lo == hi {
                return;
            }
            // `line` is the loser at `mid`: it can only win on one side of it.
            if Self::eval(line, lo) < Self::eval(cur, lo) {
                node = self.child(node, false);
                hi = mid;
            } else if Self::eval(line, hi) < Self::eval(cur, hi) {
                node = self.child(node, true);
                lo = mid + 1;
            } else {
                return;
            }
        }
    }

    /// Add the line `y = a * x + b` on the whole domain.
    pub fn add_line(&mut self, a: V, b: V) {
        self.insert(0, self.lo, self.hi, (a, b));
    }

    /// Add the line `y = a * x + b` for `x` in the half-open range `[l, r)` (clamped to the
    /// domain; an empty range adds nothing).
    pub fn add_segment(&mut self, l: i64, r: i64, a: V, b: V) {
        let (l, r) = (l.max(self.lo), r.min(self.hi.saturating_add(1)));
        if l < r {
            self.add_to_range(0, self.lo, self.hi, l, r - 1, (a, b));
        }
    }

    fn add_to_range(&mut self, node: usize, lo: i64, hi: i64, l: i64, r: i64, line: Line<V>) {
        if r < lo || hi < l {
            return;
        }
        if l <= lo && hi <= r {
            self.insert(node, lo, hi, line);
            return;
        }
        // Only descend (and allocate) where the range actually reaches.
        let mid = lo + (hi - lo) / 2;
        if l <= mid {
            let left = self.child(node, false);
            self.add_to_range(left, lo, mid, l, r, line);
        }
        if r > mid {
            let right = self.child(node, true);
            self.add_to_range(right, mid + 1, hi, l, r, line);
        }
    }

    /// The minimum over all lines covering `x`, or the `empty` value if there is none.
    ///
    /// # Panics
    /// Panics if `x` is outside the domain.
    pub fn query(&self, x: i64) -> V {
        assert!(self.lo <= x && x <= self.hi, "x = {x} outside the domain [{}, {}]", self.lo, self.hi);
        let (mut node, mut lo, mut hi) = (0usize, self.lo, self.hi);
        let mut best = self.empty;
        loop {
            if let Some(line) = self.nodes[node].line {
                best = best.min(Self::eval(line, x));
            }
            if lo == hi {
                return best;
            }
            let mid = lo + (hi - lo) / 2;
            let next = if x <= mid {
                hi = mid;
                self.nodes[node].left
            } else {
                lo = mid + 1;
                self.nodes[node].right
            };
            if next == 0 {
                return best;
            }
            node = next as usize;
        }
    }
}
