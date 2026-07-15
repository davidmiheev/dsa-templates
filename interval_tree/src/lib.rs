//! Interval tree for overlap queries on half-open `[lo, hi)` intervals.
//!
//! An `Interval<T>` is a tagged pair `(lo, hi)` with an associated
//! payload. This crate builds the simple static interval tree described by
//! the CLRS-style construction:
//!
//! * median-split on `mid = (lo + hi) / 2`
//! * store intervals whose midpoint lies on the left of `mid` in `left`
//! * store intervals whose midpoint lies on the right of `mid` in `right`
//! * store intervals that contain `mid` in `middle` (their `lo <= mid <= hi`)
//!
//! [`IntervalTree::query`] returns every interval whose `[lo, hi)` overlaps
//! a query interval `[qlo, qhi)`.
//!
//! The implementation is the canonical `O(n log n)` build / `O(k log n)`
//! query variant; it is built once on construction and immutable
//! thereafter.

/// A half-open interval `[lo, hi)` carrying a payload.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval<T> {
    /// Inclusive lower bound.
    pub lo: i64,
    /// Exclusive upper bound.
    pub hi: i64,
    /// Associated payload.
    pub value: T,
}

impl<T> Interval<T> {
    /// Construct a new interval.
    pub fn new(lo: i64, hi: i64, value: T) -> Self {
        debug_assert!(lo <= hi, "interval must satisfy lo <= hi");
        Self { lo, hi, value }
    }

    /// Return the midpoint `(lo + hi) / 2`.
    fn mid(&self) -> i64 {
        (self.lo + self.hi) / 2
    }

    /// Does this interval overlap `[qlo, qhi)`?
    fn overlaps(&self, qlo: i64, qhi: i64) -> bool {
        self.lo < qhi && qlo < self.hi
    }

    fn contains_point(&self, p: i64) -> bool {
        self.lo <= p && p <= self.hi
    }
}

/// Static interval tree.
#[derive(Debug)]
pub struct IntervalTree<T> {
    root: Option<Node<T>>,
}

/// Internal node.
#[derive(Debug)]
struct Node<T> {
    mid: i64,
    middle: Vec<Interval<T>>,
    left: Option<Box<Node<T>>>,
    right: Option<Box<Node<T>>>,
}

impl<T: Clone> IntervalTree<T> {
    /// Build an interval tree from a slice of intervals.
    ///
    /// # Panics
    ///
    /// Panics if any interval has `lo > hi`.
    pub fn new(intervals: &[Interval<T>]) -> Self {
        Self { root: build(intervals, 0, intervals.len()) }
    }

    /// Return all intervals whose `[lo, hi)` overlaps `[qlo, qhi)`.
    pub fn query(&self, qlo: i64, qhi: i64) -> Vec<Interval<T>> {
        let mut out = Vec::new();
        if let Some(root) = self.root.as_ref() {
            query_node(root, qlo, qhi, &mut out);
        }
        out
    }
}

fn build<T: Clone>(intervals: &[Interval<T>], lo: usize, hi: usize) -> Option<Node<T>> {
    if lo >= hi {
        return None;
    }
    // Split point of the slice
    let mid_idx = (lo + hi) / 2;
    let mid_point = intervals[mid_idx].mid();

    let mut middle = Vec::new();
    let mut left_intervals = Vec::new();
    let mut right_intervals = Vec::new();
    for it in &intervals[lo..hi] {
        if it.contains_point(mid_point) {
            middle.push(it.clone());
        } else if it.hi < mid_point {
            left_intervals.push(it.clone());
        } else if it.lo > mid_point {
            right_intervals.push(it.clone());
        } else {
            // Shouldn't happen: interval ends on one side, starts on other, but doesn't contain mid.
            // Fall back to the midpoint comparison.
            if it.mid() < mid_point {
                left_intervals.push(it.clone());
            } else {
                right_intervals.push(it.clone());
            }
        }
    }

    Some(Node {
        mid: mid_point,
        middle,
        left: build(&left_intervals, 0, left_intervals.len()).map(Box::new),
        right: build(&right_intervals, 0, right_intervals.len()).map(Box::new),
    })
}

fn query_node<T: Clone>(node: &Node<T>, qlo: i64, qhi: i64, out: &mut Vec<Interval<T>>) {
    // Check middle (intervals that contain the node's `mid`).
    for it in &node.middle {
        if it.overlaps(qlo, qhi) {
            out.push(it.clone());
        }
    }
    if qlo < node.mid {
        if let Some(left) = node.left.as_deref() {
            query_node(left, qlo, qhi, out);
        }
    }
    if qhi > node.mid {
        if let Some(right) = node.right.as_deref() {
            query_node(right, qlo, qhi, out);
        }
    }
}
