//! Convex layers (onion peeling) with a decremental convex hull.

use crate::Point;

const NONE: i32 = -1;

#[derive(Clone, Copy)]
struct Node {
    /// the bridge between the hulls of the two children: the last left-hull point of the left
    /// child and the first of the right child (both are point indices)
    bl: i32,
    br: i32,
    /// the point range `[lo, hi)` of the subtree
    lo: i32,
    hi: i32,
    left: i32,
    right: i32,
}

/// The left hull of a set of points sorted by `(y, x)`, maintained under deletions with a
/// segment tree whose nodes store the bridge between the hulls of their halves (Overmars and
/// van Leeuwen); a deletion costs `O(log^2 n)`.
struct LeftHull {
    ps: Vec<Point>,
    nodes: Vec<Node>,
    root: i32,
}

impl LeftHull {
    fn new(ps: Vec<Point>) -> Self {
        let n = ps.len();
        let mut h = LeftHull { ps, nodes: vec![Node { bl: 0, br: 0, lo: 0, hi: 0, left: NONE, right: NONE }; 2 * n], root: 0 };
        h.build(0, 0, n as i32);
        h
    }

    /// Exact in `i64` for coordinates up to `2^24` in absolute value.
    fn cross(&self, a: i32, b: i32, c: i32) -> i64 {
        let (pa, pb, pc) = (self.ps[a as usize], self.ps[b as usize], self.ps[c as usize]);
        (pb.0 - pa.0) * (pc.1 - pa.1) - (pb.1 - pa.1) * (pc.0 - pa.0)
    }

    fn is_leaf(&self, w: i32) -> bool {
        self.nodes[w as usize].left == NONE && self.nodes[w as usize].right == NONE
    }

    /// Recomputes the bridge of node `w` from its two children.
    fn pull(&mut self, w: i32) {
        let (mut l, mut r) = (self.nodes[w as usize].left, self.nodes[w as usize].right);
        let split_y = self.ps[self.nodes[r as usize].lo as usize].1 as i128;
        while !self.is_leaf(l) || !self.is_leaf(r) {
            let (a, b) = (self.nodes[l as usize].bl, self.nodes[l as usize].br);
            let (c, d) = (self.nodes[r as usize].bl, self.nodes[r as usize].br);
            if a != b && self.cross(a, b, c) > 0 {
                l = self.nodes[l as usize].left;
            } else if c != d && self.cross(b, c, d) > 0 {
                r = self.nodes[r as usize].right;
            } else if a == b {
                r = self.nodes[r as usize].left;
            } else if c == d {
                l = self.nodes[l as usize].right;
            } else {
                let s1 = self.cross(a, b, c) as i128;
                let s2 = self.cross(b, a, d) as i128;
                debug_assert!(s1 + s2 >= 0);
                let (yd, yc) = (self.ps[d as usize].1 as i128, self.ps[c as usize].1 as i128);
                if s1 + s2 == 0 || s1 * yd + s2 * yc < split_y * (s1 + s2) {
                    l = self.nodes[l as usize].right;
                } else {
                    r = self.nodes[r as usize].left;
                }
            }
        }
        self.nodes[w as usize].bl = self.nodes[l as usize].lo;
        self.nodes[w as usize].br = self.nodes[r as usize].lo;
    }

    fn build(&mut self, w: i32, lo: i32, hi: i32) {
        let node = &mut self.nodes[w as usize];
        node.lo = lo;
        node.hi = hi;
        if hi - lo == 1 {
            node.left = NONE;
            node.right = NONE;
            node.bl = lo;
            node.br = lo;
        } else {
            let mid = (lo + hi) / 2;
            let (l, r) = (w + 1, w + 2 * (mid - lo));
            node.left = l;
            node.right = r;
            self.build(l, lo, mid);
            self.build(r, mid, hi);
            self.pull(w);
        }
    }

    /// Deletes the points with indices in `[lo, hi)` from the subtree of `w`; returns the new
    /// root of that subtree (`NONE` if it became empty).
    fn erase(&mut self, w: i32, lo: i32, hi: i32) -> i32 {
        if w == NONE {
            return NONE;
        }
        let node = self.nodes[w as usize];
        if hi <= node.lo || lo >= node.hi {
            return w;
        }
        if lo <= node.lo && hi >= node.hi {
            return NONE;
        }
        let l = self.erase(node.left, lo, hi);
        let r = self.erase(node.right, lo, hi);
        self.nodes[w as usize].left = l;
        self.nodes[w as usize].right = r;
        if l == NONE {
            return r;
        }
        if r == NONE {
            return l;
        }
        self.pull(w);
        w
    }

    /// Appends the hull points of the subtree `w` that lie between the point indices `l` and `r`.
    fn collect(&self, w: i32, l: i32, r: i32, out: &mut Vec<i32>) {
        let node = &self.nodes[w as usize];
        if self.is_leaf(w) {
            out.push(node.lo);
        } else if r <= node.bl {
            self.collect(node.left, l, r, out);
        } else if l >= node.br {
            self.collect(node.right, l, r, out);
        } else {
            self.collect(node.left, l, node.bl, out);
            self.collect(node.right, node.br, r, out);
        }
    }

    fn hull(&self) -> Vec<i32> {
        let mut out = Vec::new();
        if self.root != NONE {
            self.collect(self.root, 0, self.ps.len() as i32 - 1, &mut out);
        }
        out
    }

    fn remove(&mut self, i: i32) {
        self.root = self.erase(self.root, i, i + 1);
    }
}

/// For each of the distinct points, the number of the iteration (starting at 1) in which it is
/// removed when the convex hull boundary (vertices and the points on its edges) of the remaining
/// points is removed repeatedly; `O(n log^2 n)`.
///
/// # Panics
/// Panics unless the points are distinct with coordinates of absolute value at most `2^24`.
pub fn convex_layers(points: &[Point]) -> Vec<usize> {
    let n = points.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| (points[i].1, points[i].0));
    assert!(order.windows(2).all(|w| points[w[0]] != points[w[1]]), "the points must be distinct");
    assert!(points.iter().all(|&(x, y)| x.abs() <= 1 << 24 && y.abs() <= 1 << 24), "coordinates must be at most 2^24");
    let sorted: Vec<Point> = order.iter().map(|&i| points[i]).collect();
    let mut left = LeftHull::new(sorted.clone());
    // the right hull is the left hull of the point reflection, listed in reverse order
    let reflected: Vec<Point> = sorted.iter().rev().map(|&(x, y)| (-x, -y)).collect();
    let mut right = LeftHull::new(reflected);
    let mut layer = vec![0usize; n];
    let (mut removed, mut round) = (0usize, 0usize);
    while removed < n {
        round += 1;
        let mut hull: Vec<i32> = left.hull();
        hull.extend(right.hull().into_iter().map(|i| (n as i32 - 1) - i));
        hull.sort_unstable();
        hull.dedup();
        for &i in &hull {
            layer[order[i as usize]] = round;
            removed += 1;
            left.remove(i);
            right.remove(n as i32 - 1 - i);
        }
    }
    layer
}
