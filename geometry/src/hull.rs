//! Convex hull and farthest pairs.

use crate::{cross, dist2, Point};

/// The vertices of the convex hull in counterclockwise order, as indices into `points`: no three
/// consecutive vertices are collinear and duplicates are reported once (by their first index).
/// One index for a single distinct point, two for a segment. `O(n log n)`.
pub fn convex_hull(points: &[Point]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..points.len()).collect();
    order.sort_by_key(|&i| (points[i], i));
    order.dedup_by_key(|&mut i| points[i]);
    let n = order.len();
    if n <= 2 {
        return order;
    }
    let mut hull: Vec<usize> = Vec::with_capacity(n + 1);
    for pass in 0..2 {
        let start = hull.len();
        let iter: Box<dyn Iterator<Item = &usize>> = if pass == 0 { Box::new(order.iter()) } else { Box::new(order.iter().rev()) };
        for &i in iter {
            while hull.len() >= start + 2 && cross(points[hull[hull.len() - 2]], points[hull[hull.len() - 1]], points[i]) <= 0 {
                hull.pop();
            }
            hull.push(i);
        }
        hull.pop();
    }
    hull
}

/// A pair of distinct indices whose points are as far apart as possible (any such pair), by
/// rotating calipers on the hull, `O(n log n)`.
///
/// # Panics
/// Panics if there are fewer than two points.
pub fn farthest_pair(points: &[Point]) -> (usize, usize) {
    assert!(points.len() >= 2, "need at least two points");
    let hull = convex_hull(points);
    let h = hull.len();
    if h == 1 {
        // all points coincide
        return (0, 1);
    }
    if h == 2 {
        return (hull[0], hull[1]);
    }
    let p = |i: usize| points[hull[i % h]];
    let (mut best, mut pair) = (0u128, (hull[0], hull[1]));
    let mut j = 1usize;
    for i in 0..h {
        // advance j while the next vertex is farther from the edge (i, i + 1)
        while cross(p(i), p(i + 1), p(j + 1)).abs() > cross(p(i), p(i + 1), p(j)).abs() {
            j += 1;
        }
        for k in [i, i + 1] {
            let d = dist2(p(k), p(j));
            if d > best {
                best = d;
                pair = (hull[k % h], hull[j % h]);
            }
        }
    }
    pair
}

/// For every vertex of a strictly convex polygon (counterclockwise), the index of a vertex
/// farthest from it, in `O(n log n)`.
///
/// The distance matrix of the doubled vertex list satisfies the quadrangle inequality (the
/// diagonals of a convex quadrilateral are longer than two opposite sides), so the leftmost
/// maximum of row `i` over the columns `i + 1 ..= i + n - 1` moves right with `i` and a divide
/// and conquer finds all of them.
///
/// # Panics
/// Panics if the polygon has fewer than 2 vertices.
pub fn all_farthest_neighbors(polygon: &[Point]) -> Vec<usize> {
    let n = polygon.len();
    assert!(n >= 2, "need at least two vertices");
    let at = |j: usize| polygon[j % n];
    let mut opt = vec![0usize; n];
    // rows lo..hi, answers within columns clo..=chi
    let mut stack = vec![(0usize, n, 1usize, 2 * n - 2)];
    while let Some((lo, hi, clo, chi)) = stack.pop() {
        if lo >= hi {
            continue;
        }
        let mid = (lo + hi) / 2;
        let from = clo.max(mid + 1);
        let to = chi.min(mid + n - 1);
        let (mut best, mut best_j) = (0u128, from);
        for j in from..=to {
            let d = dist2(at(mid), at(j));
            if d > best {
                best = d;
                best_j = j;
            }
        }
        opt[mid] = best_j;
        stack.push((lo, mid, clo, best_j));
        stack.push((mid + 1, hi, best_j, chi));
    }
    opt.iter().map(|&j| j % n).collect()
}
