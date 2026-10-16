//! Minimum enclosing circle with exact predicates.

use crate::{cross, Point};

#[derive(Clone, Copy)]
enum Circle {
    One(Point),
    Diameter(Point, Point),
    Three(Point, Point, Point),
}

/// The sign of `det` for the in-circle test of `p` against the circle through `a, b, c`
/// (positive when `p` is inside for a counterclockwise triangle).
fn in_circle(a: Point, b: Point, c: Point, p: Point) -> i128 {
    let (ax, ay) = ((a.0 - p.0) as i128, (a.1 - p.1) as i128);
    let (bx, by) = ((b.0 - p.0) as i128, (b.1 - p.1) as i128);
    let (cx, cy) = ((c.0 - p.0) as i128, (c.1 - p.1) as i128);
    let (a2, b2, c2) = (ax * ax + ay * ay, bx * bx + by * by, cx * cx + cy * cy);
    ax * (by * c2 - b2 * cy) - ay * (bx * c2 - b2 * cx) + a2 * (bx * cy - by * cx)
}

impl Circle {
    /// `> 0` strictly inside, `0` on the circle, `< 0` outside.
    fn side(&self, p: Point) -> i128 {
        match *self {
            Circle::One(a) => {
                if p == a {
                    0
                } else {
                    -1
                }
            }
            Circle::Diameter(a, b) => {
                // inside the circle with diameter ab iff the angle apb is obtuse
                -((a.0 - p.0) as i128 * (b.0 - p.0) as i128 + (a.1 - p.1) as i128 * (b.1 - p.1) as i128)
            }
            Circle::Three(a, b, c) => {
                let orient = cross(a, b, c);
                assert!(orient != 0, "degenerate circle");
                let d = in_circle(a, b, c, p);
                if orient > 0 {
                    d
                } else {
                    -d
                }
            }
        }
    }
}

/// For each point, whether it lies on the circumference of the minimum enclosing circle of the
/// set. Randomized incremental construction (expected `O(n)`) with exact integer predicates
/// (coordinates up to about `10^9` in absolute value).
///
/// # Panics
/// Panics on an empty input.
pub fn minimum_enclosing_circle(points: &[Point]) -> Vec<bool> {
    assert!(!points.is_empty(), "need at least one point");
    let n = points.len();
    // a deterministic shuffle (xorshift)
    let mut order: Vec<usize> = (0..n).collect();
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    for i in (1..n).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        order.swap(i, (state % (i as u64 + 1)) as usize);
    }
    let pts: Vec<Point> = order.iter().map(|&i| points[i]).collect();
    let mut circle = Circle::One(pts[0]);
    for i in 1..n {
        if circle.side(pts[i]) >= 0 {
            continue;
        }
        circle = Circle::One(pts[i]);
        for j in 0..i {
            if circle.side(pts[j]) >= 0 {
                continue;
            }
            circle = Circle::Diameter(pts[i], pts[j]);
            for k in 0..j {
                if circle.side(pts[k]) >= 0 {
                    continue;
                }
                circle = Circle::Three(pts[i], pts[j], pts[k]);
            }
        }
    }
    points.iter().map(|&p| circle.side(p) == 0).collect()
}
