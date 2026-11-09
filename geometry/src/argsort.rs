//! Sorting points by argument.

use crate::{cross, Point};
use std::cmp::Ordering;

fn class(p: Point) -> u8 {
    match (p.1.signum(), p.0.signum()) {
        (-1, _) => 0,
        (0, -1) => 3,
        (0, _) => 1, // positive x axis and the origin (argument 0)
        _ => 2,
    }
}

/// Sorts the points by the argument `atan2(y, x)` in `(-pi, pi]`, counterclockwise starting just
/// after the negative x axis; the origin counts as argument 0 and equal arguments keep their
/// input order. Exact (integer cross products).
pub fn argument_sort(points: &mut [Point]) {
    points.sort_by(|&a, &b| {
        let (ca, cb) = (class(a), class(b));
        if ca != cb {
            return ca.cmp(&cb);
        }
        if ca == 1 || ca == 3 {
            return Ordering::Equal;
        }
        // same open half plane: a before b iff cross(a, b) > 0
        let c = cross((0, 0), a, b);
        0.cmp(&c)
    });
}
