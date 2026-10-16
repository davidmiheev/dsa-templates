//! Computational geometry on integer points, with exact predicates.
//!
//! * [`convex_hull`], [`farthest_pair`], [`all_farthest_neighbors`] (rotating calipers /
//!   monotone divide and conquer on a convex polygon), [`closest_pair`].
//! * [`argument_sort`]: sort by angle with integer cross products.
//! * [`minimum_enclosing_circle`]: exact incremental algorithm; returns which points lie on the
//!   circle.
//! * [`manhattan_mst`] (four octant sweeps) and [`euclidean_mst`] (Borůvka over a kd-tree).
//! * [`convex_layers`]: onion peeling with a decremental convex hull (Overmars–van Leeuwen).
//! * [`TriangleCounter`]: points of a set strictly inside a triangle of another set, `O(1)` per
//!   query after `O(N^2 M)` preprocessing.
//!
//! A point is a pair `(x, y)` of `i64`.
//!
//! # Example
//!
//! ```
//! use geometry::*;
//!
//! let pts = [(0, 0), (4, 0), (4, 4), (0, 4), (2, 2), (2, 0)];
//! assert_eq!(convex_hull(&pts).len(), 4);
//! assert_eq!(farthest_pair(&pts), (0, 2));
//! ```

mod argsort;
mod closest;
mod enclosing;
mod hull;
mod layers;
mod mst;
mod triangle;

pub use argsort::argument_sort;
pub use closest::closest_pair;
pub use enclosing::minimum_enclosing_circle;
pub use hull::{all_farthest_neighbors, convex_hull, farthest_pair};
pub use layers::convex_layers;
pub use mst::{euclidean_mst, manhattan_mst};
pub use triangle::TriangleCounter;

/// A point of the plane.
pub type Point = (i64, i64);

/// The cross product `(b - a) x (c - a)`, exact for coordinates up to `2^61`.
#[inline]
pub(crate) fn cross(a: Point, b: Point, c: Point) -> i128 {
    (b.0 - a.0) as i128 * (c.1 - a.1) as i128 - (b.1 - a.1) as i128 * (c.0 - a.0) as i128
}

/// The squared distance, exact for coordinates up to `2^61`.
#[inline]
pub(crate) fn dist2(a: Point, b: Point) -> u128 {
    let (dx, dy) = ((a.0 - b.0) as i128, (a.1 - b.1) as i128);
    (dx * dx + dy * dy) as u128
}
