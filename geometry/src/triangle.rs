//! Counting points strictly inside triangles.

use crate::Point;

/// Answers "how many points of `B` lie strictly inside the triangle `A_a A_b A_c`" in `O(1)`.
///
/// For every pair of points of `A` it counts the points of `B` strictly between their
/// horizontal lines that lie strictly left of / on the segment, plus the points of `B` on the
/// horizontal line of every point of `A` left of / below it; a triangle is then the signed
/// combination of its three edges (the classic "ray to the left" decomposition, which treats
/// boundary points consistently). Preprocessing `O(|A|^2 |B|)`.
pub struct TriangleCounter {
    a: Vec<Point>,
    n: usize,
    /// `edge_left[u * n + v]`: points of `B` with `a[u].y < y < a[v].y` strictly left of the line
    /// `u -> v` ... (0 unless `a[u]` is below `a[v]`)
    edge_left: Vec<u32>,
    /// the same points exactly on the segment's line
    edge_on: Vec<u32>,
    /// points of `B` on the horizontal line through `a[u]` strictly left of it
    point_left: Vec<u32>,
    /// points of `B` equal to `a[u]`
    point_on: Vec<u32>,
}

fn below(p: Point, q: Point) -> bool {
    (p.1, p.0) < (q.1, q.0)
}

impl TriangleCounter {
    pub fn new(a: &[Point], b: &[Point]) -> Self {
        let n = a.len();
        let mut edge_left = vec![0u32; n * n];
        let mut edge_on = vec![0u32; n * n];
        let mut point_left = vec![0u32; n];
        let mut point_on = vec![0u32; n];
        for (u, &pu) in a.iter().enumerate() {
            for &q in b {
                if q.1 == pu.1 {
                    if q.0 < pu.0 {
                        point_left[u] += 1;
                    }
                    if q.0 == pu.0 {
                        point_on[u] += 1;
                    }
                }
            }
        }
        for u in 0..n {
            for v in 0..n {
                if a[u].1 >= a[v].1 {
                    continue;
                }
                let (mut left, mut on) = (0, 0);
                for &q in b {
                    if q.1 > a[u].1 && q.1 < a[v].1 {
                        let det = (a[u].0 - a[v].0) as i128 * (q.1 - a[v].1) as i128 - (a[u].1 - a[v].1) as i128 * (q.0 - a[v].0) as i128;
                        if det < 0 {
                            left += 1;
                        } else if det == 0 {
                            on += 1;
                        }
                    }
                }
                edge_left[u * n + v] = left;
                edge_on[u * n + v] = on;
            }
        }
        TriangleCounter { a: a.to_vec(), n, edge_left, edge_on, point_left, point_on }
    }

    /// The number of points of `B` strictly inside the triangle `A_a A_b A_c` (0 if it is
    /// degenerate).
    pub fn count(&self, a: usize, b: usize, c: usize) -> u32 {
        let n = self.n;
        let (mut a, mut b, mut c) = (a, b, c);
        if below(self.a[b], self.a[a]) {
            std::mem::swap(&mut a, &mut b);
        }
        if below(self.a[c], self.a[b]) {
            std::mem::swap(&mut b, &mut c);
        }
        if below(self.a[b], self.a[a]) {
            std::mem::swap(&mut a, &mut b);
        }
        let (pa, pb, pc) = (self.a[a], self.a[b], self.a[c]);
        let det = (pa.0 - pc.0) as i128 * (pb.1 - pc.1) as i128 - (pa.1 - pc.1) as i128 * (pb.0 - pc.0) as i128;
        if det == 0 {
            return 0;
        }
        let l = |u: usize, v: usize| self.edge_left[u * n + v] as i64;
        let m = |u: usize, v: usize| self.edge_on[u * n + v] as i64;
        let ans: i64 = if pa.1 == pb.1 {
            // c on top, a and b below it
            l(b, c) - (l(a, c) + m(a, c))
        } else if pb.1 == pc.1 {
            l(a, c) - (l(a, b) + m(a, b))
        } else if det < 0 {
            l(a, c) - (l(b, c) + m(b, c)) - (l(a, b) + m(a, b)) - (self.point_left[b] + self.point_on[b]) as i64
        } else {
            l(a, b) + l(b, c) + self.point_left[b] as i64 - (l(a, c) + m(a, c))
        };
        ans as u32
    }
}
