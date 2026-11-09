//! Closest pair of points.

use crate::Point;

/// A pair of distinct indices with the smallest distance, by divide and conquer, `O(n log n)`.
///
/// # Panics
/// Panics if there are fewer than two points.
pub fn closest_pair(points: &[Point]) -> (usize, usize) {
    let n = points.len();
    assert!(n >= 2, "need at least two points");
    let mut pts: Vec<(i64, i64, usize)> = points.iter().enumerate().map(|(i, &(x, y))| (x, y, i)).collect();
    pts.sort_unstable();
    let mut best = (u128::MAX, 0usize, 0usize);
    let mut scratch = pts.clone();
    rec(&mut pts, &mut scratch, &mut best);
    (best.1, best.2)
}

fn d2(a: (i64, i64, usize), b: (i64, i64, usize)) -> u128 {
    let (dx, dy) = ((a.0 - b.0) as i128, (a.1 - b.1) as i128);
    (dx * dx + dy * dy) as u128
}

/// Sorts `pts` by `y` while updating `best` with the closest pair inside it. On entry `pts` is
/// sorted by `x`.
fn rec(pts: &mut [(i64, i64, usize)], scratch: &mut [(i64, i64, usize)], best: &mut (u128, usize, usize)) {
    let n = pts.len();
    if n <= 3 {
        for i in 0..n {
            for j in i + 1..n {
                let d = d2(pts[i], pts[j]);
                if d < best.0 {
                    *best = (d, pts[i].2, pts[j].2);
                }
            }
        }
        pts.sort_unstable_by_key(|p| (p.1, p.0));
        return;
    }
    let mid = n / 2;
    let split_x = pts[mid].0;
    {
        let (l, r) = pts.split_at_mut(mid);
        let (sl, sr) = scratch.split_at_mut(mid);
        rec(l, sl, best);
        rec(r, sr, best);
    }
    // merge by y
    let (mut i, mut j, mut k) = (0, mid, 0);
    while i < mid && j < n {
        if (pts[i].1, pts[i].0) <= (pts[j].1, pts[j].0) {
            scratch[k] = pts[i];
            i += 1;
        } else {
            scratch[k] = pts[j];
            j += 1;
        }
        k += 1;
    }
    while i < mid {
        scratch[k] = pts[i];
        i += 1;
        k += 1;
    }
    while j < n {
        scratch[k] = pts[j];
        j += 1;
        k += 1;
    }
    pts.copy_from_slice(&scratch[..n]);
    // strip around the dividing line
    let mut strip: Vec<(i64, i64, usize)> = Vec::new();
    for &p in pts.iter() {
        let dx = (p.0 - split_x) as i128;
        if (dx * dx) as u128 >= best.0 {
            continue;
        }
        for &q in strip.iter().rev() {
            let dy = (p.1 - q.1) as i128;
            if (dy * dy) as u128 >= best.0 {
                break;
            }
            let d = d2(p, q);
            if d < best.0 {
                *best = (d, p.2, q.2);
            }
        }
        strip.push(p);
    }
}
