//! Stress tests of every geometry routine against brute force on small integer coordinates.
#![allow(clippy::needless_range_loop, clippy::collapsible_if, clippy::neg_multiply)]

use geometry::*;
use std::collections::HashSet;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn cross(a: Point, b: Point, c: Point) -> i64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}

fn d2(a: Point, b: Point) -> i64 {
    (a.0 - b.0).pow(2) + (a.1 - b.1).pow(2)
}

fn random_points(rng: &mut Rng, n: usize, range: i64) -> Vec<Point> {
    (0..n).map(|_| (rng.range(-range, range + 1), rng.range(-range, range + 1))).collect()
}

/// The hull vertex test: `p` is a vertex iff it is not in the closed convex hull of the others.
fn is_vertex(pts: &[Point], i: usize) -> bool {
    let p = pts[i];
    let others: Vec<Point> = pts.iter().enumerate().filter(|&(j, &q)| j != i && q != p).map(|(_, &q)| q).collect();
    let on_segment = |a: Point, b: Point| cross(a, b, p) == 0 && (a.0.min(b.0)..=a.0.max(b.0)).contains(&p.0) && (a.1.min(b.1)..=a.1.max(b.1)).contains(&p.1);
    for &a in &others {
        for &b in &others {
            if on_segment(a, b) {
                return false;
            }
            for &c in &others {
                let (s1, s2, s3) = (cross(a, b, p), cross(b, c, p), cross(c, a, p));
                if (s1 >= 0 && s2 >= 0 && s3 >= 0) || (s1 <= 0 && s2 <= 0 && s3 <= 0) {
                    if cross(a, b, c) != 0 {
                        return false;
                    }
                }
            }
        }
    }
    true
}

#[test]
fn hull_pairs_and_sorting() {
    let mut st = Stress::new("geometry/basic");
    for case in 0..4000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = 2 + rng.below(if case % 20 == 0 { 25 } else { 9 });
        let range = [1i64, 2, 5, 20][rng.below(4)];
        let pts = random_points(&mut rng, n, range);
        // hull
        let hull = convex_hull(&pts);
        let want: HashSet<Point> = (0..n).filter(|&i| is_vertex(&pts, i)).map(|i| pts[i]).collect();
        let got: HashSet<Point> = hull.iter().map(|&i| pts[i]).collect();
        assert_eq!(got, want, "hull of {pts:?}");
        assert_eq!(hull.len(), got.len(), "duplicates in the hull");
        if hull.len() >= 3 {
            for k in 0..hull.len() {
                assert!(cross(pts[hull[k]], pts[hull[(k + 1) % hull.len()]], pts[hull[(k + 2) % hull.len()]]) > 0, "not strictly ccw");
            }
        }
        // farthest and closest pairs
        let mut far = 0;
        let mut near = i64::MAX;
        for i in 0..n {
            for j in i + 1..n {
                far = far.max(d2(pts[i], pts[j]));
                near = near.min(d2(pts[i], pts[j]));
            }
        }
        let (a, b) = farthest_pair(&pts);
        assert!(a != b && d2(pts[a], pts[b]) == far, "farthest of {pts:?}: ({a},{b})");
        let (a, b) = closest_pair(&pts);
        assert!(a != b && d2(pts[a], pts[b]) == near, "closest of {pts:?}: ({a},{b})");
        // argument sort against atan2 (exact enough for tiny coordinates)
        let mut sorted = pts.clone();
        argument_sort(&mut sorted);
        let key = |p: Point| -> f64 {
            if p == (0, 0) {
                0.0
            } else {
                (p.1 as f64).atan2(p.0 as f64)
            }
        };
        assert!(sorted.windows(2).all(|w| key(w[0]) <= key(w[1]) + 1e-12), "argument order of {pts:?}: {sorted:?}");
        let mut a1 = pts.clone();
        let mut a2 = sorted.clone();
        a1.sort_unstable();
        a2.sort_unstable();
        assert_eq!(a1, a2);
        st.case();
        st.ops(n * n);
    }
    // all farthest neighbours of strictly convex polygons
    for case in 0..2000 * scale() {
        let mut rng = Rng::new(50_000 + case as u64);
        let cnt = 3 + rng.below(40);
        let pts = random_points(&mut rng, cnt, 30);
        let hull = convex_hull(&pts);
        if hull.len() < 3 {
            continue;
        }
        let poly: Vec<Point> = hull.iter().map(|&i| pts[i]).collect();
        let f = all_farthest_neighbors(&poly);
        for i in 0..poly.len() {
            let best = (0..poly.len()).map(|j| d2(poly[i], poly[j])).max().unwrap();
            assert_eq!(d2(poly[i], poly[f[i]]), best, "polygon {poly:?} vertex {i}");
        }
    }
    st.done();
}

#[test]
fn enclosing_circle_and_triangles() {
    let mut st = Stress::new("geometry/circle");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(60_000 + case as u64);
        let n = 2 + rng.below(8);
        let mut pts = random_points(&mut rng, n, 15);
        pts.sort_unstable();
        pts.dedup();
        if pts.len() < 2 {
            continue;
        }
        let got = minimum_enclosing_circle(&pts);
        // brute force in f64: the smallest circle through two or three points containing all
        let mut best_r = f64::INFINITY;
        let mut best: Option<(f64, f64, f64)> = None;
        let inside = |c: (f64, f64, f64)| pts.iter().all(|&p| ((p.0 as f64 - c.0).powi(2) + (p.1 as f64 - c.1).powi(2)).sqrt() <= c.2 + 1e-9);
        for i in 0..pts.len() {
            for j in i + 1..pts.len() {
                let c = ((pts[i].0 + pts[j].0) as f64 / 2.0, (pts[i].1 + pts[j].1) as f64 / 2.0, (d2(pts[i], pts[j]) as f64).sqrt() / 2.0);
                if c.2 < best_r && inside(c) {
                    best_r = c.2;
                    best = Some(c);
                }
                for k in j + 1..pts.len() {
                    let (a, b, cc) = (pts[i], pts[j], pts[k]);
                    let det = 2.0 * cross(a, b, cc) as f64;
                    if det == 0.0 {
                        continue;
                    }
                    let (a2, b2, c2) = (a.0 * a.0 + a.1 * a.1, b.0 * b.0 + b.1 * b.1, cc.0 * cc.0 + cc.1 * cc.1);
                    let ux = (a2 * (b.1 - cc.1) + b2 * (cc.1 - a.1) + c2 * (a.1 - b.1)) as f64 / (-det) * -1.0;
                    let uy = (a2 * (cc.0 - b.0) + b2 * (a.0 - cc.0) + c2 * (b.0 - a.0)) as f64 / (-det) * -1.0;
                    let c = (ux, uy, ((a.0 as f64 - ux).powi(2) + (a.1 as f64 - uy).powi(2)).sqrt());
                    if c.2 < best_r - 1e-12 && inside(c) {
                        best_r = c.2;
                        best = Some(c);
                    }
                }
            }
        }
        let c = best.unwrap();
        for (idx, &p) in pts.iter().enumerate() {
            let on = (((p.0 as f64 - c.0).powi(2) + (p.1 as f64 - c.1).powi(2)).sqrt() - c.2).abs() < 1e-7;
            assert_eq!(got[idx], on, "points {pts:?}: circle {c:?} point {p:?}");
        }
        st.case();
    }
    // points in triangles
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(70_000 + case as u64);
        let (n, m) = (1 + rng.below(6), 1 + rng.below(12));
        let range = [2i64, 4, 8][rng.below(3)];
        let a = random_points(&mut rng, n, range);
        let b = random_points(&mut rng, m, range);
        let counter = TriangleCounter::new(&a, &b);
        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    let (pa, pb, pc) = (a[x], a[y], a[z]);
                    let orient = cross(pa, pb, pc);
                    let want = if orient == 0 {
                        0
                    } else {
                        b.iter()
                            .filter(|&&p| {
                                let (s1, s2, s3) = (cross(pa, pb, p), cross(pb, pc, p), cross(pc, pa, p));
                                (s1 > 0 && s2 > 0 && s3 > 0) || (s1 < 0 && s2 < 0 && s3 < 0)
                            })
                            .count() as u32
                    };
                    assert_eq!(counter.count(x, y, z), want, "A={a:?} B={b:?} triangle {x} {y} {z}");
                }
            }
        }
        st.case();
    }
    st.done();
}

fn prim_squared(pts: &[Point]) -> Vec<i64> {
    let n = pts.len();
    let (mut in_tree, mut dist) = (vec![false; n], vec![i64::MAX; n]);
    dist[0] = 0;
    let mut out = Vec::new();
    for _ in 0..n {
        let u = (0..n).filter(|&i| !in_tree[i]).min_by_key(|&i| dist[i]).unwrap();
        in_tree[u] = true;
        out.push(dist[u]);
        for v in 0..n {
            if !in_tree[v] {
                dist[v] = dist[v].min(d2(pts[u], pts[v]));
            }
        }
    }
    out.remove(0);
    out.sort_unstable();
    out
}

#[test]
fn spanning_trees_and_layers() {
    let mut st = Stress::new("geometry/mst");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(80_000 + case as u64);
        let n = 1 + rng.below(if case % 10 == 0 { 60 } else { 12 });
        let rg = [3i64, 10, 100][rng.below(3)];
        let pts = random_points(&mut rng, n, rg);
        // Euclidean
        let tree = euclidean_mst(&pts);
        assert_eq!(tree.len(), n - 1);
        let mut w: Vec<i64> = tree.iter().map(|&(u, v)| d2(pts[u], pts[v])).collect();
        w.sort_unstable();
        assert_eq!(w, if n == 1 { Vec::new() } else { prim_squared(&pts) }, "points {pts:?}");
        // Manhattan
        let man = |a: Point, b: Point| (a.0 - b.0).abs() + (a.1 - b.1).abs();
        let (total, tree) = manhattan_mst(&pts);
        assert_eq!(tree.len(), n - 1);
        let (mut in_tree, mut dist) = (vec![false; n], vec![i64::MAX; n]);
        dist[0] = 0;
        let mut want = 0i64;
        for _ in 0..n {
            let u = (0..n).filter(|&i| !in_tree[i]).min_by_key(|&i| dist[i]).unwrap();
            in_tree[u] = true;
            want += dist[u];
            for v in 0..n {
                if !in_tree[v] {
                    dist[v] = dist[v].min(man(pts[u], pts[v]));
                }
            }
        }
        assert_eq!(total as i64, want, "manhattan of {pts:?}");
        assert_eq!(tree.iter().map(|&(u, v)| man(pts[u], pts[v])).sum::<i64>(), want);
        // spanning: n - 1 edges connecting everything
        let mut comp: Vec<usize> = (0..n).collect();
        for &(u, v) in &tree {
            let (cu, cv) = (comp[u], comp[v]);
            assert_ne!(cu, cv, "cycle");
            for c in comp.iter_mut() {
                if *c == cv {
                    *c = cu;
                }
            }
        }
        st.case();
    }
    // convex layers against naive peeling
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(90_000 + case as u64);
        let range = [2i64, 3, 5, 12][rng.below(4)];
        let cnt = 1 + rng.below(if case % 10 == 0 { 40 } else { 12 });
        let mut pts = random_points(&mut rng, cnt, range);
        pts.sort_unstable();
        pts.dedup();
        let n = pts.len();
        // naive: repeatedly remove all points on the boundary (vertices and edge points)
        let mut layer = vec![0usize; n];
        let mut left = n;
        let mut round = 0;
        while left > 0 {
            round += 1;
            let alive: Vec<usize> = (0..n).filter(|&i| layer[i] == 0).collect();
            let mut hull_pts: Vec<usize> = Vec::new();
            for &i in &alive {
                // on the boundary iff some line through it has all alive points on one closed side
                let boundary = alive.iter().all(|_| true) && {
                    let mut ok = false;
                    // candidate directions: towards every other alive point (and the 4 axis directions)
                    let mut dirs: Vec<Point> = alive.iter().filter(|&&j| j != i).map(|&j| (pts[j].0 - pts[i].0, pts[j].1 - pts[i].1)).collect();
                    dirs.extend([(1, 0), (0, 1), (-1, 0), (0, -1)]);
                    for &(dx, dy) in &dirs {
                        let line = |j: usize| dx * (pts[j].1 - pts[i].1) - dy * (pts[j].0 - pts[i].0);
                        if alive.iter().all(|&j| line(j) >= 0) || alive.iter().all(|&j| line(j) <= 0) {
                            ok = true;
                            break;
                        }
                    }
                    ok
                };
                if boundary {
                    hull_pts.push(i);
                }
            }
            for i in hull_pts {
                layer[i] = round;
                left -= 1;
            }
        }
        assert_eq!(convex_layers(&pts), layer, "points {pts:?}");
        st.case();
    }
    st.done();
}
