//! Minimum spanning trees of point sets under the Manhattan and Euclidean metrics.

use crate::Point;
use std::collections::BTreeMap;


/// A minimum spanning tree of the complete graph on `points` with Manhattan distances: `(total
/// weight, edges)` with `n - 1` index pairs. Candidate edges come from four sweeps over the
/// eight octants (each point is joined to its nearest neighbour in every octant) and Kruskal
/// picks the tree, `O(n log n)`.
pub fn manhattan_mst(points: &[Point]) -> (u64, Vec<(usize, usize)>) {
    let n = points.len();
    let mut ps: Vec<Point> = points.to_vec();
    let mut order: Vec<usize> = (0..n).collect();
    let mut edges: Vec<(i64, usize, usize)> = Vec::new();
    for k in 0..4 {
        order.sort_by_key(|&i| (ps[i].0 + ps[i].1, i));
        let mut sweep: BTreeMap<i64, usize> = BTreeMap::new();
        for &i in &order {
            while let Some((&key, &j)) = sweep.range(-ps[i].1..).next() {
                let (dx, dy) = (ps[i].0 - ps[j].0, ps[i].1 - ps[j].1);
                if dy > dx {
                    break;
                }
                edges.push((dx + dy, i, j));
                sweep.remove(&key);
            }
            sweep.insert(-ps[i].1, i);
        }
        for p in ps.iter_mut() {
            if k & 1 == 1 {
                p.0 = -p.0;
            } else {
                std::mem::swap(&mut p.0, &mut p.1);
            }
        }
    }
    edges.sort_unstable();
    let mut uf = union_find::plain(n);
    let (mut total, mut tree) = (0u64, Vec::with_capacity(n.saturating_sub(1)));
    for (w, u, v) in edges {
        if uf.merge(u, v) {
            total += w as u64;
            tree.push((u, v));
        }
    }
    (total, tree)
}

const LEAF: usize = 8;

struct KdNode {
    lo: usize,
    hi: usize,
    left: usize,
    right: usize,
    bbox: [i64; 4], // min x, max x, min y, max y
}

/// A kd-tree whose points are stored contiguously in tree order (`tx`, `ty`, `order[k]` is the
/// original index of the `k`-th stored point).
struct KdTree {
    order: Vec<usize>,
    tx: Vec<i64>,
    ty: Vec<i64>,
    nodes: Vec<KdNode>,
}

impl KdTree {
    fn new(pts: &[Point]) -> Self {
        let mut t = KdTree { order: (0..pts.len()).collect(), tx: Vec::new(), ty: Vec::new(), nodes: Vec::new() };
        t.build(pts, 0, pts.len());
        t.tx = t.order.iter().map(|&i| pts[i].0).collect();
        t.ty = t.order.iter().map(|&i| pts[i].1).collect();
        t
    }

    fn build(&mut self, pts: &[Point], lo: usize, hi: usize) -> usize {
        let mut bbox = [i64::MAX, i64::MIN, i64::MAX, i64::MIN];
        for &i in &self.order[lo..hi] {
            let (x, y) = pts[i];
            bbox = [bbox[0].min(x), bbox[1].max(x), bbox[2].min(y), bbox[3].max(y)];
        }
        let id = self.nodes.len();
        self.nodes.push(KdNode { lo, hi, left: usize::MAX, right: usize::MAX, bbox });
        if hi - lo > LEAF {
            let by_x = bbox[1] - bbox[0] >= bbox[3] - bbox[2];
            let mid = (lo + hi) / 2;
            self.order[lo..hi].select_nth_unstable_by_key(mid - lo, |&i| if by_x { pts[i].0 } else { pts[i].1 });
            let l = self.build(pts, lo, mid);
            let r = self.build(pts, mid, hi);
            self.nodes[id].left = l;
            self.nodes[id].right = r;
        }
        id
    }

    fn box_dist2(&self, node: usize, p: Point) -> u64 {
        let b = &self.nodes[node].bbox;
        let dx = (b[0] - p.0).max(p.0 - b[1]).max(0) as u64;
        let dy = (b[2] - p.1).max(p.1 - b[3]).max(0) as u64;
        dx * dx + dy * dy
    }

    /// Fills `uniform[node]` with the common component (`tc` is indexed in tree order) of the
    /// node's points or `usize::MAX`.
    fn mark(&self, node: usize, tc: &[u32], uniform: &mut [u32]) -> u32 {
        let n = &self.nodes[node];
        let c = if n.left == usize::MAX {
            let first = tc[n.lo];
            if tc[n.lo..n.hi].iter().all(|&c| c == first) {
                first
            } else {
                u32::MAX
            }
        } else {
            let (a, b) = (self.mark(n.left, tc, uniform), self.mark(n.right, tc, uniform));
            if a == b {
                a
            } else {
                u32::MAX
            }
        };
        uniform[node] = c;
        c
    }

    /// The nearest stored point outside component `c` of the point `qp`, as `(dist2, tree
    /// position)` improving `best` (only strictly closer points replace it).
    fn nearest(&self, node: usize, qp: Point, c: u32, tc: &[u32], uniform: &[u32], best: &mut (u64, usize)) {
        if uniform[node] == c {
            return;
        }
        let n = &self.nodes[node];
        if self.box_dist2(node, qp) >= best.0 {
            return;
        }
        if n.left == usize::MAX {
            for k in n.lo..n.hi {
                if tc[k] != c {
                    let (dx, dy) = (qp.0 - self.tx[k], qp.1 - self.ty[k]);
                    let d = (dx * dx + dy * dy) as u64;
                    if d < best.0 {
                        *best = (d, k);
                    }
                }
            }
            return;
        }
        let (dl, dr) = (self.box_dist2(n.left, qp), self.box_dist2(n.right, qp));
        let (first, second) = if dl <= dr { (n.left, n.right) } else { (n.right, n.left) };
        self.nearest(first, qp, c, tc, uniform, best);
        self.nearest(second, qp, c, tc, uniform, best);
    }
}

/// A minimum spanning tree of the complete Euclidean graph on `points` (coordinates up to about
/// `2^30`): `n - 1` index pairs. Borůvka rounds, each finding for every point the nearest point of
/// another component with a kd-tree whose nodes know when all their points share a component
/// (ties between equal edges are harmless because every added edge is checked against the
/// union-find).
#[allow(clippy::needless_range_loop)]
pub fn euclidean_mst(points: &[Point]) -> Vec<(usize, usize)> {
    let n = points.len();
    let mut tree = Vec::with_capacity(n.saturating_sub(1));
    if n <= 1 {
        return tree;
    }
    let kd = KdTree::new(points);
    let mut uf = union_find::plain(n);
    let mut uniform = vec![u32::MAX; kd.nodes.len()];
    let mut components = n;
    let mut tc = vec![0u32; n];
    while components > 1 {
        for (slot, &i) in tc.iter_mut().zip(&kd.order) {
            *slot = uf.find(i) as u32;
        }
        kd.mark(0, &tc, &mut uniform);
        // the lightest edge leaving each component (ties arbitrary; cycles are skipped below)
        let mut best: Vec<(u64, usize, usize)> = vec![(u64::MAX, usize::MAX, usize::MAX); n];
        for k in 0..n {
            let c = tc[k] as usize;
            let mut found = (best[c].0, usize::MAX);
            kd.nearest(0, (kd.tx[k], kd.ty[k]), c as u32, &tc, &uniform, &mut found);
            if found.1 != usize::MAX {
                best[c] = (found.0, kd.order[k], kd.order[found.1]);
            }
        }
        for c in 0..n {
            if best[c].1 != usize::MAX {
                let (_, u, v) = best[c];
                if uf.merge(u, v) {
                    tree.push((u, v));
                    components -= 1;
                }
            }
        }
    }
    tree
}
