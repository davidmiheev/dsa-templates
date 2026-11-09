//! Segment tree beats (Ji driver segment tree) over `i64`.
//!
//! Supports on a range `[l, r)`:
//!
//! * [`chmin`](SegmentTreeBeats::chmin) `a_i = min(a_i, x)`,
//! * [`chmax`](SegmentTreeBeats::chmax) `a_i = max(a_i, x)`,
//! * [`add`](SegmentTreeBeats::add) `a_i += x`,
//!
//! and the queries [`sum`](SegmentTreeBeats::sum), [`min`](SegmentTreeBeats::min),
//! [`max`](SegmentTreeBeats::max). Every node keeps the two largest and two smallest distinct
//! values with their multiplicities, so a `chmin` only has to "dig deeper" where it changes
//! more than the maximum; the total cost is `O((n + q) log^2 n)` amortized (`O(log n)` without
//! `add`).
//!
//! Values are `i64`; the structure itself uses `i64::MIN` / `i64::MAX` as "absent second
//! value" sentinels, so stored values must stay strictly between them.
//!
//! # Example
//!
//! ```
//! use segment_tree_beats::SegmentTreeBeats;
//!
//! let mut t = SegmentTreeBeats::new(&[5, 1, 9, 3]);
//! t.chmin(0, 4, 4); // [4, 1, 4, 3]
//! t.chmax(1, 3, 2); // [4, 2, 4, 3]
//! t.add(2, 4, 10);  // [4, 2, 14, 13]
//! assert_eq!(t.sum(0, 4), 33);
//! assert_eq!((t.min(0, 4), t.max(1, 3)), (2, 14));
//! ```

const NEG: i64 = i64::MIN;
const POS: i64 = i64::MAX;

#[derive(Clone, Copy)]
struct Node {
    sum: i64,
    max1: i64,
    max2: i64,
    max_cnt: u32,
    min1: i64,
    min2: i64,
    min_cnt: u32,
    /// Pending addition for the children.
    add: i64,
}

impl Node {
    fn leaf(x: i64) -> Node {
        Node { sum: x, max1: x, max2: NEG, max_cnt: 1, min1: x, min2: POS, min_cnt: 1, add: 0 }
    }

    fn apply_add(&mut self, len: i64, x: i64) {
        self.sum += x * len;
        self.max1 += x;
        if self.max2 != NEG {
            self.max2 += x;
        }
        self.min1 += x;
        if self.min2 != POS {
            self.min2 += x;
        }
        self.add += x;
    }

    /// Requires `self.max2 < x < self.max1`.
    fn apply_chmin(&mut self, x: i64) {
        self.sum -= (self.max1 - x) * self.max_cnt as i64;
        if self.min1 == self.max1 {
            self.min1 = x;
        }
        if self.min2 == self.max1 {
            self.min2 = x;
        }
        self.max1 = x;
    }

    /// Requires `self.min1 < x < self.min2`.
    fn apply_chmax(&mut self, x: i64) {
        self.sum += (x - self.min1) * self.min_cnt as i64;
        if self.max1 == self.min1 {
            self.max1 = x;
        }
        if self.max2 == self.min1 {
            self.max2 = x;
        }
        self.min1 = x;
    }

    fn merge(l: &Node, r: &Node) -> Node {
        let (max1, max2, max_cnt) = if l.max1 > r.max1 {
            (l.max1, l.max2.max(r.max1), l.max_cnt)
        } else if l.max1 < r.max1 {
            (r.max1, r.max2.max(l.max1), r.max_cnt)
        } else {
            (l.max1, l.max2.max(r.max2), l.max_cnt + r.max_cnt)
        };
        let (min1, min2, min_cnt) = if l.min1 < r.min1 {
            (l.min1, l.min2.min(r.min1), l.min_cnt)
        } else if l.min1 > r.min1 {
            (r.min1, r.min2.min(l.min1), r.min_cnt)
        } else {
            (l.min1, l.min2.min(r.min2), l.min_cnt + r.min_cnt)
        };
        Node { sum: l.sum + r.sum, max1, max2, max_cnt, min1, min2, min_cnt, add: 0 }
    }
}

/// See the [module documentation](self).
pub struct SegmentTreeBeats {
    n: usize,
    nodes: Vec<Node>,
}

impl SegmentTreeBeats {
    /// Build over `a` in `O(n)`.
    ///
    /// # Panics
    /// Panics if `a` is empty.
    pub fn new(a: &[i64]) -> Self {
        assert!(!a.is_empty(), "the array must not be empty");
        let mut t = SegmentTreeBeats { n: a.len(), nodes: vec![Node::leaf(0); 4 * a.len()] };
        t.build(1, 0, a.len(), a);
        t
    }

    fn build(&mut self, u: usize, lo: usize, hi: usize, a: &[i64]) {
        if hi - lo == 1 {
            self.nodes[u] = Node::leaf(a[lo]);
            return;
        }
        let mid = (lo + hi) / 2;
        self.build(2 * u, lo, mid, a);
        self.build(2 * u + 1, mid, hi, a);
        self.nodes[u] = Node::merge(&self.nodes[2 * u], &self.nodes[2 * u + 1]);
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    fn check(&self, l: usize, r: usize) {
        assert!(l <= r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
    }

    fn push_down(&mut self, u: usize, lo: usize, hi: usize) {
        let mid = (lo + hi) / 2;
        let node = self.nodes[u];
        for (c, len) in [(2 * u, mid - lo), (2 * u + 1, hi - mid)] {
            if node.add != 0 {
                self.nodes[c].apply_add(len as i64, node.add);
            }
            if self.nodes[c].max1 > node.max1 {
                self.nodes[c].apply_chmin(node.max1);
            }
            if self.nodes[c].min1 < node.min1 {
                self.nodes[c].apply_chmax(node.min1);
            }
        }
        self.nodes[u].add = 0;
    }

    fn pull(&mut self, u: usize) {
        let add = self.nodes[u].add;
        self.nodes[u] = Node::merge(&self.nodes[2 * u], &self.nodes[2 * u + 1]);
        self.nodes[u].add = add;
    }

    /// `a_i = min(a_i, x)` for `i` in `[l, r)`.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn chmin(&mut self, l: usize, r: usize, x: i64) {
        self.check(l, r);
        if l < r {
            self.chmin_rec(1, 0, self.n, l, r, x);
        }
    }

    fn chmin_rec(&mut self, u: usize, lo: usize, hi: usize, l: usize, r: usize, x: i64) {
        if r <= lo || hi <= l || self.nodes[u].max1 <= x {
            return;
        }
        if l <= lo && hi <= r && self.nodes[u].max2 < x {
            self.nodes[u].apply_chmin(x);
            return;
        }
        self.push_down(u, lo, hi);
        let mid = (lo + hi) / 2;
        self.chmin_rec(2 * u, lo, mid, l, r, x);
        self.chmin_rec(2 * u + 1, mid, hi, l, r, x);
        self.pull(u);
    }

    /// `a_i = max(a_i, x)` for `i` in `[l, r)`.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn chmax(&mut self, l: usize, r: usize, x: i64) {
        self.check(l, r);
        if l < r {
            self.chmax_rec(1, 0, self.n, l, r, x);
        }
    }

    fn chmax_rec(&mut self, u: usize, lo: usize, hi: usize, l: usize, r: usize, x: i64) {
        if r <= lo || hi <= l || self.nodes[u].min1 >= x {
            return;
        }
        if l <= lo && hi <= r && self.nodes[u].min2 > x {
            self.nodes[u].apply_chmax(x);
            return;
        }
        self.push_down(u, lo, hi);
        let mid = (lo + hi) / 2;
        self.chmax_rec(2 * u, lo, mid, l, r, x);
        self.chmax_rec(2 * u + 1, mid, hi, l, r, x);
        self.pull(u);
    }

    /// `a_i += x` for `i` in `[l, r)`.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn add(&mut self, l: usize, r: usize, x: i64) {
        self.check(l, r);
        if l < r {
            self.add_rec(1, 0, self.n, l, r, x);
        }
    }

    fn add_rec(&mut self, u: usize, lo: usize, hi: usize, l: usize, r: usize, x: i64) {
        if r <= lo || hi <= l {
            return;
        }
        if l <= lo && hi <= r {
            self.nodes[u].apply_add((hi - lo) as i64, x);
            return;
        }
        self.push_down(u, lo, hi);
        let mid = (lo + hi) / 2;
        self.add_rec(2 * u, lo, mid, l, r, x);
        self.add_rec(2 * u + 1, mid, hi, l, r, x);
        self.pull(u);
    }

    /// Sum of `a[l..r]` (0 for an empty range).
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn sum(&mut self, l: usize, r: usize) -> i64 {
        self.check(l, r);
        if l == r {
            return 0;
        }
        self.sum_rec(1, 0, self.n, l, r)
    }

    fn sum_rec(&mut self, u: usize, lo: usize, hi: usize, l: usize, r: usize) -> i64 {
        if r <= lo || hi <= l {
            return 0;
        }
        if l <= lo && hi <= r {
            return self.nodes[u].sum;
        }
        self.push_down(u, lo, hi);
        let mid = (lo + hi) / 2;
        self.sum_rec(2 * u, lo, mid, l, r) + self.sum_rec(2 * u + 1, mid, hi, l, r)
    }

    /// Minimum of `a[l..r]`.
    ///
    /// # Panics
    /// Panics unless `l < r <= len()`.
    pub fn min(&mut self, l: usize, r: usize) -> i64 {
        assert!(l < r && r <= self.n, "invalid non-empty range [{l}, {r}) for length {}", self.n);
        self.min_rec(1, 0, self.n, l, r)
    }

    fn min_rec(&mut self, u: usize, lo: usize, hi: usize, l: usize, r: usize) -> i64 {
        if r <= lo || hi <= l {
            return POS;
        }
        if l <= lo && hi <= r {
            return self.nodes[u].min1;
        }
        self.push_down(u, lo, hi);
        let mid = (lo + hi) / 2;
        self.min_rec(2 * u, lo, mid, l, r).min(self.min_rec(2 * u + 1, mid, hi, l, r))
    }

    /// Maximum of `a[l..r]`.
    ///
    /// # Panics
    /// Panics unless `l < r <= len()`.
    pub fn max(&mut self, l: usize, r: usize) -> i64 {
        assert!(l < r && r <= self.n, "invalid non-empty range [{l}, {r}) for length {}", self.n);
        self.max_rec(1, 0, self.n, l, r)
    }

    fn max_rec(&mut self, u: usize, lo: usize, hi: usize, l: usize, r: usize) -> i64 {
        if r <= lo || hi <= l {
            return NEG;
        }
        if l <= lo && hi <= r {
            return self.nodes[u].max1;
        }
        self.push_down(u, lo, hi);
        let mid = (lo + hi) / 2;
        self.max_rec(2 * u, lo, mid, l, r).max(self.max_rec(2 * u + 1, mid, hi, l, r))
    }

    /// Element `i`.
    pub fn get(&mut self, i: usize) -> i64 {
        assert!(i < self.n, "index {i} out of range for length {}", self.n);
        self.sum(i, i + 1)
    }
}
