//! Range frequency queries with point assignments: how many times does `x` occur in
//! `a[l..r]`, while `a[p] = v` assignments keep arriving?
//!
//! The structure is *offline in the assignments only*: the constructor needs the list of all
//! `(position, value)` pairs that will ever be assigned (the queries themselves stay online).
//! It then keeps, for every value, a Fenwick tree over the sorted positions where that value
//! can ever appear, all stored in two flat arrays. An assignment is two Fenwick updates
//! and a count is two Fenwick prefix sums plus three binary searches, `O(log n)` each, with
//! `O(n + q)` memory in total.
//!
//! # Example
//!
//! ```
//! use point_set_range_frequency::PointSetRangeFrequency;
//!
//! // a = [5, 7, 5, 5], and we will later assign a[1] = 5 and a[3] = 7.
//! let mut f = PointSetRangeFrequency::new(&[5, 7, 5, 5], &[(1, 5), (3, 7)]);
//! assert_eq!(f.count(0, 4, 5), 3);
//! f.set(1, 5);
//! assert_eq!(f.count(0, 4, 5), 4);
//! f.set(3, 7);
//! assert_eq!(f.count(2, 4, 7), 1);
//! assert_eq!(f.count(0, 4, 9), 0); // values never present are fine
//! ```

/// See the [module documentation](self).
pub struct PointSetRangeFrequency<T> {
    a: Vec<T>,
    /// Sorted distinct values that are ever present.
    values: Vec<T>,
    /// Group `g` owns `pos[start[g]..start[g + 1]]` and the same range of `tree`.
    start: Vec<usize>,
    /// Sorted positions where the group's value can appear.
    pos: Vec<u32>,
    /// Fenwick trees (1-based inside each group) over the occurrence indicators.
    tree: Vec<i32>,
}

impl<T: Ord + Copy> PointSetRangeFrequency<T> {
    /// Start from `initial`; `updates` lists every `(position, value)` that [`set`](Self::set)
    /// may be called with later.
    ///
    /// # Panics
    /// Panics if an update position is out of range or `initial.len() >= 2^32`.
    pub fn new(initial: &[T], updates: &[(usize, T)]) -> Self {
        let n = initial.len();
        assert!(n < u32::MAX as usize, "array too long");
        let mut pairs: Vec<(T, u32)> = Vec::with_capacity(n + updates.len());
        pairs.extend(initial.iter().enumerate().map(|(i, &v)| (v, i as u32)));
        for &(p, v) in updates {
            assert!(p < n, "update position {p} out of range for length {n}");
            pairs.push((v, p as u32));
        }
        pairs.sort_unstable();
        pairs.dedup();
        let mut values = Vec::new();
        let mut start = Vec::new();
        for (i, &(v, _)) in pairs.iter().enumerate() {
            if i == 0 || pairs[i - 1].0 != v {
                values.push(v);
                start.push(i);
            }
        }
        start.push(pairs.len());
        let pos: Vec<u32> = pairs.iter().map(|&(_, p)| p).collect();
        let mut s = PointSetRangeFrequency { a: initial.to_vec(), values, start, pos, tree: vec![0; pairs.len()] };
        // Linear Fenwick build per group.
        for (i, &v) in initial.iter().enumerate() {
            let g = s.group(&v).expect("initial values have a group");
            let j = s.index_of(g, i);
            s.tree[s.start[g] + j] = 1;
        }
        for g in 0..s.values.len() {
            let (lo, hi) = (s.start[g], s.start[g + 1]);
            for i in 1..=hi - lo {
                let parent = i + (i & i.wrapping_neg());
                if parent <= hi - lo {
                    s.tree[lo + parent - 1] += s.tree[lo + i - 1];
                }
            }
        }
        s
    }

    pub fn len(&self) -> usize {
        self.a.len()
    }

    pub fn is_empty(&self) -> bool {
        self.a.is_empty()
    }

    /// Current value at position `i`.
    pub fn get(&self, i: usize) -> T {
        assert!(i < self.a.len(), "index {i} out of range for length {}", self.a.len());
        self.a[i]
    }

    fn group(&self, v: &T) -> Option<usize> {
        self.values.binary_search(v).ok()
    }

    /// Index of position `p` inside group `g`'s sorted positions (which must contain it).
    fn index_of(&self, g: usize, p: usize) -> usize {
        let slice = &self.pos[self.start[g]..self.start[g + 1]];
        slice.binary_search(&(p as u32)).expect("(position, value) pair was not declared in `new`")
    }

    fn add(&mut self, g: usize, idx: usize, delta: i32) {
        let (lo, m) = (self.start[g], self.start[g + 1] - self.start[g]);
        let mut i = idx + 1;
        while i <= m {
            self.tree[lo + i - 1] += delta;
            i += i & i.wrapping_neg();
        }
    }

    /// Sum of the first `k` indicators of group `g`.
    fn prefix(&self, g: usize, mut k: usize) -> i32 {
        let lo = self.start[g];
        let mut s = 0;
        while k > 0 {
            s += self.tree[lo + k - 1];
            k &= k - 1;
        }
        s
    }

    /// Assign `a[p] = v`.
    ///
    /// # Panics
    /// Panics if `p` is out of range or `(p, v)` was not declared in [`new`](Self::new).
    pub fn set(&mut self, p: usize, v: T) {
        assert!(p < self.a.len(), "index {p} out of range for length {}", self.a.len());
        let old = self.a[p];
        if old == v {
            return;
        }
        let (go, gn) = (self.group(&old).unwrap(), self.group(&v).expect("(position, value) pair was not declared in `new`"));
        let (io, i_n) = (self.index_of(go, p), self.index_of(gn, p));
        self.add(go, io, -1);
        self.add(gn, i_n, 1);
        self.a[p] = v;
    }

    /// The number of `i` in `l..r` with `a[i] == x`.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn count(&self, l: usize, r: usize, x: T) -> usize {
        assert!(l <= r && r <= self.a.len(), "invalid range [{l}, {r}) for length {}", self.a.len());
        let Some(g) = self.group(&x) else { return 0 };
        let slice = &self.pos[self.start[g]..self.start[g + 1]];
        let (lo, hi) = (slice.partition_point(|&p| (p as usize) < l), slice.partition_point(|&p| (p as usize) < r));
        (self.prefix(g, hi) - self.prefix(g, lo)) as usize
    }
}
