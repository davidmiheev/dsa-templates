//! A generic segment tree.
//!
//! `T` is the type of the data in the array; `F` is the merge function
//! (typically a sum, min, max, or xor closure).
//!
//! Both `query` and `update` run in `O(log n)`. The tree is stored in
//! a flat `Vec<T>` of size `4 * n`, which is the standard bound for an
//! implicit binary-tree representation.

/// Generic segment tree.
pub struct SegmentTree<T, F> {
    size: usize,
    tree: Vec<T>,
    default: T,
    merge: F,
}

impl<T, F> SegmentTree<T, F>
where
    T: Clone + Copy,
    F: Fn(T, T) -> T,
{
    /// Build a segment tree over `arr` using `default` as the identity
    /// for `merge` and `merge` itself as the associative combining op.
    pub fn new(arr: &[T], default: T, merge: F) -> Self {
        let size = arr.len();
        #[allow(unused_mut)]
        let mut tree = vec![default; 4 * size];
        let mut seg_tree = SegmentTree {
            size,
            tree,
            default,
            merge,
        };
        seg_tree.build(arr, 0, 0, size - 1);
        seg_tree
    }

    fn build(&mut self, arr: &[T], node: usize, start: usize, end: usize) {
        if start == end {
            self.tree[node] = arr[start];
            return;
        }
        let mid = (start + end) / 2;
        self.build(arr, 2 * node + 1, start, mid);
        self.build(arr, 2 * node + 2, mid + 1, end);
        self.tree[node] = (self.merge)(self.tree[2 * node + 1], self.tree[2 * node + 2]);
    }

    /// Query `merge` over the closed range `[l, r]`.
    pub fn query(&self, l: usize, r: usize) -> T {
        self.query_recursive(0, 0, self.size - 1, l, r)
    }

    fn query_recursive(&self, node: usize, start: usize, end: usize, l: usize, r: usize) -> T {
        if r < start || end < l {
            return self.default;
        }
        if l <= start && end <= r {
            return self.tree[node];
        }
        let mid = (start + end) / 2;
        let p1 = self.query_recursive(2 * node + 1, start, mid, l, r);
        let p2 = self.query_recursive(2 * node + 2, mid + 1, end, l, r);
        (self.merge)(p1, p2)
    }

    /// Point update: set `arr[idx] = val` and recompute affected nodes.
    pub fn update(&mut self, idx: usize, val: T) {
        self.update_recursive(0, 0, self.size - 1, idx, val);
    }

    fn update_recursive(&mut self, node: usize, start: usize, end: usize, idx: usize, val: T) {
        if start == end {
            self.tree[node] = val;
            return;
        }
        let mid = (start + end) / 2;
        if start <= idx && idx <= mid {
            self.update_recursive(2 * node + 1, start, mid, idx, val);
        } else {
            self.update_recursive(2 * node + 2, mid + 1, end, idx, val);
        }
        self.tree[node] = (self.merge)(self.tree[2 * node + 1], self.tree[2 * node + 2]);
    }
}