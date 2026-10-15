//! A sequence of keyed items with
//!
//! * [`sort_ascending`](RangeSortSequence::sort_ascending) /
//!   [`sort_descending`](RangeSortSequence::sort_descending) of a range by key,
//! * [`set`](RangeSortSequence::set): replace one item,
//! * [`prod`](RangeSortSequence::prod): fold of the items of a range with a (possibly
//!   non-commutative) [`Monoid`], in sequence order.
//!
//! Every item has a distinct `u64` key below `2^bits` (map your keys to ranks or any monotone
//! `u64` first) and a monoid value. The sequence is kept as maximal *runs* that are each
//! sorted ascending or descending; a run is a dynamic segment tree over the key universe whose
//! nodes store the fold of their keys in ascending **and** in descending key order. Sorting a
//! range cuts the two boundary runs (split by rank) and merges all runs inside by segment-tree
//! merging, which is amortized `O(bits)` per original item. A second segment tree over the
//! *positions* stores each run's fold at the run's start, so a fold cuts its two boundaries and
//! multiplies the runs in between in `O(log n)`, however many there are. Total
//! `O((n + q) * (bits + log n))` amortized.
//!
//! Keys must stay distinct at all times (merging two runs with a common key panics).
//!
//! # Example
//!
//! ```
//! use algebra::Monoid;
//! use range_sort_sequence::RangeSortSequence;
//!
//! // Concatenation of digits: a non-commutative monoid.
//! enum Concat {}
//! impl Monoid for Concat {
//!     type S = Vec<u8>;
//!     fn identity() -> Vec<u8> { vec![] }
//!     fn op(a: &Vec<u8>, b: &Vec<u8>) -> Vec<u8> { [a.as_slice(), b.as_slice()].concat() }
//! }
//!
//! // (key, value) per position.
//! let items = [(30, vec![3]), (10, vec![1]), (20, vec![2]), (40, vec![4])];
//! let mut s = RangeSortSequence::<Concat>::new(6, &items);
//! assert_eq!(s.prod(0, 4), vec![3, 1, 2, 4]);
//! s.sort_ascending(0, 3);
//! assert_eq!(s.prod(0, 4), vec![1, 2, 3, 4]);
//! s.sort_descending(1, 4);
//! assert_eq!(s.prod(0, 4), vec![1, 4, 3, 2]);
//! s.set(0, 5, vec![9]);
//! assert_eq!(s.prod(0, 2), vec![9, 4]);
//! ```

use algebra::Monoid;
use std::collections::BTreeMap;

const NIL: u32 = 0;

struct Node<S> {
    left: u32,
    right: u32,
    cnt: u32,
    /// Fold of the subtree's values in ascending key order.
    asc: S,
    /// Fold of the subtree's values in descending key order.
    desc: S,
}

struct Run {
    end: usize,
    root: u32,
    descending: bool,
}

/// Minimal iterative segment tree for a possibly non-commutative monoid.
struct PositionTree<M: Monoid> {
    size: usize,
    tree: Vec<M::S>,
}

impl<M: Monoid> PositionTree<M> {
    fn new(n: usize) -> Self {
        let size = n.max(1).next_power_of_two();
        PositionTree { size, tree: vec![M::identity(); 2 * size] }
    }

    fn set(&mut self, i: usize, v: M::S) {
        let mut i = i + self.size;
        self.tree[i] = v;
        while i > 1 {
            i >>= 1;
            self.tree[i] = M::op(&self.tree[2 * i], &self.tree[2 * i + 1]);
        }
    }

    fn prod(&self, l: usize, r: usize) -> M::S {
        let (mut l, mut r) = (l + self.size, r + self.size);
        let (mut left, mut right) = (M::identity(), M::identity());
        while l < r {
            if l & 1 == 1 {
                left = M::op(&left, &self.tree[l]);
                l += 1;
            }
            if r & 1 == 1 {
                r -= 1;
                right = M::op(&self.tree[r], &right);
            }
            l >>= 1;
            r >>= 1;
        }
        M::op(&left, &right)
    }
}

/// See the [module documentation](self).
pub struct RangeSortSequence<M: Monoid> {
    n: usize,
    bits: u32,
    /// `nodes[0]` is the empty tree: count 0 and identity folds.
    nodes: Vec<Node<M::S>>,
    free: Vec<u32>,
    /// Run start -> run; the runs tile `0..n`.
    runs: BTreeMap<usize, Run>,
    /// The fold of the run starting at each position (the identity at non-starts).
    by_position: PositionTree<M>,
}

impl<M: Monoid> RangeSortSequence<M> {
    /// A sequence of the given `(key, value)` items; keys are distinct and `< 2^bits`.
    ///
    /// # Panics
    /// Panics if `bits` is not in `1..=63`, a key is `>= 2^bits`, or keys repeat.
    pub fn new(bits: u32, items: &[(u64, M::S)]) -> Self {
        assert!((1..=63).contains(&bits), "bits must be in 1..=63");
        assert!(items.len() < u32::MAX as usize / (bits as usize + 1), "too many items");
        let nil = Node { left: NIL, right: NIL, cnt: 0, asc: M::identity(), desc: M::identity() };
        let mut s = RangeSortSequence { n: items.len(), bits, nodes: vec![nil], free: Vec::new(), runs: BTreeMap::new(), by_position: PositionTree::new(items.len()) };
        let mut keys: Vec<u64> = items.iter().map(|it| it.0).collect();
        keys.sort_unstable();
        assert!(keys.windows(2).all(|w| w[0] != w[1]), "keys must be distinct");
        for (i, (key, value)) in items.iter().enumerate() {
            let root = s.new_path(*key, value.clone());
            s.runs.insert(i, Run { end: i + 1, root, descending: false });
            s.by_position.tree[s.by_position.size + i] = value.clone();
        }
        for i in (1..s.by_position.size).rev() {
            s.by_position.tree[i] = M::op(&s.by_position.tree[2 * i], &s.by_position.tree[2 * i + 1]);
        }
        s
    }

    /// The fold of the whole run rooted at `root`, in the run's own order.
    fn run_fold(&self, root: u32, descending: bool) -> M::S {
        let node = &self.nodes[root as usize];
        if descending { node.desc.clone() } else { node.asc.clone() }
    }

    /// Record in `by_position` that a run starting at `start` now has the given tree/direction.
    fn index_run(&mut self, start: usize, root: u32, descending: bool) {
        let fold = self.run_fold(root, descending);
        self.by_position.set(start, fold);
    }

    /// Number of items.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    fn alloc(&mut self, node: Node<M::S>) -> u32 {
        if let Some(id) = self.free.pop() {
            self.nodes[id as usize] = node;
            id
        } else {
            self.nodes.push(node);
            (self.nodes.len() - 1) as u32
        }
    }

    /// A chain of nodes from the root down to the leaf of `key`.
    fn new_path(&mut self, key: u64, value: M::S) -> u32 {
        assert!(key >> self.bits == 0, "key {key} does not fit in {} bits", self.bits);
        let mut u = self.alloc(Node { left: NIL, right: NIL, cnt: 1, asc: value.clone(), desc: value });
        for level in 1..=self.bits {
            let (left, right) = if (key >> (level - 1)) & 1 == 0 { (u, NIL) } else { (NIL, u) };
            let parent = self.alloc(Node { left, right, cnt: 0, asc: M::identity(), desc: M::identity() });
            self.pull(parent);
            u = parent;
        }
        u
    }

    fn pull(&mut self, u: u32) {
        let (l, r) = (self.nodes[u as usize].left as usize, self.nodes[u as usize].right as usize);
        let cnt = self.nodes[l].cnt + self.nodes[r].cnt;
        let asc = M::op(&self.nodes[l].asc, &self.nodes[r].asc);
        let desc = M::op(&self.nodes[r].desc, &self.nodes[l].desc);
        let node = &mut self.nodes[u as usize];
        node.cnt = cnt;
        node.asc = asc;
        node.desc = desc;
    }

    /// Split into the `k` smallest keys and the rest.
    fn split(&mut self, u: u32, k: u32) -> (u32, u32) {
        let cnt = self.nodes[u as usize].cnt;
        if k == 0 {
            return (NIL, u);
        }
        if k >= cnt {
            return (u, NIL);
        }
        // 0 < k < cnt: `u` is an internal node.
        let (l, r) = (self.nodes[u as usize].left, self.nodes[u as usize].right);
        let cl = self.nodes[l as usize].cnt;
        let blank = || Node { left: NIL, right: NIL, cnt: 0, asc: M::identity(), desc: M::identity() };
        if k <= cl {
            let (a, b) = self.split(l, k);
            let new = self.alloc(blank());
            self.nodes[new as usize].left = a;
            self.pull(new);
            self.nodes[u as usize].left = b;
            self.pull(u);
            (new, u)
        } else {
            let (a, b) = self.split(r, k - cl);
            let new = self.alloc(blank());
            self.nodes[new as usize].right = b;
            self.pull(new);
            self.nodes[u as usize].right = a;
            self.pull(u);
            (u, new)
        }
    }

    /// Union of two key-disjoint trees of height `level`.
    fn merge(&mut self, a: u32, b: u32, level: u32) -> u32 {
        if a == NIL {
            return b;
        }
        if b == NIL {
            return a;
        }
        assert!(level > 0, "duplicate key");
        let (al, ar) = (self.nodes[a as usize].left, self.nodes[a as usize].right);
        let (bl, br) = (self.nodes[b as usize].left, self.nodes[b as usize].right);
        let l = self.merge(al, bl, level - 1);
        let r = self.merge(ar, br, level - 1);
        self.nodes[a as usize].left = l;
        self.nodes[a as usize].right = r;
        self.free.push(b);
        self.pull(a);
        a
    }

    /// Make `pos` a run boundary (`0 < pos < n`).
    fn cut(&mut self, pos: usize) {
        if pos == 0 || pos >= self.n {
            return;
        }
        let (&start, run) = self.runs.range(..=pos).next_back().expect("runs tile the sequence");
        if start == pos {
            return;
        }
        let (end, root, descending) = (run.end, run.root, run.descending);
        let total = self.nodes[root as usize].cnt;
        let first_len = (pos - start) as u32;
        let (first, second) = if descending {
            // The first positions hold the largest keys.
            let (low, high) = self.split(root, total - first_len);
            (high, low)
        } else {
            self.split(root, first_len)
        };
        self.runs.insert(start, Run { end: pos, root: first, descending });
        self.runs.insert(pos, Run { end, root: second, descending });
        self.index_run(start, first, descending);
        self.index_run(pos, second, descending);
    }

    fn check_range(&self, l: usize, r: usize) {
        assert!(l <= r && r <= self.n, "invalid range [{l}, {r}) for length {}", self.n);
    }

    fn sort(&mut self, l: usize, r: usize, descending: bool) {
        self.check_range(l, r);
        if r - l < 2 {
            return;
        }
        self.cut(l);
        self.cut(r);
        let starts: Vec<usize> = self.runs.range(l..r).map(|(&s, _)| s).collect();
        let mut root = NIL;
        for s in starts {
            let run = self.runs.remove(&s).expect("run exists");
            root = self.merge(root, run.root, self.bits);
            if s != l {
                self.by_position.set(s, M::identity());
            }
        }
        self.runs.insert(l, Run { end: r, root, descending });
        self.index_run(l, root, descending);
    }

    /// Sort `[l, r)` by ascending key.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn sort_ascending(&mut self, l: usize, r: usize) {
        self.sort(l, r, false);
    }

    /// Sort `[l, r)` by descending key.
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn sort_descending(&mut self, l: usize, r: usize) {
        self.sort(l, r, true);
    }

    /// Replace the item at position `i` by `(key, value)`. The key must differ from all other
    /// current keys.
    ///
    /// # Panics
    /// Panics if `i >= len()` or the key does not fit in `bits` bits.
    pub fn set(&mut self, i: usize, key: u64, value: M::S) {
        assert!(i < self.n, "index {i} out of range for length {}", self.n);
        assert!(key >> self.bits == 0, "key {key} does not fit in {} bits", self.bits);
        self.cut(i);
        self.cut(i + 1);
        let old = self.runs.remove(&i).expect("a run starts at i");
        self.release(old.root);
        let root = self.new_path(key, value.clone());
        self.runs.insert(i, Run { end: i + 1, root, descending: false });
        self.by_position.set(i, value);
    }

    /// Return the nodes of a (single-item) tree to the free list.
    fn release(&mut self, mut u: u32) {
        while u != NIL {
            self.free.push(u);
            let node = &self.nodes[u as usize];
            u = if node.left != NIL { node.left } else { node.right };
        }
    }

    /// The fold of the values at positions `l..r`, in sequence order (the identity for an
    /// empty range).
    ///
    /// # Panics
    /// Panics unless `l <= r <= len()`.
    pub fn prod(&mut self, l: usize, r: usize) -> M::S {
        self.check_range(l, r);
        if l == r {
            return M::identity();
        }
        self.cut(l);
        self.cut(r);
        self.by_position.prod(l, r)
    }
}
