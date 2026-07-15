//! A treap-based sorted list.
//!
//! Supports the following operations, each `O(log n)`:
//! - [`SortedList::insert`] / [`SortedList::remove`]
//! - [`SortedList::bisect_left`] / [`SortedList::bisect_right`]
//! - [`SortedList::get`]
//!
//! [`SortedList::len`] and [`SortedList::is_empty`] are `O(1)`.
//!
//! Implementation note: nodes are stored in a flat `Vec<Node>` addressed by
//! 1-based indices (0 = none). Removed nodes go to a `free_stack` so the
//! underlying storage can be reused without growing the array.

use std::cmp::Ordering;

#[derive(Debug)]
struct Node {
    key: i64,
    priority: u32,
    left: usize,
    right: usize,
    size: usize,
}

impl Node {
    fn new(key: i64, seed: &mut u32) -> Self {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 17;
        *seed ^= *seed << 5;
        Self {
            key,
            priority: *seed,
            left: 0,
            right: 0,
            size: 1,
        }
    }

    fn update(idx: usize, nodes: &mut [Node]) {
        if idx == 0 {
            return;
        }
        let l = nodes[idx - 1].left;
        let r = nodes[idx - 1].right;
        let l_size = if l == 0 { 0 } else { nodes[l - 1].size };
        let r_size = if r == 0 { 0 } else { nodes[r - 1].size };
        nodes[idx - 1].size = 1 + l_size + r_size;
    }
}

/// Treap-backed sorted list of `i64` keys with multiplicity.
pub struct SortedList {
    nodes: Vec<Node>,
    root: usize,
    seed: u32,
    free_stack: Vec<usize>,
}

impl SortedList {
    /// Construct an empty sorted list.
    pub fn new() -> Self {
        Self {
            nodes: Vec::with_capacity(1024 * 8),
            root: 0,
            seed: 12345,
            free_stack: Vec::new(),
        }
    }

    /// Split by key. If `inclusive` is true, keys == key go to the left tree.
    /// Returns `(left_idx, right_idx)` where indices are 1-based (0 means None).
    fn split(&mut self, node_idx: usize, key: i64, inclusive: bool) -> (usize, usize) {
        if node_idx == 0 {
            return (0, 0);
        }

        let n_key = self.nodes[node_idx - 1].key;
        let go_left = if inclusive { n_key >= key } else { n_key > key };

        if go_left {
            let left_child = self.nodes[node_idx - 1].left;
            let (l, r) = self.split(left_child, key, inclusive);
            self.nodes[node_idx - 1].left = r;
            Node::update(node_idx, &mut self.nodes);
            (l, node_idx)
        } else {
            let right_child = self.nodes[node_idx - 1].right;
            let (l, r) = self.split(right_child, key, inclusive);
            self.nodes[node_idx - 1].right = l;
            Node::update(node_idx, &mut self.nodes);
            (node_idx, r)
        }
    }

    /// Merge two trees where all keys in `l` are less than all keys in `r`.
    fn merge(&mut self, l: usize, r: usize) -> usize {
        if l == 0 || r == 0 {
            return if l == 0 { r } else { l };
        }

        if self.nodes[l - 1].priority > self.nodes[r - 1].priority {
            let l_right = self.nodes[l - 1].right;
            self.nodes[l - 1].right = self.merge(l_right, r);
            Node::update(l, &mut self.nodes);
            l
        } else {
            let r_left = self.nodes[r - 1].left;
            self.nodes[r - 1].left = self.merge(l, r_left);
            Node::update(r, &mut self.nodes);
            r
        }
    }

    /// Insert a key. Duplicates are allowed.
    pub fn insert(&mut self, key: i64) {
        let (l, r) = self.split(self.root, key, true);
        let node_idx = if let Some(idx) = self.free_stack.pop() {
            self.nodes[idx - 1] = Node::new(key, &mut self.seed);
            idx
        } else {
            self.nodes.push(Node::new(key, &mut self.seed));
            self.nodes.len()
        };
        let merged_left = self.merge(l, node_idx);
        self.root = self.merge(merged_left, r);
    }

    /// Remove one occurrence of `key`.
    pub fn remove(&mut self, key: i64) {
        let (l, mid_r) = self.split(self.root, key, true);
        let (mid, r) = self.split(mid_r, key + 1, true);
        if mid != 0 {
            let m_left = self.nodes[mid - 1].left;
            let m_right = self.nodes[mid - 1].right;
            let new_mid = self.merge(m_left, m_right);
            self.free_stack.push(mid);
            let merged_left = self.merge(l, new_mid);
            self.root = self.merge(merged_left, r);
        } else {
            self.root = self.merge(l, r);
        }
    }

    /// Number of keys strictly less than `key`.
    pub fn bisect_left(&self, key: i64) -> usize {
        let mut curr = self.root;
        let mut rank = 0;
        while curr != 0 {
            let node = &self.nodes[curr - 1];
            if node.key >= key {
                curr = node.left;
            } else {
                let l_size = if node.left == 0 { 0 } else { self.nodes[node.left - 1].size };
                rank += 1 + l_size;
                curr = node.right;
            }
        }
        rank
    }

    /// Number of keys less than or equal to `key`.
    pub fn bisect_right(&self, key: i64) -> usize {
        let mut curr = self.root;
        let mut rank = 0;
        while curr != 0 {
            let node = &self.nodes[curr - 1];
            if node.key > key {
                curr = node.left;
            } else {
                let l_size = if node.left == 0 { 0 } else { self.nodes[node.left - 1].size };
                rank += 1 + l_size;
                curr = node.right;
            }
        }
        rank
    }

    /// Get the element at `index` (0-based), or `None` if out of bounds.
    pub fn get(&self, index: usize) -> Option<i64> {
        let mut curr = self.root;
        let mut i = index;
        while curr != 0 {
            let node = &self.nodes[curr - 1];
            let left_size = if node.left == 0 { 0 } else { self.nodes[node.left - 1].size };
            match i.cmp(&left_size) {
                Ordering::Less => curr = node.left,
                Ordering::Equal => return Some(node.key),
                Ordering::Greater => {
                    i -= left_size + 1;
                    curr = node.right;
                }
            }
        }
        None
    }

    /// Number of elements in the list (with multiplicities).
    pub fn len(&self) -> usize {
        if self.root == 0 {
            0
        } else {
            self.nodes[self.root - 1].size
        }
    }

    /// Whether the list is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}