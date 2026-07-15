//! Generic binary tree with iterative traversals.
//!
//! This crate provides a small `BinaryTree<T>` type backed by `Vec<Option<T>>`
//! for use as a quick reference / heap-storage representation of a complete
//! binary tree (the same convention used by competitive-programming templates
//! where the root lives at index 1).
//!
//! Provided operations:
//!
//! * [`BinaryTree::new`] — allocate an empty tree with a fixed capacity
//! * [`BinaryTree::set`] — write a value at an index
//! * [`BinaryTree::get`] — read a value at an index
//! * [`BinaryTree::inorder`], [`BinaryTree::preorder`], [`BinaryTree::postorder`] —
//!   iterative traversals over a 1-indexed tree stored at indices `1..=n`
//! * [`BinaryTree::level_order`] — breadth-first traversal
//!
//! Indices follow the standard heap convention: for a node at index `i`,
//! the left child is `2i` and the right child is `2i + 1`.

/// A binary tree stored in an indexed `Vec`, following the heap convention
/// where node `i` has children at `2i` and `2i + 1`.
#[derive(Debug, Clone)]
pub struct BinaryTree<T> {
    data: Vec<Option<T>>,
}

impl<T: Clone> BinaryTree<T> {
    /// Create a tree that can hold nodes at indices `1..=size`.
    pub fn new(size: usize) -> Self {
        let mut data = Vec::with_capacity(size + 1);
        data.push(None);
        for _ in 0..size {
            data.push(None);
        }
        Self { data }
    }

    /// Write `value` at index `idx`.
    ///
    /// # Panics
    ///
    /// Panics if `idx == 0` or `idx >= data.len()`.
    pub fn set(&mut self, idx: usize, value: T) {
        assert!(idx > 0 && idx < self.data.len(), "index out of range");
        self.data[idx] = Some(value);
    }

    /// Read the value at index `idx`.
    ///
    /// # Panics
    ///
    /// Panics if `idx == 0` or `idx >= data.len()`.
    pub fn get(&self, idx: usize) -> Option<&T> {
        assert!(idx > 0 && idx < self.data.len(), "index out of range");
        self.data[idx].as_ref()
    }

    /// In-order traversal (left, root, right) of the nodes at `1..=n`.
    ///
    /// Only nodes that have been set (via [`set`](Self::set)) are yielded.
    pub fn inorder(&self, n: usize) -> Vec<T> {
        let mut out = Vec::new();
        self.inorder_rec(1, n, &mut out);
        out
    }

    fn inorder_rec(&self, idx: usize, n: usize, out: &mut Vec<T>) {
        if idx > n {
            return;
        }
        self.inorder_rec(idx * 2, n, out);
        if let Some(v) = &self.data[idx] {
            out.push(v.clone());
        }
        self.inorder_rec(idx * 2 + 1, n, out);
    }

    /// Pre-order traversal (root, left, right) of the nodes at `1..=n`.
    ///
    /// Only nodes that have been set are yielded.
    pub fn preorder(&self, n: usize) -> Vec<T> {
        let mut out = Vec::new();
        self.preorder_rec(1, n, &mut out);
        out
    }

    fn preorder_rec(&self, idx: usize, n: usize, out: &mut Vec<T>) {
        if idx > n {
            return;
        }
        if let Some(v) = &self.data[idx] {
            out.push(v.clone());
        }
        self.preorder_rec(idx * 2, n, out);
        self.preorder_rec(idx * 2 + 1, n, out);
    }

    /// Post-order traversal (left, right, root) of the nodes at `1..=n`.
    ///
    /// Only nodes that have been set are yielded.
    pub fn postorder(&self, n: usize) -> Vec<T> {
        let mut out = Vec::new();
        self.postorder_rec(1, n, &mut out);
        out
    }

    fn postorder_rec(&self, idx: usize, n: usize, out: &mut Vec<T>) {
        if idx > n {
            return;
        }
        self.postorder_rec(idx * 2, n, out);
        self.postorder_rec(idx * 2 + 1, n, out);
        if let Some(v) = &self.data[idx] {
            out.push(v.clone());
        }
    }

    /// Level-order (breadth-first) traversal of the nodes at `1..=n`.
    ///
    /// Only nodes that have been set are yielded.
    pub fn level_order(&self, n: usize) -> Vec<T> {
        let mut out = Vec::new();
        for idx in 1..=n {
            if let Some(v) = &self.data[idx] {
                out.push(v.clone());
            }
        }
        out
    }
}
