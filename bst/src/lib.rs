//! Binary Search Tree (BST) for ordered integer keys.
//!
//! A simple, non-balancing BST implemented with `Box`-linked nodes.
//! Provides:
//!
//! * [`Bst::new`] — create an empty BST
//! * [`Bst::insert`] — insert a key (duplicates are placed in the right subtree)
//! * [`Bst::contains`] — test membership
//! * [`Bst::inorder`] — return keys in sorted (ascending) order
//! * [`Bst::min`], [`Bst::max`] — read the smallest / largest key
//!
//! This is intentionally non-balancing; it is meant as a teaching and
//! reference snippet, not as a high-performance structure. For practical
//! use prefer [`std::collections::BTreeSet`].

use std::cmp::Ordering;

/// Node of the BST.
#[derive(Debug)]
struct Node {
    key: i64,
    left: Option<Box<Node>>,
    right: Option<Box<Node>>,
}

impl Node {
    fn new(key: i64) -> Self {
        Self { key, left: None, right: None }
    }

    /// Insert `key` into the subtree rooted at `self`. Returns `true` if
    /// the key was newly inserted (i.e. was not already present).
    fn insert(&mut self, key: i64) -> bool {
        match key.cmp(&self.key) {
            Ordering::Less => {
                if let Some(left) = self.left.as_mut() {
                    left.insert(key)
                } else {
                    self.left = Some(Box::new(Node::new(key)));
                    true
                }
            }
            Ordering::Greater => {
                if let Some(right) = self.right.as_mut() {
                    right.insert(key)
                } else {
                    self.right = Some(Box::new(Node::new(key)));
                    true
                }
            }
            Ordering::Equal => false,
        }
    }

    fn contains(&self, key: i64) -> bool {
        match key.cmp(&self.key) {
            Ordering::Less => self.left.as_ref().is_some_and(|n| n.contains(key)),
            Ordering::Greater => self.right.as_ref().is_some_and(|n| n.contains(key)),
            Ordering::Equal => true,
        }
    }

    fn inorder(&self, out: &mut Vec<i64>) {
        if let Some(left) = &self.left {
            left.inorder(out);
        }
        out.push(self.key);
        if let Some(right) = &self.right {
            right.inorder(out);
        }
    }
}

/// Binary search tree holding ordered `i64` keys.
#[derive(Debug, Default)]
pub struct Bst {
    root: Option<Box<Node>>,
}

impl Bst {
    /// Create an empty BST.
    pub fn new() -> Self {
        Self { root: None }
    }

    /// Insert `key` into the tree. Returns `true` if the key was newly
    /// inserted, `false` if it was already present.
    pub fn insert(&mut self, key: i64) -> bool {
        if let Some(root) = self.root.as_mut() {
            root.insert(key)
        } else {
            self.root = Some(Box::new(Node::new(key)));
            true
        }
    }

    /// Test whether `key` is present.
    pub fn contains(&self, key: i64) -> bool {
        self.root.as_ref().is_some_and(|n| n.contains(key))
    }

    /// Return all keys in ascending (sorted) order.
    pub fn inorder(&self) -> Vec<i64> {
        let mut out = Vec::new();
        if let Some(root) = &self.root {
            root.inorder(&mut out);
        }
        out
    }

    /// Return the smallest key, or `None` for an empty tree.
    pub fn min(&self) -> Option<i64> {
        let mut node = self.root.as_deref()?;
        while let Some(left) = node.left.as_deref() {
            node = left;
        }
        Some(node.key)
    }

    /// Return the largest key, or `None` for an empty tree.
    pub fn max(&self) -> Option<i64> {
        let mut node = self.root.as_deref()?;
        while let Some(right) = node.right.as_deref() {
            node = right;
        }
        Some(node.key)
    }
}
