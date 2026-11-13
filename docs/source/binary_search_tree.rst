Binary Search Tree
==================

A binary search tree (BST) keeps its keys ordered: every key in the left subtree of a node is smaller than the node's key, and every key in the right subtree is not smaller.
Search, insert and delete walk a single root-to-leaf path, so they take time proportional to the height of the tree, ``O(log n)`` on average and ``O(n)`` for a degenerate (sorted-input) tree.
The implementation here is a simple, non-balancing BST over integer keys.

Rust Implementation
-------------------

.. rust-api:: bst
   :description: Generated API reference for the Binary Search Tree crate.

.. literalinclude:: ../../bst/src/lib.rs
   :language: rust
   :linenos:
