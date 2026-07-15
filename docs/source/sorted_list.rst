Sorted List
============

A treap-based sorted list data structure.

Supports:
- O(log n) insert, remove, bisect_left, bisect_right, get
- O(1) len, is_empty

Uses a flat array representation (indices instead of boxed nodes) with a free stack for memory reuse.

Rust Implementation
-------------------

.. rust-api:: sorted_list
   :description: Generated API reference for the Sorted List crate.

.. literalinclude:: ../../sorted_list/src/lib.rs
   :language: rust
   :linenos: