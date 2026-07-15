Union Find (Disjoint Set)
=========================

Union Find is a data structure that keeps track of elements which are partitioned into disjoint sets (for example, connected components of a graph).
It provides operations to union existing sets, and find the set to which a particular element belongs.
Each set is represented by a unique identifier, and this identifier is root of tree with all elements in the set as its vertices.
It is commonly used in network connectivity, image processing, and clustering algorithms.

Rust Implementation
-------------------

.. rust-api:: union_find
   :title: Union Find
   :description: Generated API reference for the Union Find crate.

.. literalinclude:: ../../union_find/src/lib.rs
   :language: rust
   :linenos:

Python Implementation
---------------------
.. automodule:: union_find.union_find
    :members:
    :undoc-members:
    :show-inheritance:
    :inherited-members: