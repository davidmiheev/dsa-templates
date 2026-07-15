Centroid Decomposition
======================

Centroid Decomposition of a tree.

Centroid decomposition recursively splits a tree at its centroid — the
vertex whose removal leaves components of size at most ``n/2`` — and
repeats on each component. The resulting recursion tree has height
O(log n), and any two original vertices share an ancestor in the
decomposition at depth at most O(log n).

This is the foundation for "small-to-large" tricks on trees, for
answering distance queries offline (count paths of length ≤ k), and for
many divide-and-conquer-on-tree problems.

This crate exposes ``build_centroid_tree``, which returns the centroid
tree: for each vertex, the centroid chosen at the next level up (its
"centroid parent") and the depth of the vertex in the centroid tree.

The Rust crate lives at ``centroid_decomposition/src/lib.rs``; build & run it with:

.. code-block:: console

    $ cd centroid_decomposition && cargo run
    $ cd centroid_decomposition && cargo test

.. rust-api:: centroid_decomposition
   :description: Generated API reference for the Centroid Decomposition crate.

.. literalinclude:: ../../centroid_decomposition/src/lib.rs
   :language: rust
   :linenos: