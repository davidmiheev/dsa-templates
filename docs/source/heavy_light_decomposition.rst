Heavy-Light Decomposition
=========================

Heavy-Light Decomposition (HLD).

Heavy-Light Decomposition splits a tree into "heavy" paths (chains of
vertices where each parent has its largest subtree as its first child)
and "light" edges (connecting chains). Once decomposed, any path from
``u`` to ``v`` can be split into O(log n) segments, each lying within a
single chain. By storing per-chain data in a Fenwick tree or segment
tree, you can answer path queries / updates in O(log² n) overall.

This crate ships the decomposition itself: it computes for each vertex
its parent, depth, subtree size, the head of the chain it belongs to,
and its position in the linearised chain order. You can then plug in
any range-query data structure of your choice on top of that
linearisation.

The Rust crate lives at ``heavy_light_decomposition/src/lib.rs``; build & run it with:

.. code-block:: console

    $ cd heavy_light_decomposition && cargo run
    $ cd heavy_light_decomposition && cargo test

.. rust-api:: heavy_light_decomposition
   :description: Generated API reference for the Heavy-Light Decomposition crate.

.. literalinclude:: ../../heavy_light_decomposition/src/lib.rs
   :language: rust
   :linenos: