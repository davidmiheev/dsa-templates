Dijkstra's Algorithm
====================

Dijkstra's single-source shortest path on a weighted graph.

Given a directed or undirected graph with non-negative edge weights,
Dijkstra's algorithm computes the shortest distance from a source vertex
to every other reachable vertex in ``O((V + E) log V)`` time using a
binary heap.

The Rust crate lives at ``dijkstra/src/lib.rs``; build & run it with:

.. code-block:: console

    $ cd dijkstra && cargo run
    $ cd dijkstra && cargo test

.. rust-api:: dijkstra
   :description: Generated API reference for the Dijkstra crate.

.. literalinclude:: ../../dijkstra/src/lib.rs
   :language: rust
   :linenos: