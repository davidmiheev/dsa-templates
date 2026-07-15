Breadth-First Search
====================

Breadth-First Search (BFS) on a graph.

BFS visits neighbours of the current node before going deeper. It is
typically used to compute shortest paths on unweighted graphs, to find
connected components, or to perform level-order traversal.

The Rust crate lives at ``bfs/src/lib.rs``; build & run it with:

.. code-block:: console

    $ cd bfs && cargo run
    $ cd bfs && cargo test

.. rust-api:: bfs
   :description: Generated API reference for the Breadth-First Search crate.

.. literalinclude:: ../../bfs/src/lib.rs
   :language: rust
   :linenos: