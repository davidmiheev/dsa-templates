Depth-First Search
==================

Depth-First Search (DFS) on a graph.

DFS explores as deep as possible along each branch before backtracking.
It is the natural building block for topological order, cycle detection
on directed/undirected graphs, strongly connected components, and tree
traversals.

The Rust crate lives at ``dfs/src/lib.rs``; build & run it with:

.. code-block:: console

    $ cd dfs && cargo run
    $ cd dfs && cargo test

.. rust-api:: dfs
   :description: Generated API reference for the Depth-First Search crate.

.. literalinclude:: ../../dfs/src/lib.rs
   :language: rust
   :linenos: