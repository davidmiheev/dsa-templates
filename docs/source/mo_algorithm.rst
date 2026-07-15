Mo's Algorithm
==============

Mo's algorithm for offline range queries.

Mo's algorithm solves queries of the form "given a value array and a
list of range queries, answer each query using only O(|hi - lo|) moves
between consecutive queries (with O(1) amortized per move) once the
queries have been reordered by Mo's block ordering.

The classic instance is counting distinct elements in a range. This
crate ships ``count_distinct_in_range`` which handles that case
directly. For more elaborate query types build your own solver using
``solve_blocked``.

The Rust crate lives at ``mo_algorithm/src/lib.rs``; build & run it with:

.. code-block:: console

    $ cd mo_algorithm && cargo run
    $ cd mo_algorithm && cargo test

.. rust-api:: mo_algorithm
   :description: Generated API reference for the Mo's algorithm crate.

.. literalinclude:: ../../mo_algorithm/src/lib.rs
   :language: rust
   :linenos: