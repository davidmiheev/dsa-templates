Sparse Table
============

Sparse Table for idempotent range queries.

A Sparse Table answers range queries on a static array in O(1) query
time after O(n log n) preprocessing, *as long as* the query operation
is idempotent and associative. The classic examples are ``min``,
``max``, and ``gcd``; sum and xor work too.

The trick: precompute ``st[k][i] = op(a[i], a[i+1], ..., a[i+2^k-1])``
and overlap two intervals of length ``2^k`` to cover any range.

The Rust crate lives at ``sparse_table/src/lib.rs``; build & run it with:

.. code-block:: console

    $ cd sparse_table && cargo run
    $ cd sparse_table && cargo test

.. rust-api:: sparse_table
   :description: Generated API reference for the Sparse Table crate.

.. literalinclude:: ../../sparse_table/src/lib.rs
   :language: rust
   :linenos: