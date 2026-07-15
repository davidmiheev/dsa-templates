Rolling Hash
============

Rolling (Rabin-Karp-style) hash for strings and byte slices.

Provides two flavours:

* ``RollingHash`` — single-modulus rolling hash over a byte slice. The
  current hash value is maintained incrementally as the window slides,
  giving O(1) updates per shift.
* ``DoubleRollingHash`` — same idea, with two different moduli for very
  low collision rates.

Both flavours precompute the base-power table for substrings up to some
``max_len``. To compare hashes of arbitrary substrings of a long string,
use ``compute_hash`` (single) or build a ``DoubleRollingHash`` and query
substrings using ``hash_range(lo, hi)``.

The Rust crate lives at ``rolling_hash/src/lib.rs``; build & run it with:

.. code-block:: console

    $ cd rolling_hash && cargo run
    $ cd rolling_hash && cargo test

.. rust-api:: rolling_hash
   :description: Generated API reference for the rolling hash crate.

.. literalinclude:: ../../rolling_hash/src/lib.rs
   :language: rust
   :linenos: