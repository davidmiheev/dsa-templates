Bit Manipulation
================

Bit manipulation utilities and a tiny fixed-size bitset.

The free-standing functions cover the operations that show up over and
over in competitive programming: ``popcount``, ``trailing_zeros`` /
``trailing_ones``, ``highest_bit`` / ``lowest_bit``, ``is_power_of_two``,
``next_power_of_two`` / ``log2_floor``, and enumeration of all subsets
of a bit mask.

The ``Bitset`` type wraps ``Vec<u64>`` and exposes the operations you'd
expect: ``get``, ``set``, ``reset``, ``flip``, plus ``count_ones`` over
the whole set.

The Rust crate lives at ``bit_manipulation/src/lib.rs``; build & run it with:

.. code-block:: console

    $ cd bit_manipulation && cargo run
    $ cd bit_manipulation && cargo test

.. rust-api:: bit_manipulation
   :description: Generated API reference for the bit manipulation crate.

.. literalinclude:: ../../bit_manipulation/src/lib.rs
   :language: rust
   :linenos: