Interval Tree
========================

An interval tree is a tree data structure that holds intervals and
supports stabbing queries — given a point ``x``, find all intervals
that contain ``x`` — in ``O(log n + k)`` time, where ``k`` is the
number of reported intervals.

The canonical implementation is the *augmented* interval tree: a
binary search tree keyed on interval start, where each node stores
the maximum endpoint in its subtree. A stabbing query visits only
those subtrees whose maximum endpoint is at least ``x``.

The crate here exposes a builder that constructs the tree from a
list of intervals and a query method that returns the indices (or
any associated payload) of all matching intervals.

Rust Implementation
-------------------

The Rust crate lives at ``interval_tree/src/lib.rs``; build & run it
with:

.. code-block:: console

    $ cd interval_tree && cargo run
    $ cd interval_tree && cargo test

.. rust-api:: interval_tree
   :description: Generated API reference for the Interval Tree crate.

.. literalinclude:: ../../interval_tree/src/lib.rs
   :language: rust
   :linenos:

Applications
------------

* Stabbing queries on event intervals
* Overlap detection in scheduling problems
* Window-based aggregations (e.g. "which bookings are active now?")