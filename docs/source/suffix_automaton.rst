Suffix Automaton
================

Suffix Automaton (SAM).

A suffix automaton is a directed acyclic automaton that recognises
exactly the set of substrings of a given string. It has at most
``2n - 1`` states for an input of length ``n``, and every state
corresponds to a distinct set of end positions in the original string —
hence to a distinct set of substrings.

Once built, you can answer:

* Is ``p`` a substring of ``s``? (O(|p|))
* How many distinct substrings does ``s`` have? (sum over states of
  ``len(state) - len(link(state))``)
* Number of occurrences of a pattern (when ``state.occ`` is computed)
* Longest common substring between two strings

The Rust crate lives at ``suffix_automaton/src/lib.rs``; build & run it with:

.. code-block:: console

    $ cd suffix_automaton && cargo run
    $ cd suffix_automaton && cargo test

.. rust-api:: suffix_automaton
   :description: Generated API reference for the Suffix Automaton crate.

.. literalinclude:: ../../suffix_automaton/src/lib.rs
   :language: rust
   :linenos: