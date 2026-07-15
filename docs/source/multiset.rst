MultiSet
========

A multiset is a set that allows for multiple occurrences of the same element. It is similar to a list, but it does not allow for duplicate elements.
We will use multiset that has keys in sorted order

Rust Implementation
-------------------
In Rust, we can implement a multiset using a `BTreeMap` to store the elements and their counts.
The `BTreeMap` will maintain the sorted order of the elements, allowing us to efficiently get the smallest and largest elements.
The `BTreeMap` isn't prone to collisions blowup (`HashMap` is prone), so it is a good choice for implementing a multiset.

.. rust-api:: multiset
   :description: Generated API reference for the Multiset crate.

.. literalinclude:: ../../multiset/src/lib.rs
   :language: rust
   :linenos: