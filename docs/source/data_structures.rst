Data Structures
===============

Containers, range-query structures and persistent structures.

.. toctree::

   stack
   queue
   multiset
   sorted_list
   binary_tree
   binary_search_tree
   trie
   union_find
   fenwick_tree
   segment_tree
   sparse_table
   interval_tree

Algebra Traits
--------------

.. rust-api:: algebra
   :description: Monoid, group and action traits shared by the data-structure crates.

Lazy Segment Tree
-----------------

.. rust-api:: lazy_segment_tree
   :description: Lazy segment tree with range update and range query.

Persistent Lazy Segment Tree
----------------------------

.. rust-api:: persistent_lazy_segment_tree
   :description: Fully persistent lazy segment tree.

Dynamic Segment Tree
--------------------

.. rust-api:: dynamic_segment_tree
   :description: Segment trees over a huge index range that allocate only the nodes they touch.

Segment Tree Beats
------------------

.. rust-api:: segment_tree_beats
   :description: Segment tree beats (Ji driver segment tree) over ``i64``.

Linear Add, Range Min
---------------------

.. rust-api:: linear_add_range_min
   :description: Range add of an arithmetic progression and range minimum.

Range Sort Sequence
-------------------

.. rust-api:: range_sort_sequence
   :description: Sequence supporting range sort by key, point set and range monoid product.

Point Set, Range Frequency
--------------------------

.. rust-api:: point_set_range_frequency
   :description: Range frequency queries with point assignments.

Wavelet Matrix
--------------

.. rust-api:: wavelet_matrix
   :description: Wavelet matrix: order statistics and frequencies on ranges of a static array.

Implicit Treap
--------------

.. rust-api:: implicit_treap
   :description: Implicit treap: a dynamic sequence with ``O(log n)`` expected time per operation.

k-d Tree
--------

.. rust-api:: kd_tree
   :description: 2D kd-tree over a fixed point set with a value per point.

Li Chao Tree
------------

.. rust-api:: li_chao_tree
   :description: Li Chao tree: minimum of a set of lines or segments at integer points.

Binary Trie
-----------

.. rust-api:: binary_trie
   :description: Binary trie multiset over fixed-width unsigned integers, for xor queries.

Fast Set
--------

.. rust-api:: fast_set
   :description: Set of integers in ``0..n`` as a 64-ary bitset tree.

Queue Aggregation
-----------------

.. rust-api:: queue_aggregation
   :description: Queue and deque that keep the fold of a monoid over their elements.

Persistent Queue
----------------

.. rust-api:: persistent_queue
   :description: Fully persistent queue.

Ring-Buffer Deque
-----------------

.. rust-api:: deque
   :description: Growable ring-buffer deque with a power-of-two capacity.

Fast Hash Map
-------------

.. rust-api:: hash_map
   :description: ``FastMap``: fast hash map from ``u64`` keys.
