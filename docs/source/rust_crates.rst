Rust Crates
===========

Every crate of the Rust workspace that has no dedicated topic page, grouped by area.
Each entry links to the ``cargo doc`` API reference; the source is available there through the *source* links.

Data structures
---------------

Sequences, range queries, persistent and specialised containers.

algebra
~~~~~~~

.. rust-api:: algebra
   :description: Monoid, group and action traits shared by the data-structure crates.

lazy_segment_tree
~~~~~~~~~~~~~~~~~

.. rust-api:: lazy_segment_tree
   :description: Lazy segment tree with range update and range query.

persistent_lazy_segment_tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rust-api:: persistent_lazy_segment_tree
   :description: Fully persistent lazy segment tree.

dynamic_segment_tree
~~~~~~~~~~~~~~~~~~~~

.. rust-api:: dynamic_segment_tree
   :description: Segment trees over a huge index range that allocate only the nodes they touch.

segment_tree_beats
~~~~~~~~~~~~~~~~~~

.. rust-api:: segment_tree_beats
   :description: Segment tree beats (Ji driver segment tree) over ``i64``.

linear_add_range_min
~~~~~~~~~~~~~~~~~~~~

.. rust-api:: linear_add_range_min
   :description: Range add of an arithmetic progression and range minimum.

range_sort_sequence
~~~~~~~~~~~~~~~~~~~

.. rust-api:: range_sort_sequence
   :description: Sequence supporting range sort by key, point set and range monoid product.

point_set_range_frequency
~~~~~~~~~~~~~~~~~~~~~~~~~

.. rust-api:: point_set_range_frequency
   :description: Range frequency queries with point assignments.

wavelet_matrix
~~~~~~~~~~~~~~

.. rust-api:: wavelet_matrix
   :description: Wavelet matrix: order statistics and frequencies on ranges of a static array.

implicit_treap
~~~~~~~~~~~~~~

.. rust-api:: implicit_treap
   :description: Implicit treap: a dynamic sequence with ``O(log n)`` expected time per operation.

kd_tree
~~~~~~~

.. rust-api:: kd_tree
   :description: 2D kd-tree over a fixed point set with a value per point.

li_chao_tree
~~~~~~~~~~~~

.. rust-api:: li_chao_tree
   :description: Li Chao tree: minimum of a set of lines or segments at integer points.

binary_trie
~~~~~~~~~~~

.. rust-api:: binary_trie
   :description: Binary trie multiset over fixed-width unsigned integers, for xor queries.

fast_set
~~~~~~~~

.. rust-api:: fast_set
   :description: Set of integers in ``0..n`` as a 64-ary bitset tree.

queue_aggregation
~~~~~~~~~~~~~~~~~

.. rust-api:: queue_aggregation
   :description: Queue and deque that keep the fold of a monoid over their elements.

persistent_queue
~~~~~~~~~~~~~~~~

.. rust-api:: persistent_queue
   :description: Fully persistent queue.

deque
~~~~~

.. rust-api:: deque
   :description: Growable ring-buffer deque with a power-of-two capacity.

hash_map
~~~~~~~~

.. rust-api:: hash_map
   :description: ``FastMap``: fast hash map from ``u64`` keys.

binary_tree
~~~~~~~~~~~

.. rust-api:: binary_tree
   :description: Generic binary tree with iterative traversals.

bst
~~~

.. rust-api:: bst
   :description: Binary search tree for ordered integer keys.

Trees
-----

LCA, dynamic forests and tree decompositions.

lca
~~~

.. rust-api:: lca
   :description: Lowest common ancestor and level ancestor of a rooted tree.

tree_algos
~~~~~~~~~~

.. rust-api:: tree_algos
   :description: Small self-contained tree algorithms.

link_cut_tree
~~~~~~~~~~~~~

.. rust-api:: link_cut_tree
   :description: Link-cut tree over a forest with monoid-valued vertices.

euler_tour_tree
~~~~~~~~~~~~~~~

.. rust-api:: euler_tour_tree
   :description: Euler tour tree: a dynamic forest with monoid-valued vertices.

tree_path_affine_sum
~~~~~~~~~~~~~~~~~~~~

.. rust-api:: tree_path_affine_sum
   :description: Sum over all vertices of a function composed along the path from a chosen root.

permutation_tree
~~~~~~~~~~~~~~~~

.. rust-api:: permutation_tree
   :description: Permutation tree (common interval decomposition tree).

tree_decomposition
~~~~~~~~~~~~~~~~~~

.. rust-api:: tree_decomposition
   :description: Tree decomposition of width at most 2.

Graphs
------

Connectivity, orderings, spanning trees, colouring, flows and matchings.

csr_graph
~~~~~~~~~

.. rust-api:: csr_graph
   :description: Compressed sparse row adjacency in two flat arrays.

scc
~~~

.. rust-api:: scc
   :description: Strongly connected components by an iterative Tarjan search.

lowlink
~~~~~~~

.. rust-api:: lowlink
   :description: Lowlink values of an undirected multigraph: bridges, articulation points, biconnected components.

three_edge_cc
~~~~~~~~~~~~~

.. rust-api:: three_edge_cc
   :description: 3-edge-connected components of an undirected multigraph.

st_numbering
~~~~~~~~~~~~

.. rust-api:: st_numbering
   :description: st-numbering of a graph (Tarjan).

top_sort
~~~~~~~~

.. rust-api:: top_sort
   :description: Topological sort (Kahn's algorithm).

cycle_detection
~~~~~~~~~~~~~~~

.. rust-api:: cycle_detection
   :description: Find a simple cycle in a directed or undirected multigraph.

eulerian_trail
~~~~~~~~~~~~~~

.. rust-api:: eulerian_trail
   :description: Eulerian trails and circuits by an iterative Hierholzer search.

two_sat
~~~~~~~

.. rust-api:: two_sat
   :description: 2-SAT via strongly connected components.

dominator_tree
~~~~~~~~~~~~~~

.. rust-api:: dominator_tree
   :description: Dominator tree of a directed graph (Lengauer-Tarjan).

directed_mst
~~~~~~~~~~~~

.. rust-api:: directed_mst
   :description: Minimum-cost arborescence (directed MST) by Chu-Liu/Edmonds.

min_span_tree
~~~~~~~~~~~~~

.. rust-api:: min_span_tree
   :description: Kruskal's minimum spanning forest.

incremental_msf
~~~~~~~~~~~~~~~

.. rust-api:: incremental_msf
   :description: Online minimum spanning forest under edge insertions.

incremental_scc
~~~~~~~~~~~~~~~

.. rust-api:: incremental_scc
   :description: Offline incremental strongly connected components.

offline_connectivity
~~~~~~~~~~~~~~~~~~~~

.. rust-api:: offline_connectivity
   :description: Offline dynamic connectivity.

complement_graph
~~~~~~~~~~~~~~~~

.. rust-api:: complement_graph
   :description: Connected components of a graph's complement without building it.

chordal_graph
~~~~~~~~~~~~~

.. rust-api:: chordal_graph
   :description: Chordal graphs and perfect elimination orderings.

subgraph_counting
~~~~~~~~~~~~~~~~~

.. rust-api:: subgraph_counting
   :description: Small-subgraph enumeration on sparse graphs by degree ordering.

graph_coloring
~~~~~~~~~~~~~~

.. rust-api:: graph_coloring
   :description: Exact graph colouring for small graphs by inclusion-exclusion.

edge_coloring
~~~~~~~~~~~~~

.. rust-api:: edge_coloring
   :description: Optimal edge colouring of a bipartite multigraph.

independent_set
~~~~~~~~~~~~~~~

.. rust-api:: independent_set
   :description: Maximum independent set and maximum clique on graphs with at most 64 vertices.

steiner_tree
~~~~~~~~~~~~

.. rust-api:: steiner_tree
   :description: Minimum Steiner tree (Dreyfus-Wagner).

graph_center
~~~~~~~~~~~~

.. rust-api:: graph_center
   :description: Minimum diameter spanning tree of a weighted graph (Hassin-Tamir).

extreme_sets
~~~~~~~~~~~~

.. rust-api:: extreme_sets
   :description: Laminar tree of extreme sets of a weighted undirected graph.

matching
~~~~~~~~

.. rust-api:: matching
   :description: Matching: Hopcroft-Karp, blossom, Hungarian, weighted blossom.

k_shortest_walks
~~~~~~~~~~~~~~~~

.. rust-api:: k_shortest_walks
   :description: The ``k`` shortest walks between two vertices of a directed graph.

min_cost_flow
~~~~~~~~~~~~~

.. rust-api:: min_cost_flow
   :description: Minimum-cost b-flow with bounds, negative costs and dual potentials.

Strings and sequences
---------------------

Suffix structures, palindromes, subsequences.

strings
~~~~~~~

.. rust-api:: strings
   :description: String algorithms: suffix array, Z-function, Manacher, Aho-Corasick, eertree, runs.

subsequence
~~~~~~~~~~~

.. rust-api:: subsequence
   :description: Longest increasing subsequence in ``O(n log n)``.

seaweed
~~~~~~~

.. rust-api:: seaweed
   :description: Seaweed (sticky braid) machinery and LIS of subarrays.

Number theory and combinatorics
-------------------------------

Modular arithmetic, sieves, big integers, counting.

modint
~~~~~~

.. rust-api:: modint
   :description: ``StaticModInt``: integers modulo a compile-time constant.

number_theory
~~~~~~~~~~~~~

.. rust-api:: number_theory
   :description: 64-bit number theory: Miller-Rabin, Pollard's rho, primitive roots.

sublinear
~~~~~~~~~

.. rust-api:: sublinear
   :description: Sums over floor quotients in sublinear time.

prime_sieve
~~~~~~~~~~~

.. rust-api:: prime_sieve
   :description: Sieves of Eratosthenes.

stern_brocot
~~~~~~~~~~~~

.. rust-api:: stern_brocot
   :description: The Stern-Brocot tree of positive rationals.

bigint
~~~~~~

.. rust-api:: bigint
   :description: Arbitrary-precision integers.

combinatorics
~~~~~~~~~~~~~

.. rust-api:: combinatorics
   :description: Enumerative combinatorics.

divisor_transform
~~~~~~~~~~~~~~~~~

.. rust-api:: divisor_transform
   :description: Zeta and Mobius transforms along the divisibility order.

gf2_64
~~~~~~

.. rust-api:: gf2_64
   :description: Arithmetic and fast polynomial multiplication in ``GF(2^64)``.

Polynomials, convolution and linear algebra
-------------------------------------------

NTT-based convolutions, power series, matrices and recurrences.

convolution
~~~~~~~~~~~

.. rust-api:: convolution
   :description: Convolutions by the number-theoretic transform.

polynomial
~~~~~~~~~~

.. rust-api:: polynomial
   :description: Polynomials and formal power series over an NTT-friendly prime field.

poly_mod_p
~~~~~~~~~~

.. rust-api:: poly_mod_p
   :description: Polynomials over ``F_p`` for a runtime prime.

set_power_series
~~~~~~~~~~~~~~~~

.. rust-api:: set_power_series
   :description: Subset transforms, subset convolution, exp, log and composition.

multiplicative_convolution
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rust-api:: multiplicative_convolution
   :description: Convolutions where the indices multiply.

multivariate_convolution
~~~~~~~~~~~~~~~~~~~~~~~~

.. rust-api:: multivariate_convolution
   :description: Products of multivariate polynomials stored as dense arrays.

min_plus
~~~~~~~~

.. rust-api:: min_plus
   :description: Min-plus convolution of convex, concave and arbitrary sequences.

linear_recurrence
~~~~~~~~~~~~~~~~~

.. rust-api:: linear_recurrence
   :description: Linear recurrences over ``Z/MZ``.

modmatrix
~~~~~~~~~

.. rust-api:: modmatrix
   :description: Dense matrices over ``Z/MZ``.

bit_matrix
~~~~~~~~~~

.. rust-api:: bit_matrix
   :description: Matrices over ``F_2`` packed into ``u64`` words.

spanning_count
~~~~~~~~~~~~~~

.. rust-api:: spanning_count
   :description: Counting spanning structures with determinants.

Geometry
--------

geometry
~~~~~~~~

.. rust-api:: geometry
   :description: Computational geometry on integer points with exact predicates.
