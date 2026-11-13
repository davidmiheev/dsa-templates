Graphs
======

Traversals, shortest paths, connectivity, spanning trees, colouring, flows and matchings.

.. toctree::

   graph
   bfs
   dfs
   dijkstra
   minimum_spanning_tree
   topological_sort

CSR Graph
---------

.. rust-api:: csr_graph
   :description: Compressed sparse row adjacency in two flat arrays.

Strongly Connected Components
-----------------------------

.. rust-api:: scc
   :description: Strongly connected components by an iterative Tarjan search.

Lowlink: Bridges and Articulation Points
----------------------------------------

.. rust-api:: lowlink
   :description: Lowlink values of an undirected multigraph: bridges, articulation points, biconnected components.

3-Edge-Connected Components
---------------------------

.. rust-api:: three_edge_cc
   :description: 3-edge-connected components of an undirected multigraph.

st-Numbering
------------

.. rust-api:: st_numbering
   :description: st-numbering of a graph (Tarjan).

Cycle Detection
---------------

.. rust-api:: cycle_detection
   :description: Find a simple cycle in a directed or undirected multigraph.

Eulerian Trail
--------------

.. rust-api:: eulerian_trail
   :description: Eulerian trails and circuits by an iterative Hierholzer search.

2-SAT
-----

.. rust-api:: two_sat
   :description: 2-SAT via strongly connected components.

Dominator Tree
--------------

.. rust-api:: dominator_tree
   :description: Dominator tree of a directed graph (Lengauer-Tarjan).

Directed MST
------------

.. rust-api:: directed_mst
   :description: Minimum-cost arborescence (directed MST) by Chu-Liu/Edmonds.

Incremental MSF
---------------

.. rust-api:: incremental_msf
   :description: Online minimum spanning forest under edge insertions.

Incremental SCC
---------------

.. rust-api:: incremental_scc
   :description: Offline incremental strongly connected components.

Offline Dynamic Connectivity
----------------------------

.. rust-api:: offline_connectivity
   :description: Offline dynamic connectivity.

Complement Graph Components
---------------------------

.. rust-api:: complement_graph
   :description: Connected components of a graph's complement without building it.

Chordal Graphs
--------------

.. rust-api:: chordal_graph
   :description: Chordal graphs and perfect elimination orderings.

Small Subgraph Counting
-----------------------

.. rust-api:: subgraph_counting
   :description: Small-subgraph enumeration on sparse graphs by degree ordering.

Graph Coloring
--------------

.. rust-api:: graph_coloring
   :description: Exact graph colouring for small graphs by inclusion-exclusion.

Bipartite Edge Coloring
-----------------------

.. rust-api:: edge_coloring
   :description: Optimal edge colouring of a bipartite multigraph.

Independent Set and Clique
--------------------------

.. rust-api:: independent_set
   :description: Maximum independent set and maximum clique on graphs with at most 64 vertices.

Steiner Tree
------------

.. rust-api:: steiner_tree
   :description: Minimum Steiner tree (Dreyfus-Wagner).

Graph Center
------------

.. rust-api:: graph_center
   :description: Minimum diameter spanning tree of a weighted graph (Hassin-Tamir).

Extreme Sets
------------

.. rust-api:: extreme_sets
   :description: Laminar tree of extreme sets of a weighted undirected graph.

Matching
--------

.. rust-api:: matching
   :description: Matching: Hopcroft-Karp, blossom, Hungarian, weighted blossom.

k Shortest Walks
----------------

.. rust-api:: k_shortest_walks
   :description: The ``k`` shortest walks between two vertices of a directed graph.

Min-Cost Flow
-------------

.. rust-api:: min_cost_flow
   :description: Minimum-cost b-flow with bounds, negative costs and dual potentials.
