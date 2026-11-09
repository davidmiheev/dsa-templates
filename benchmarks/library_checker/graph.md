# Library Checker `graph`: dsa-templates vs ac-library-rs

All 38 problems of the [Library Checker](https://judge.yosupo.jp/) `graph` category are solved
with the Rust crates of this repository and every solution passes the problem's own checker on
every generated test. This page compares run time with
[ac-library-rs](https://github.com/rust-lang-ja/ac-library-rs) (`ac`) wherever ac has something
to build on.

## Method

* Same I/O layer for both (`FastScan` / `Out` from `template/rust`), `cargo build --release`.
* For each problem every test is run once per implementation and checked with the problem's
  checker; the heaviest test is then re-run 5 times, alternating implementations. The table
  reports the **minimum** wall-clock time of those runs, in milliseconds (reading the input and
  writing the output included).
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10%.
* The only graph tools in ac-library-rs are `Dsu`, `SccGraph`, `MfGraph`, `MinCostFlowGraph` and
  `TwoSat` (plus its segment / Fenwick trees and modular arithmetic). The "ac" column is therefore
  **—** for most rows.
  Where ac's containers can carry the numeric part of the solution, the ac column is the same
  algorithm with a *hand-written* layer on top (flow networks for matching and assignment, a
  divide and conquer on `SccGraph` + `Dsu`, Gaussian elimination on `ModInt998244353`) and the
  comment says so. ac's flow and SCC routines recurse, so those solutions run on a 512 MB stack
  thread, like contest solutions do.
* Every solution is iterative (million-vertex paths are fine) and uses at least one crate of this
  repository.
* The solutions and the benchmark runner are deliberately **not** part of this repository.

## Summary

Compared against ac on 8 problems: geometric mean of ours/ac = **0.51**; 8 within 10% or faster, 0 slower than that.

New or reworked crates driven by this category:

* graph basics: `csr_graph` (CSR adjacency with edge ids), `scc` (iterative Tarjan), `lowlink`
  (bridges, articulation points, 2-edge-connected and biconnected components), `dijkstra`
  (rewritten on CSR, path recovery), `min_span_tree`, `cycle_detection`, `eulerian_trail`,
  `dominator_tree`, `three_edge_cc`, `st_numbering`, `complement_graph`;
* matching and flow: `matching` (Hopcroft-Karp, blossom, Hungarian, weighted blossom),
  `min_cost_flow` (b-flow with bounds and duals), `edge_coloring`, `directed_mst`;
* counting and exact algorithms: `modint`, `modmatrix`, `spanning_count` (matrix-tree and BEST
  theorems), `subgraph_counting` (triangles, cliques, per-edge 4-cycles), `graph_coloring`,
  `independent_set`, `steiner_tree`, `chordal_graph`, `tree_decomposition`, `graph_center`,
  `k_shortest_walks`, `extreme_sets`;
* dynamic and offline: `offline_connectivity` with `union_find::RollbackSumUnionFind`,
  `incremental_scc`, `incremental_msf` (on `link_cut_tree`);
* containers: `hash_map` (open addressing), `deque`.

Every crate has a randomized stress test against a brute-force reference (`cargo test --release
--workspace`). The conventions of the repository hold: structs are generic where it makes sense, and
there is no `Option` in public APIs (preconditions are asserted; "none" results are empty vectors
or sentinels such as `NONE`).

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `assignment` | 5 | `matching` | 171 | 1275 | 0.13× | `hungarian` (potentials, O(n^3)). ac column: min-cost flow on a bipartite network (`MinCostFlowGraph`, costs shifted to be non-negative). |
| `biconnected_components` | 5 | `csr_graph`, `lowlink` | 618 | — | — | `Lowlink::biconnected_components` (iterative DFS, vertex sets, isolated vertices as singletons). ac has no lowlink. |
| `bipartite_edge_coloring` | 10 | `edge_coloring` | 109 | — | — | Regularise (pack vertices into groups of degree <= Delta, pad with dummy edges), Euler-split for even degree, `hopcroft_karp` perfect matching for odd degree. ac has no edge colouring. |
| `bipartitematching` | 5 | `matching` | 439 | 2157 | 0.20× | `hopcroft_karp`. ac column: `MfGraph` max flow on the source/sink network (needs a 512 MB stack thread). |
| `chordal_graph_recognition` | 5 | `chordal_graph` | 86 | — | — | Maximum cardinality search + O(n+m) perfect-elimination check; on failure a shortest path between the two non-adjacent later neighbours closes an induced cycle. ac: nothing. |
| `chromatic_number` | 10 | `graph_coloring` | 99 | — | — | Subset inclusion-exclusion over independent-set counts, O(2^n chi), checked modulo 998244353 and 2^64. ac: nothing. |
| `chromatic_polynomial` | 10 | `graph_coloring` | 1696 | — | — | Per-subset independent-set counts by size, then [y^n] u^i for all i from a cubic recurrence, binomial basis -> monomials; O(2^n n^3) worst case, about 1.7 s at n = 20 (limit 10 s). ac: nothing. |
| `connected_components_of_complement_graph` | 5 | `complement_graph` | 97 | — | — | BFS over the complement with a `FastSet` of unvisited vertices: O((n+m) log_64 n). ac: nothing. |
| `counting_c4` | 10 | `subgraph_counting` | 2259 | — | — | Degree ordering + wedge counting per top vertex, per-edge counts for multigraphs, O(m sqrt m). ac: nothing. |
| `counting_eulerian_circuits` | 5 | `spanning_count` | 81 | 79 | 1.03× | BEST theorem: `count_arborescences` (determinant) times factorials, `union_find` connectivity check. ac column: hand-written Gaussian elimination on `ModInt998244353` plus `Dsu`. |
| `counting_spanning_tree_directed` | 5 | `spanning_count` | 89 | 86 | 1.03× | Matrix-tree theorem, `modmatrix::determinant` (O(n^3)). ac column: the same elimination hand-written on ac's `ModInt998244353`; both are I/O-bound at this size. |
| `counting_spanning_tree_undirected` | 5 | `spanning_count` | 81 | 84 | 0.96× | Matrix-tree theorem, `modmatrix::determinant`. ac column: hand-written elimination on ac's `ModInt998244353`. |
| `cycle_detection` | 5 | `cycle_detection` | 133 | — | — | Iterative three-colour DFS returning edge ids. ac's `SccGraph` finds components but not a cycle. |
| `cycle_detection_undirected` | 5 | `cycle_detection` | 192 | — | — | Iterative DFS skipping the parent *edge* (multigraph safe). ac: nothing. |
| `directedmst` | 5 | `directed_mst` | 69 | — | — | Chu-Liu/Edmonds with skew heaps (lazy subtract tag, iterative merge) and `RollbackUnionFind` for contraction; expansion by rolling back. ac: nothing. |
| `dominatortree` | 5 | `dominator_tree` | 61 | — | — | Lengauer-Tarjan with path-compressed `eval`, fully iterative. ac: nothing. |
| `dynamic_graph_vertex_add_component_sum` | 10 | `offline_connectivity`, `union_find` | 191 | — | — | Offline divide and conquer over time (`offline_dynamic_connectivity`) on `RollbackSumUnionFind` (component sums that survive rollback; point updates walk the union tree). ac's `Dsu` cannot roll back. |
| `enumerate_cliques` | 5 | `subgraph_counting` | 2 | — | — | Degree-ordered clique enumeration. The input is tiny (n, m <= 100): a correctness problem, not a speed one. ac: nothing. |
| `enumerate_triangles` | 5 | `subgraph_counting` | 59 | — | — | Degree-oriented triangle enumeration, O(m sqrt m). ac: nothing. |
| `eulerian_trail_directed` | 5 | `eulerian_trail` | 55 | — | — | Hierholzer on CSR, iterative. ac: nothing. |
| `eulerian_trail_undirected` | 5 | `eulerian_trail` | 51 | — | — | Hierholzer with per-edge used flags, iterative. ac: nothing. |
| `general_matching` | 5 | `matching` | 6 | — | — | `blossom` (Edmonds, O(n^3), BFS with base compression). ac: nothing. |
| `general_weighted_matching` | 5 | `matching` | 538 | — | — | `max_weight_matching`: O(n^3) weighted blossom with dual variables. ac: nothing. |
| `global_minimum_cut_of_dynamic_star_augmented_graph` | 10 | `extreme_sets`, `heavy_light_decomposition`, `lazy_segment_tree` | 1245 | — | — | `extreme_set_tree` (laminar family of extreme sets by repeated contraction), then `Hld` path-to-root add on a `LazySegmentTree` (range add / global min). ac: nothing. |
| `incremental_minimum_spanning_forest` | 10 | `incremental_msf` | 6065 | — | — | `IncrementalMsf`: `LinkCutTree` over vertices plus one node per edge, path-max query decides which edge leaves. Same time as the reference C++ solution on the same machine (about 6 s on the largest test, limit 10 s). ac has no link-cut tree. |
| `incremental_scc` | 5 | `incremental_scc`, `union_find` | 1094 | 3913 | 0.28× | `incremental_scc_times`: offline divide and conquer over time using `scc`, then `union_find` merges in time order maintaining sum x_i x_j. ac column: the same divide and conquer hand-written on `SccGraph` + `Dsu`. |
| `k_shortest_walk` | 5 | `k_shortest_walks` | 223 | — | — | Eppstein: reverse `dijkstra`, persistent leftist heaps of sidetrack edges merged along the shortest-path tree, best-first expansion. ac: nothing. |
| `maximum_independent_set` | 5 | `independent_set` | 2 | — | — | Bitmask branching (take degree <= 1 vertices, else branch on the maximum degree vertex). n <= 40, a few milliseconds. ac: nothing. |
| `min_cost_b_flow` | 10 | `min_cost_flow` | 13 | — | — | `min_cost_b_flow`: capacity-scaling successive shortest paths with lower bounds, negative costs and dual potentials; passes the `anti_ssp` test that defeats plain SSP. ac's `MinCostFlowGraph` returns no potentials. |
| `minimum_diameter_spanning_tree` | 5 | `graph_center` | 354 | — | — | n Dijkstra runs, absolute centre by the skyline sweep over every edge, shortest-path tree from the (possibly mid-edge) centre. ac: nothing. |
| `minimum_spanning_tree` | 5 | `min_span_tree` | 116 | 144 | 0.81× | Kruskal on `union_find`. ac column: sort + `Dsu`. |
| `minimum_steiner_tree` | 5 | `steiner_tree` | 137 | — | — | Dreyfus-Wagner with edge recovery, O(3^k n + 2^k m log m). ac: nothing. |
| `scc` | 5 | `scc` | 154 | 210 | 0.73× | Iterative Tarjan with components in topological order. ac column: `SccGraph::scc` (recursive, needs a 512 MB stack). |
| `shortest_path` | 5 | `csr_graph`, `dijkstra` | 327 | — | — | `dijkstra` with parent edges and path recovery. ac has no shortest-path routine. |
| `st_numbering` | 5 | `st_numbering` | 93 | — | — | Tarjan's lowpoint list insertion (iterative), biconnectivity checked on the way. ac: nothing. |
| `three_edge_connected_components` | 5 | `three_edge_cc` | 110 | — | — | Randomised cycle-space hashing on a DFS tree (fixed seed): classes of equal edge signatures, fingerprint per vertex. ac: nothing. |
| `tree_decomposition_width_2` | 10 | `tree_decomposition` | 740 | — | — | Series-parallel reduction (delete degree <= 2 vertices, fill in an edge) with a `hash_map::FastMap` edge set; bags hang below the first-eliminated neighbour. ac: nothing. |
| `two_edge_connected_components` | 5 | `csr_graph`, `lowlink` | 79 | — | — | `Lowlink::two_edge_connected_components`. ac: nothing. |

