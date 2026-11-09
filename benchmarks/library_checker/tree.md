# Library Checker `tree`: dsa-templates vs ac-library-rs

All 20 problems of the [Library Checker](https://judge.yosupo.jp/) `tree` category are solved
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
* ac-library-rs has **no tree support at all** (no LCA, HLD, centroid, link-cut, ...). The "ac"
  column is therefore **—** for most rows. Where ac's `Segtree`, `FenwickTree` or `convolution`
  can carry the numeric part of the solution, the ac column is the same algorithm with a
  *hand-written* tree layer (HLD, Euler preorder, centroid decomposition) on top of ac's
  containers, and the comment says so; those rows compare our tree layer + containers with
  hand-rolled code + ac's containers.
* The solutions and the benchmark runner are deliberately **not** part of this repository.

## Summary

Compared against ac on 5 problems: geometric mean of ours/ac = **0.61**; 5 within 10% or faster, 0 slower than that.

New or reworked crates driven by this category: `heavy_light_decomposition` (iterative rewrite,
CSR input, `lca` / `jump` / `path_vertex`, allocation-free ordered path pieces),
`centroid_decomposition` (iterative `CentroidTree` with ancestor distances, `ContourSum`,
`ContourAdd`, `pair_counts_by_distance`), `tree_algos` (Cartesian tree, diameter, rooted
isomorphism classes, rerooting, tree-constrained scheduling), `link_cut_tree`,
`euler_tour_tree`, `tree_path_affine_sum`, `permutation_tree`, `convolution` (NTT). Each has a
randomized stress test against a naive model.

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `cartesian_tree` | 5 | `tree_algos` | 76 | — | — | Stack construction (`cartesian_tree`), I/O bound. Nothing in ac applies. |
| `common_interval_decomposition_tree` | 5 | `permutation_tree` | 609 | — | — | `PermutationTree`: monotone stack of strong intervals + lazy segment tree finding the nearest common interval, O(n log n). ac could supply the `LazySegtree` but the algorithm is not in ac; no ac solution was written. |
| `dynamic_tree_subtree_add_subtree_sum` | 5 | `euler_tour_tree` | 2295 | — | — | `EulerTourTree` with a lazy whole-component add (the hard case for link-cut trees). ac has none. |
| `dynamic_tree_vertex_add_path_sum` | 5 | `link_cut_tree` | 748 | — | — | `LinkCutTree` over a sum monoid. ac has no link-cut tree. |
| `dynamic_tree_vertex_add_subtree_sum` | 5 | `euler_tour_tree` | 2284 | — | — | `EulerTourTree`: cut the edge, read the component, link back. ac has none. |
| `dynamic_tree_vertex_set_path_composite` | 5 | `link_cut_tree` | 1060 | — | — | `LinkCutTree` over a non-commutative affine monoid (reversal-aware). ac has none. |
| `frequency_table_of_tree_distance` | 10 | `centroid_decomposition`, `convolution` | 2234 | 2271 | 0.98× | `CentroidTree::pair_counts_by_distance` driving our `convolution_u64`. ac column: a hand-written centroid decomposition driving ac's `convolution_i64`. |
| `jump_on_tree` | 5 | `lca` | 849 | — | — | `lca::Lca` for the LCA and `lca::LevelAncestor` (jump pointers + ladders) for an `O(1)` k-th ancestor. Same run, heaviest test: 769 ms vs 1160 ms for the previous `Hld::path_vertex` solution and 2655 ms for binary lifting; summed over all tests 10.2 s vs 13.7 s (HLD) vs 18.5 s (binary lifting). Nothing in ac applies. |
| `lca` | 5 | `lca` | 278 | 379 | 0.73× | `lca::Lca`: `O(1)` query from a sparse table over blocks of 16 preorder slots plus in-block prefix / suffix minima, built by sequential passes over a BFS relabelling (`O(n)` memory). Same run, heaviest test: 312 ms vs 625 ms for the previous `Hld::lca` solution and 559 ms for binary lifting (`lca::BinaryLifting`); summed over all 25 tests 5.6 s vs 8.4 s (HLD) vs 18.0 s (binary lifting, 3-5x slower on deep trees: every step is a cache miss). ac has no tree support: the ac column is an Euler-preorder + ac `Segtree` range-max formulation. |
| `point_set_tree_path_composite_sum` | 10 | `tree_path_affine_sum` | 2107 | — | — | `TreePathAffineSum::sum_from(r)`: dynamic DP over heavy chains plus a second segment tree for the "outside the subtree" maps; O(log^2 n) per query for any root. ac has none. |
| `point_set_tree_path_composite_sum_fixed_root` | 10 | `tree_path_affine_sum` | 1436 | — | — | Same structure, `sum_from_root()` only (no outside part). ac has none. |
| `rooted_tree_isomorphism_classification` | 5 | `tree_algos` | 403 | — | — | Bottom-up canonical ids (sorted child-id vectors in a `HashMap`). Nothing in ac applies. |
| `rooted_tree_topological_order_with_minimum_inversions` | 5 | `tree_algos` | 374 | — | — | `min_inversion_order`: repeatedly merge the group with the smallest D/C ratio into its parent's group (heap + linked lists). Nothing in ac applies. |
| `tree_diameter` | 5 | `tree_algos` | 169 | — | — | Two BFS sweeps (`tree_diameter`). Nothing in ac applies. |
| `tree_path_composite_sum` | 5 | `tree_algos` (rerooting) | 186 | — | — | Generic `rerooting` (prefix/suffix merges, no inverse needed). Nothing in ac applies. |
| `vertex_add_path_sum` | 5 | `lca`, `fenwick` | 451 | 1372 | 0.33× | Preorder numbering of `lca::Lca`: a vertex value is added on its subtree range of a `FenwickTree`, so the prefix sum at `tin[v]` is the root path sum and a path is `P(u) + P(v) - 2 P(lca) + a[lca]` (three prefix sums, `O(1)` LCA). Same run, heaviest test: 496 ms vs 988 ms for the previous HLD + path-segment solution; summed over all tests 5.5 s vs 9.7 s. ac column: hand-written HLD (same algorithm) + ac `FenwickTree`. |
| `vertex_add_range_contour_sum_on_tree` | 10 | `centroid_decomposition` | 935 | — | — | `ContourSum`: centroid decomposition with two Fenwick trees per centroid. ac has no tree support. |
| `vertex_add_subtree_sum` | 5 | `lca`, `fenwick` | 263 | 452 | 0.58× | Preorder subtree ranges of `lca::Lca` + `FenwickTree` (no HLD needed). Same run, heaviest test: 316 ms vs 492 ms for the previous HLD-based solution; summed over all tests 2.4 s vs 3.9 s. ac column: hand-written HLD + ac `FenwickTree`. |
| `vertex_get_range_contour_add_on_tree` | 10 | `centroid_decomposition` | 561 | — | — | `ContourAdd`: dual of the above (range add by distance, point read). ac has no tree support. |
| `vertex_set_path_composite` | 5 | `heavy_light_decomposition`, `segment_tree` | 391 | 623 | 0.63× | HLD ordered path pieces (`for_each_path_segment_ordered`) over a `SegmentTree` of (forward, backward) compositions; HLD stays because the path composite changes with point updates, which binary lifting or an `O(1)` LCA cannot maintain. ac column: hand-written HLD + ac `Segtree`. |

