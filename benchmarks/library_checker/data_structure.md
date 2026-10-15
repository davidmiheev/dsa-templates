# Library Checker `data_structure`: dsa-templates vs ac-library-rs

All 47 problems of the [Library Checker](https://judge.yosupo.jp/) `data_structure` category are
solved with the Rust crates of this repository and every solution passes the problem's own
checker on every generated test. This page compares run time with
[ac-library-rs](https://github.com/rust-lang-ja/ac-library-rs) (`ac`) wherever ac has something
to build on.

## Method

* Same I/O layer for both (`FastScan` / `Out` from `template/rust`), same algorithm family where
  possible, `cargo build --release`.
* For each problem every test is run once per implementation and checked with the problem's
  checker; the heaviest test is then re-run 5 times, alternating implementations. The table
  reports the **minimum** wall-clock time of those runs, in milliseconds (includes reading the
  input and writing the output).
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10%.
* "ac" column legend: **—** means ac offers nothing to build the solution on (no ordered set, no
  persistence, no 2D / tree structure, ...); where a workaround exists it is described in the
  comment (offline coordinate compression, a hand-written CDQ on ac's `FenwickTree`, ...). Those
  workarounds are *offline*, ours are mostly *online*, which is why a few rows are slower.
* The solutions and the benchmark runner are deliberately **not** part of this repository.

## Summary

Compared against ac on 25 problems: geometric mean of ours/ac = **0.91**; 20 within 10% or faster, 5 slower than that (`persistent_queue` 1.50×, `point_set_range_composite_large_array` 2.76×, `range_affine_range_sum` 1.13×, `range_affine_range_sum_large_array` 2.86×, `static_range_count_distinct` 1.12×).

New or reworked crates driven by this category: `algebra`, `lazy_segment_tree`, `segment_tree`
(non-commutative, iterative), `sparse_table` (block-decomposed), `sorted_list`, `union_find`
(plain / potential / rollback / range-parallel), `fenwick` (linear build, `Fenwick2D`,
`dominance_sums`), `li_chao_tree`, `binary_trie`, `fast_set`, `queue_aggregation`,
`wavelet_matrix`, `implicit_treap`, `persistent_queue`, `persistent_lazy_segment_tree`,
`dynamic_segment_tree`, `point_set_range_frequency`, `segment_tree_beats`,
`linear_add_range_min`, `kd_tree`, `range_sort_sequence`, `seaweed`, `mo_algorithm`
(rollback Mo, second-order Mo), `hash_map` (open-addressing map) and `deque` (ring buffer). Each has a randomized stress test against a naive model.

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `area_of_union_of_rectangles` | 5 | `lazy_segment_tree` | 1252 | 1206 | 1.04× | Sweep with a lazy segtree over (min, count-of-min). Parity with ac's `LazySegtree`. |
| `associative_array` | 5 | `hash_map` | 281 | — | — | Open-addressing `FastMap` (linear probing, multiplicative hash, backward-shift deletion). ac has no hash map. |
| `deque` | 5 | `deque` | 28 | — | — | Power-of-two ring buffer `Deque` with O(1) indexing. ac has no deque. |
| `deque_operate_all_composite` | 5 | `queue_aggregation` | 37 | — | — | Two-stack deque aggregation (`DequeAggregation`). ac has no such structure. |
| `double_ended_priority_queue` | 10 | `multiset` | 339 | — | — | Counting `multiset`. ac has no ordered containers. |
| `dynamic_point_set_rectangle_affine_rectangle_sum` | 10 | `kd_tree` | 3627 | — | — | kd-tree with lazy affine tags (leaf buckets of 32). Slowest test is an adversarial case from a Library Checker issue. ac has no 2D structure. |
| `dynamic_sequence_range_affine_range_sum` | 10 | `implicit_treap` | 3178 | — | — | Implicit treap (insert / erase / lazy affine). ac has no balanced sequence. |
| `line_add_get_min` | 5 | `li_chao_tree` | 125 | — | — | Dynamic Li Chao tree. ac has none. |
| `majority_voting` | 5 | `point_set_range_frequency`, `segment_tree` | 658 | 801 | 0.82× | Boyer-Moore monoid on `SegmentTree` + `PointSetRangeFrequency` to verify. ac column: ac `Segtree` plus hand-built per-value `FenwickTree`s. |
| `ordered_set` | 10 | `sorted_list` | 220 | 649 | 0.34× | Bucketed `SortedList`. ac has no ordered set: its column is an offline workaround (compress values, 0/1 presence `Segtree`, `max_right` descent). |
| `persistent_queue` | 5 | `persistent_queue` | 156 | 104 | 1.50× | Online skew-binary level-ancestor queue. ac column is a std-only *offline* DFS over the version tree with a `VecDeque` (ac has nothing persistent), so it is not like for like. |
| `persistent_range_affine_range_sum` | 5 | `persistent_lazy_segment_tree` | 569 | — | — | Path-copying lazy segtree with `copy_range` between versions. ac has no persistence. |
| `persistent_unionfind` | 5 | `union_find` (rollback) | 41 | — | — | Offline: version tree DFS with `RollbackUnionFind`. ac's `Dsu` cannot undo. |
| `point_add_range_sum` | 5 | `fenwick` | 86 | 96 | 0.90× | `FenwickTree` vs ac `FenwickTree`. |
| `point_add_rectangle_sum` | 5 | `fenwick` | 243 | 277 | 0.88× | `fenwick::dominance_sums` (CDQ). ac column: the same CDQ written by hand on ac's `FenwickTree` (ac has no 2D structure). |
| `point_set_range_composite` | 5 | `segment_tree` | 253 | 257 | 0.99× | `SegmentTree` (non-commutative) vs ac `Segtree`. |
| `point_set_range_composite_large_array` | 10 | `dynamic_segment_tree` | 860 | 312 | 2.76× | Online `DynamicSegmentTree` over 10^9 slots. ac column: offline coordinate compression + `Segtree`; the online pointer-based tree is ~2.8x slower (cache misses on ~190 MB of nodes). |
| `point_set_range_frequency` | 10 | `point_set_range_frequency` | 121 | 198 | 0.61× | `PointSetRangeFrequency` (one Fenwick per value in flat arrays). ac column: per-value ac `FenwickTree`s built by hand. |
| `point_set_range_sort_range_composite` | 10 | `range_sort_sequence` | 1501 | — | — | `RangeSortSequence`: sorted runs as mergeable segment trees + position segtree. ac has nothing comparable. |
| `predecessor_problem` | 10 | `fast_set` | 137 | — | — | 64-ary `FastSet`. ac has no ordered set. |
| `queue_operate_all_composite` | 5 | `queue_aggregation` | 37 | — | — | Two-stack queue aggregation. ac has no such structure. |
| `range_add_range_min` | 5 | `lazy_segment_tree` | 307 | 310 | 0.99× | `LazySegmentTree` vs ac `LazySegtree`: parity. |
| `range_affine_point_get` | 5 | `lazy_segment_tree` | 541 | 505 | 1.07× | `LazySegmentTree` vs ac `LazySegtree`: parity. |
| `range_affine_range_sum` | 5 | `lazy_segment_tree` | 662 | 585 | 1.13× | `LazySegmentTree` vs ac `LazySegtree`: within ~13%; ac's `apply_range` avoids some pushes we do. |
| `range_affine_range_sum_large_array` | 10 | `dynamic_segment_tree` | 507 | 177 | 2.86× | Online `DynamicLazySegmentTree` over 10^9 slots. ac column: offline compression of query boundaries + `LazySegtree`; ~2.9x slower for the same pointer-chasing reason. |
| `range_chmin_chmax_add_range_sum` | 10 | `segment_tree_beats` | 400 | — | — | `SegmentTreeBeats`. ac's `LazySegtree` cannot express chmin/chmax with sums. |
| `range_kth_smallest` | 5 | `wavelet_matrix` | 178 | — | — | `WaveletMatrix::kth_smallest`. ac has none. |
| `range_linear_add_range_min` | 10 | `linear_add_range_min` | 1310 | — | — | Sqrt decomposition with a convex hull per block (`LinearAddRangeMin`). ac has none. |
| `range_parallel_unionfind` | 5 | `union_find` | 822 | 1204 | 0.68× | `RangeParallelUnionFind` (doubling levels with lazy range links). ac column: the same doubling idea with one ac `Dsu` per level. |
| `range_reverse_range_sum` | 5 | `implicit_treap` | 651 | — | — | `ImplicitTreap` with reversal. ac has no balanced sequence. |
| `range_set_range_composite` | 5 | `lazy_segment_tree` | 1562 | 1738 | 0.90× | `LazySegmentTree` vs ac `LazySegtree`: parity. |
| `rectangle_add_point_get` | 10 | `fenwick` | 321 | 383 | 0.84× | `fenwick::dominance_sums` (CDQ). ac column: same CDQ by hand on ac's `FenwickTree`. |
| `rectangle_sum` | 5 | `fenwick` | 237 | 259 | 0.92× | Offline sweep with `FenwickTree` vs the same sweep with ac's. |
| `segment_add_get_min` | 5 | `li_chao_tree` | 1386 | — | — | `LiChaoTree::add_segment`. ac has none. |
| `set_xor_min` | 5 | `binary_trie` | 397 | — | — | `BinaryTrie`. ac has none. |
| `static_range_count_distinct` | 5 | `fenwick` | 266 | 236 | 1.12× | Offline sweep + `FenwickTree` vs ac `FenwickTree`: ~12% slower in this run, but run-to-run noise on this machine is about that large (earlier runs showed parity). |
| `static_range_frequency` | 5 | `wavelet_matrix` | 586 | — | — | `WaveletMatrix::range_freq`. ac has none. |
| `static_range_inversions_query` | 5 | `mo_algorithm` | 139 | 750 | 0.18× | Second-order Mo (`range_inversions`): O(1) work per pointer move. ac column: plain Mo with an ac `FenwickTree` (O(log n) per move), ~5x slower. |
| `static_range_lis_query` | 5 | `seaweed` | 721 | — | — | Seaweed / unit-Monge product (`seaweed::range_lis`), O((n+q) log^2 n). ac has nothing for this. |
| `static_range_mode_query` | 5 | `mo_algorithm` | 115 | 112 | 1.03× | Rollback Mo (`solve_rollback`). ac column is hand-written rollback Mo (no ac component applies): the lib wrapper costs nothing. |
| `static_range_sum` | 5 | `fenwick` | 66 | 72 | 0.92× | `FenwickTree::new` (O(n) build) and prefix differences; ac column is `FenwickTree` filled by n `add` calls. Both are dominated by I/O. |
| `static_range_sum_with_upper_bound` | 5 | `fenwick` | 345 | 381 | 0.91× | Offline `FenwickTree` vs ac `FenwickTree`; a `WaveletMatrix` variant (`..._wm`) was also checked. |
| `static_rectangle_add_rectangle_sum` | 10 | `fenwick` | 1125 | 1070 | 1.05× | Offline sweep with a 4-channel `FenwickTree` over `Mod4`; ac column is the same sweep on ac's `FenwickTree`. |
| `staticrmq` | 5 | `sparse_table` | 134 | 163 | 0.82× | Block-decomposed `SparseTable` (O(1) query, ~3n memory) vs ac `Segtree` (`prod`). |
| `unionfind` | 5 | `union_find` | 12 | 16 | 0.73× | `union_find::plain` (single i32 array, iterative find) vs ac `Dsu`. |
| `unionfind_with_potential` | 5 | `union_find` | 27 | — | — | `PotentialUnionFind` over a group. ac's `Dsu` has no potentials. |
| `unionfind_with_potential_non_commutative_group` | 5 | `union_find` | 74 | — | — | `PotentialUnionFind` over a non-commutative group (2x2 matrices). ac has none. |

