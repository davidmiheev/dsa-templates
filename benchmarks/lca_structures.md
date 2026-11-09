# LCA problems: HLD vs binary lifting vs `O(1)` structures

Library Checker problems that need a lowest common ancestor or a k-th ancestor, solved with the
heavy-light decomposition (`heavy_light_decomposition::Hld`, the previous solutions), textbook
binary lifting and the new [`lca`](../lca) crate. All three give identical answers on every test
(checked against the expected outputs); times are wall-clock milliseconds with input parsing and
output printing included, minimum of 3-7 runs, measured in one session on one machine.

| Problem | what is asked | HLD (before) | `lca` crate | binary lifting |
|---|---|---:|---:|---:|
| `lca` | 5e5 LCA queries on a given parent array | 8425 | **5599** | 17985 |
| `jump_on_tree` | k-th vertex on a path | 13732 | **10232** | 18463 |
| `vertex_add_path_sum` | point add, path sum | 9743 | **5473** | not applicable |
| `vertex_add_subtree_sum` | point add, subtree sum | 3921 | **2363** | not applicable |

Sums of the per-test minima over all generated tests of the problem (25 / 21 / 19 / 16 tests).
On the single heaviest test of each problem:

| Problem (test) | HLD | `lca` crate | binary lifting |
|---|---:|---:|---:|
| `lca` (`max_random_02`) | 625 | **312** | 559 |
| `jump_on_tree` (`almost_line_00`) | 1160 | **769** | 2655 |
| `vertex_add_path_sum` (`max_random_02`) | 988 | **496** | |
| `vertex_add_subtree_sum` (`max_random_01`) | 492 | **316** | |

## What the data says

* **Binary lifting does not beat HLD here.** A query touches `log n` entries of a 40 MB table
  (`n = 5e5`, 19 levels), each a cache miss, so it is 3-5x slower than HLD on deep trees (lines,
  almost-lines) and only on par on shallow random trees. Its one advantage is simplicity.
* **HLD is great on paths and lines** (a single chain, almost no hops) and weak on bushy trees:
  a complete binary tree makes a query hop through ~18 chains. The `O(1)` structure is ~2x faster
  there (`binary_*`, `max_random_*`: 520-750 ms vs 250-330 ms) and ties on lines.
* **`O(1)` LCA** ([`lca::Lca`](../lca/src/lib.rs)): a sparse table over blocks of 16 preorder slots
  plus in-block prefix / suffix minima of "preorder number of the parent" (the LCA of `u != v` is the
  parent with the smallest number on `(tin[u], tin[v]]`). `O(n)` memory, the table fits in cache.
* **`O(1)` level ancestor** ([`lca::LevelAncestor`](../lca/src/lib.rs)): jump pointers + ladders of the
  long-path decomposition (one jump-pointer read, one ladder read). With the block LCA it is 25-30%
  faster than HLD on `jump_on_tree`.
* **Path sums need no decomposition**: add a vertex value on its subtree range of a Fenwick tree over
  the preorder, then `P(u) + P(v) - 2 P(lca) + a[lca]` is three prefix sums plus an `O(1)` LCA (the
  HLD solution does two prefix sums per chain piece). `vertex_add_subtree_sum` needs only the
  preorder ranges, not an HLD.
* **BFS relabelling**: every table is built over BFS labels (ancestors have smaller labels, the
  parents of consecutive labels never decrease), so the jump-pointer and sparse tables are filled by
  sequential passes even when the input labels are random; queries translate the original labels
  with one array read.
* **HLD stays where it is needed**: `vertex_set_path_composite` (non-commutative path composite with
  point updates) cannot use a static `O(1)` LCA or binary lifting - the aggregate changes with every
  update - so it keeps `Hld::for_each_path_segment_ordered` over a segment tree.

Reproduce: `lc-bench` solutions `*_hld.rs` / `*_ours.rs` / `*_bl.rs` for each problem (not in the
repository), `cmp.py tree <problem> <binaries...>`.
