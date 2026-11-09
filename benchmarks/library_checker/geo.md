# Library Checker `geo`: dsa-templates vs ac-library-rs

All 10 problems of the [Library Checker](https://judge.yosupo.jp/) `geo` category are solved with
the Rust crates of this repository and every solution passes the problem's own checker on every
generated test. ac-library-rs has **no geometry routines**, so the comparison here is against the
reference C++ solution that ships with each problem (`sol/correct.cpp`, timed on the same test and
machine), quoted in the comments.

## Method

* Same I/O layer (`FastScan` / `Out` from `template/rust`), `cargo build --release`.
* For each problem every test is run once and checked with the problem's checker; the heaviest
  test is then re-run 5 times. The table reports the **minimum** wall-clock time in milliseconds,
  input parsing and output printing included.
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10%.

## Summary

No ac column (ac-library-rs has nothing for this category). Against the reference C++ solutions: within 10% or faster on 8 of 10 problems; geometric mean of ours / reference = **0.75**.

New crate `geometry` (exact integer predicates throughout): convex hull, farthest pair, all
farthest neighbours of a convex polygon, closest pair, argument sort, minimum enclosing circle,
Manhattan and Euclidean MST, convex layers (decremental hull) and `TriangleCounter`. Each has a
randomized stress test against brute force.

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `all_furthest_neighbors` | 5 | `geometry` | 130 | — | — | `geometry::all_farthest_neighbors`: the doubled distance matrix of a convex polygon satisfies the quadrangle inequality, so the leftmost row maximum is monotone and a divide and conquer finds all of them in `O(n log n)`. ac: nothing. Reference C++: 310 ms. |
| `closest_pair` | 5 | `geometry` | 139 | — | — | `geometry::closest_pair`: divide and conquer on x with a merge by y and a strip scan, exact `u128` distances. ac: nothing. Reference C++: 240 ms. |
| `convex_layers` | 5 | `geometry` | 1905 | — | — | `geometry::convex_layers`: Overmars-van Leeuwen decremental hull (a segment tree of bridges, `O(log^2 n)` per deletion) for the left chain and, on the reflected points, the right chain; layers keep the points on hull edges. ac: nothing. Reference C++: 2.05 s (the same data structure). |
| `count_points_in_triangle` | 5 | `geometry` | 316 | — | — | `geometry::TriangleCounter`: for every ordered pair of `A` the points of `B` strictly left of / on the line within the open horizontal strip, plus the points on the horizontal lines of the vertices; a triangle is a signed combination of three such counts, `O(N^2 M + Q)`. ac: nothing. Reference C++: 281 ms (same method). |
| `euclidean_mst` | 5 | `geometry` | 403 | — | — | `geometry::euclidean_mst`: Boruvka rounds with a kd-tree (points stored contiguously in tree order, subtree component flags prune whole boxes, the component's best edge seeds each search). Reference C++: 244 ms (Delaunay triangulation). ac: nothing. |
| `furthest_pair` | 5 | `geometry` | 104 | — | — | `geometry::farthest_pair`: monotone-chain hull and rotating calipers. ac: nothing. Reference C++: 180 ms. |
| `manhattanmst` | 5 | `geometry` | 294 | — | — | `geometry::manhattan_mst`: four sweeps over the octants with an ordered map give `O(n)` candidate edges, Kruskal with `union_find`. ac: nothing. Reference C++: 521 ms. |
| `minimum_enclosing_circle` | 5 | `geometry` | 17 | — | — | `geometry::minimum_enclosing_circle`: randomized incremental construction with exact predicates (in-circle determinant, Thales test); the points on the final circle are those with a zero predicate. ac: nothing. Reference C++: 29 ms. |
| `sort_points_by_argument` | 5 | `geometry` | 53 | — | — | `geometry::argument_sort`: four classes (lower half plane, positive axis, upper half plane, negative axis), exact cross products inside a half plane. ac: nothing. Reference C++: 69 ms. |
| `static_convex_hull` | 5 | `geometry` | 185 | — | — | `geometry::convex_hull`: Andrew's monotone chain with strict turns, duplicates dropped. ac: nothing. Reference C++: 215 ms. |

