# Library Checker `string`: dsa-templates vs ac-library-rs

All 12 problems of the [Library Checker](https://judge.yosupo.jp/) `string` category are solved
with the Rust crates of this repository and every solution passes the problem's own checker on
every generated test. This page compares run time with
[ac-library-rs](https://github.com/rust-lang-ja/ac-library-rs) (`ac`) wherever ac has something to
build on, and notes the run time of the reference C++ solution on the same test.

## Method

* Same I/O layer for both (`FastScan` / `Out` from `template/rust`), `cargo build --release`.
* For each problem every test is run once per implementation and checked with the problem's
  checker; the heaviest test is then re-run 5 times, alternating implementations. The table
  reports the **minimum** wall-clock time in milliseconds, input parsing and output printing
  included.
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10% (a few percent more
  for the memory-bound `aho_corasick`).
* ac-library-rs offers suffix array, LCP array, Z-algorithm and `convolution_i64`; the ac column is
  filled for the five problems that reduce to them and is **—** otherwise.
* "Reference C++" is the `sol/correct.cpp` shipped with the problem, timed on the same test.

## Summary

Compared against ac on 5 problems: geometric mean of ours/ac = **0.77**; 5 within 10% or faster, 0 slower than that.

New crate `strings`: SA-IS suffix array on `u32` (faster than ac's), Kasai LCP, distinct
substrings, longest common substring, Z-function, Manacher, Duval's Lyndon factorisation,
Aho–Corasick, the eertree, the **double-ended eertree** (push / pop at both ends), runs and
wildcard matching; `seaweed::prefix_substring_lcs`. Every function has a randomized stress test
against brute force (the double-ended eertree against 6000 random operation sequences).

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `aho_corasick` | 5 | `strings` | 316 | — | — | `strings::AhoCorasick`: trie in creation order (the vertex numbering of the statement), breadth-first suffix links on a full `26`-symbol transition table. ac: nothing. Reference C++: 479 ms. |
| `eertree` | 5 | `strings` | 174 | — | — | `strings::Eertree`: online palindromic tree, amortised `O(1)` per character; ids shifted so the roots are `-1` / `0`. ac: nothing. Reference C++: 289 ms. |
| `enumerate_palindromes` | 5 | `strings` | 37 | — | — | `strings::manacher` over the string with separators, the `2n-1` centres in the order of the statement. ac: nothing. Reference C++: 51 ms. |
| `longest_common_substring` | 5 | `strings` | 155 | 186 | 0.83× | `strings::longest_common_substring`: SA-IS and Kasai on `S 0 T`, best LCP between neighbouring suffixes from different strings. Ac column: the same with ac's `suffix_array_arbitrary` / `lcp_array_arbitrary`. Reference C++: 100 ms (suffix automaton). |
| `lyndon_factorization` | 5 | `strings` | 15 | — | — | `strings::lyndon_factorization`: Duval's algorithm, `O(n)`. ac: nothing. Reference C++: 27 ms. |
| `number_of_substrings` | 5 | `strings` | 69 | 83 | 0.83× | `strings::count_distinct_substrings`: `n(n+1)/2 - sum lcp`. Ac column: the same with ac's suffix array. Reference C++: 56 ms. |
| `palindromes_in_deque` | 10 | `strings` | 124 | — | — | `strings::DequeEertree`: eertree with push / pop at both ends via surface occurrences (`pre_surf` / `suf_surf` per position, vertex reference counts) and quick links, `O(log n)` per operation. ac: nothing. Reference C++: 192 ms. |
| `prefix_substring_lcs` | 5 | std | 82 | — | — | `seaweed::prefix_substring_lcs`: seaweeds combed row by row, then for each prefix length one offline dominance count (`#{start >= b, end < c}`) with a Fenwick tree. ac: nothing. Reference C++: 161 ms. |
| `runenumerate` | 5 | `strings` | 135 | — | — | `strings::runs`: longest Lyndon words from the next-smaller-suffix of two suffix arrays (both alphabet orders), extended with suffix-array LCE (`sparse_table`), Bannai et al. ac: nothing. Reference C++: 88 ms. |
| `suffixarray` | 5 | `strings` | 69 | 78 | 0.88× | `strings::suffix_array_bytes`: SA-IS on `u32`. Ac column: `ac_library::suffix_array`. Reference C++: 86 ms. |
| `wildcard_pattern_matching` | 5 | `strings` | 410 | 713 | 0.58× | `strings::wildcard_match`: `sum (s - t)^2 s t` as three exact convolutions (`convolution_u64`, three NTT primes). Ac column: the same with `convolution_i64`. Reference C++: 529 ms. |
| `zalgorithm` | 5 | `strings` | 15 | 19 | 0.79× | `strings::z_algorithm`. Ac column: `ac_library::z_algorithm`. Reference C++: 34 ms. |

