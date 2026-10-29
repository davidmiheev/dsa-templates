# Library Checker `other`: dsa-templates vs ac-library-rs

All 7 problems of the [Library Checker](https://judge.yosupo.jp/) `other` category are solved with
the Rust crates of this repository and every solution passes the problem's own checker on every
generated test. ac-library-rs only has a 2-SAT, so the comparison is against the reference C++
solution that ships with each problem (`sol/correct.cpp`, timed on the same test and machine),
quoted in the comments; for `two_sat` the ac column is shown too.

## Method

* Same I/O layer (`FastScan` / `Out` from `template/rust`), `cargo build --release`.
* For each problem every test is run once and checked with the problem's checker; the heaviest
  test is then re-run 5 times. The table reports the **minimum** wall-clock time in milliseconds,
  input parsing and output printing included.
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10%.

## Summary

ac-library-rs only has `two_sat` (where it overflows the stack on one test) and nothing for the rest. Against the reference C++ solutions: within 10% or faster on 6 of 7 problems; geometric mean of ours / reference = **0.64**.

> **Note on the timings.** They were measured with the portable radix-4 transforms. The `convolution` crate now runs its NTT on AVX2 when the CPU has it (identical results): the NTT-bound solutions of this table are 1.4-2.6x faster, see [`../ntt_engines.md`](../ntt_engines.md) for the profile, the engines compared and per-problem before / after numbers.

New code: `linear_recurrence::{kth_term, consecutive_terms}` (spectral Bostan–Mori; the terms at a
huge index via `x^k mod g`), `combinatorics::{sum_exponential_times_polynomial,
sum_exponential_times_polynomial_limit}`, and the new crates `two_sat` and `subsequence`
(longest increasing subsequence). Each has a randomized stress test against brute force.

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `consecutive_terms_of_linear_recurrent_sequence` | 10 | `linear_recurrence`, `modint` | 1905 | — | — | `linear_recurrence::consecutive_terms`: `x^k mod g` by left-to-right squaring (a squaring, a reciprocal-based reduction, and a free shift for each bit), then the `M` terms as one correlation of the remainder with the first `2d` terms of the series and a final power-series quotient for the terms past `d`. ac: nothing. Reference C++: 2.21 s (Bostan–Mori with polynomial output). |
| `find_linear_recurrence` | 5 | `linear_recurrence`, `modint` | 52 | — | — | `linear_recurrence::berlekamp_massey`, `O(N^2)`. ac: nothing. Reference C++: 199 ms. |
| `kth_term_of_linearly_recurrent_sequence` | 10 | `linear_recurrence`, `modint` | 640 | — | — | `linear_recurrence::kth_term`: Bostan–Mori on spectra, two forward transforms of size `2^18` and two inverse of size `2^17` per bit; `Q(-x)` is a neighbour swap in the bit-reversed spectrum and the even/odd part is read off at the squared points, so no full-size inverse transform is needed. ac: nothing. Reference C++: 400 ms (a hand-tuned NTT; ours is 1.6x slower here, the one problem of the group where we lose). |
| `longest_increasing_subsequence` | 5 | `subsequence` | 46 | — | — | `subsequence::longest_increasing_subsequence`: patience sorting with `partition_point` and back pointers, `O(n log n)`. ac: nothing. Reference C++: 67 ms. |
| `sum_of_exponential_times_polynomial` | 5 | `combinatorics` | 517 | — | — | `combinatorics::sum_exponential_times_polynomial`: `S(n) = r^n F(n) - c` with `F` of degree `d`; the constant `c` comes from the vanishing `(d+1)`-th difference of `(S(j) + c) r^{-j}`, `F` is then interpolated at `n` in `O(d)` (`r = 0, 1` special cases); the powers `i^d` come from a linear sieve. ac: nothing. Reference C++: 778 ms. |
| `sum_of_exponential_times_polynomial_limit` | 5 | `combinatorics` | 433 | — | — | `combinatorics::sum_exponential_times_polynomial_limit`: the same constant, the limit is `-c`. ac: nothing. Reference C++: 605 ms. |
| `two_sat` | 5 | `two_sat` | 256 | FAIL | — | `two_sat::TwoSat` on `scc::Scc` (an iterative SCC; assignment from the topological order of the components). ac: `ac_library::TwoSat` runs the same time but its recursive SCC **overflows the 8 MB stack** on the long-chain test `cycle_unsat_00` (FAIL; the timing shown is the heaviest test that ran). Reference C++: 693 ms. |

