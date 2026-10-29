# Library Checker `enumerative_combinatorics`: dsa-templates vs ac-library-rs

All 17 problems of the [Library Checker](https://judge.yosupo.jp/) `enumerative_combinatorics`
category are solved with the Rust crates of this repository and every solution passes the
problem's own checker on every generated test. This page compares run time with
[ac-library-rs](https://github.com/rust-lang-ja/ac-library-rs) (`ac`) wherever ac has something
to build on, and notes the run time of the reference C++ solution on the same test.

## Method

* Same I/O layer for both (`FastScan` / `Out` from `template/rust`), `cargo build --release`.
* For each problem every test is run once per implementation and checked with the problem's
  checker; the heaviest test is then re-run 3 times, alternating implementations. The table
  reports the **minimum** wall-clock time in milliseconds, input parsing and output printing
  included.
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10%.
* ac-library-rs has no combinatorics or power-series routines. The ac column is **—** except for
  four problems that reduce to Newton iterations / convolutions, where it is *the same algorithm
  hand-written on top of ac's `convolution`* (`lcbench::acpoly`), isolating the cost of the
  convolution layer.
* "Reference C++" is the `sol/correct.cpp` shipped with the problem, timed on the same test.

## Summary

Compared against ac on 4 problems: geometric mean of ours/ac = **0.65**; 4 within 10% or faster, 0 slower than that.

> **Note on the timings.** They were measured with the portable radix-4 transforms. The `convolution` crate now runs its NTT on AVX2 when the CPU has it (identical results): the NTT-bound solutions of this table are 1.4-2.6x faster, see [`../ntt_engines.md`](../ntt_engines.md) for the profile, the engines compared and per-problem before / after numbers.

New crate `combinatorics`: `BlockFactorials` (huge factorials), `BinomialMod` / `BinomialPrime` /
`QBinomialPrime`, Stirling numbers of both kinds (rows, fixed `k`, and the small-prime / huge-`n`
queries), Bell and partition numbers, subset-sum counts, derangements modulo any integer,
distinct subsequences and `count_increasing_sequences` (lattice paths between two staircases by
divide and conquer). Each has a randomized stress test against dynamic programming or brute
force.

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `bell_number` | 10 | `combinatorics` | 614 | 1020 | 0.60× | `combinatorics::bell_numbers`: `B_k = k! [x^k] exp(e^x - 1)`, one Newton exponential. ac column: the same with ac's `convolution` (`lcbench::acpoly`). Reference C++: 188 ms (a faster online exp). |
| `binomial_coefficient` | 10 | `combinatorics` | 614 | — | — | `BinomialMod`: factor `m`; for every prime power `p^e` the factorial without `p` (table of products of integers coprime to `p`, `+-1` per full period, recursion on `n / p`) and the valuation, then CRT. Reference C++: 238 ms. |
| `binomial_coefficient_prime_mod` | 10 | `combinatorics` | 281 | — | — | `BinomialPrime`: factorial and inverse-factorial tables up to the largest `n` of the input, Barrett multiplication for the runtime modulus. ac: nothing. Reference C++: 421 ms. |
| `factorial` | 10 | `combinatorics` | 184 | — | — | `BlockFactorials`: values of `prod_{i=1}^{1024} (1024 t + i)` at `t = 0..M/1024` by degree doubling and chunked sample shifts, prefix products give `(1024 j)!`, then up to 512 multiplications. ac: nothing. Reference C++: 61 ms (a tuned `O(sqrt(M) log M)` single-query algorithm). |
| `many_factorials` | 10 | `combinatorics` | 229 | — | — | Same table, 1e5 queries of at most 512 multiplications each. ac: nothing. Reference C++: 527 ms. |
| `montmort_number_mod` | 10 | `combinatorics` | 46 | — | — | `derangements_mod`: `a_k = k a_{k-1} + (-1)^k`, no division, any modulus. ac: nothing. |
| `number_of_increasing_sequences_between_two_sequences` | 10 | `combinatorics` | 945 | — | — | `count_increasing_sequences`: lattice paths between the staircases of the prefix maxima of `a` and suffix minima of `b`; the grid is cut along its longer side down to 32x32 cells, blocks inside the region are crossed with four binomial convolutions, `O((n+m) log^2)`. ac: nothing. Reference C++: 1500 ms. |
| `number_of_subsequences` | 5 | `combinatorics` | 35 | — | — | `distinct_subsequences`: `dp_i = 2 dp_{i-1} - dp_{last(a_i) - 1}` with `hash_map::FastMap` for the last occurrences. ac: nothing. Reference C++: 489 ms (ordered map). |
| `partition_function` | 10 | `combinatorics` | 236 | 253 | 0.93× | `partition_numbers`: Newton inverse of the pentagonal-number series. Ac column: the same with ac's `convolution`. Reference C++: 842 ms. |
| `q_binomial_coefficient_prime_mod` | 10 | `combinatorics` | 553 | — | — | `QBinomialPrime`: order `d` of `q` from the factorisation of `m - 1`, tables of q-factorials below `d`, q-Lucas `C(n/d, k/d) [n%d, k%d]_q`. ac: nothing. Reference C++: 620 ms. |
| `sharp_p_subset_sum` | 10 | `combinatorics` | 538 | 827 | 0.65× | `subset_sum_counts`: `exp(sum_s c_s log(1 + x^s))`, the logarithm assembled in `O(T log T)` from the multiplicities. Ac column: the same with ac's `convolution`. Reference C++: 955 ms. |
| `stirling_number_of_the_first_kind` | 10 | `combinatorics` | 148 | — | — | `stirling1_row`: `P_{2m}(x) = P_m(x) P_m(x - m)` with a Taylor shift, `O(n log n)`. ac: nothing. Reference C++: 185 ms. |
| `stirling_number_of_the_first_kind_fixed_k` | 10 | `combinatorics` | 779 | — | — | `stirling1_fixed_k`: `s(n, k) = n!/k! [x^n] log(1+x)^k`, one `pow` of `log(1+x)/x`. ac: nothing. Reference C++: 1177 ms. |
| `stirling_number_of_the_first_kind_small_p_large_n` | 10 | `combinatorics` | 229 | — | — | `stirling1_small_prime`: `x(x-1)...(x-n+1) = (x^p - x)^a x(x-1)...(x-b+1)` modulo `p` for `n = a p + b`, so at most two Lucas binomials times a row entry; queries grouped by `b` while the Stirling row is advanced, `O(p^2 + T log_p n)`. ac: nothing. Reference C++: 166 ms. |
| `stirling_number_of_the_second_kind` | 10 | `combinatorics` | 96 | 191 | 0.50× | `stirling2_row`: one convolution of `j^n / j!` (powers by a linear sieve) with `(-1)^i / i!`. Ac column: the same with ac's `convolution`. Reference C++: 173 ms. |
| `stirling_number_of_the_second_kind_fixed_k` | 10 | `combinatorics` | 891 | — | — | `stirling2_fixed_k`: `S(n, k) = n!/k! [x^n] (e^x - 1)^k`, one `pow` of `(e^x - 1)/x`. ac: nothing. Reference C++: 1270 ms. |
| `stirling_number_of_the_second_kind_small_p_large_n` | 10 | `combinatorics` | 132 | — | — | `stirling2_small_prime`: with `k = c p + r`, `1 / prod_{j<=k}(1 - j x) = (1 - x^{p-1})^{-c} / prod_{j<=r}(1 - j x)` modulo `p`, so one coefficient of `prod_{j=r+1}^{p-1}(1 - j x)` (periodic) times one Lucas binomial; queries grouped by `r`. ac: nothing. Reference C++: 159 ms. |

