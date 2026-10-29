# NTT / FFT: profiling and engine comparison

Every polynomial, power-series, combinatorics and big-integer solution is a pile of number-theoretic
transforms (NTT), so the transform is where the time goes. This note records the profile, the engines
that were written and compared, and what ended up in the `convolution` crate.

**Result:** the `convolution` crate now runs its transforms on AVX2 (eight Montgomery lanes, detected at
run time, bit-identical output to the portable code). Over the 66 NTT-based Library Checker problems
the summed worst-test time falls from 47.7 s to 29.3 s (1.63x; polynomial `exp` / `pow` / `composition` /
`sqrt` 2.2-2.6x, `consecutive_terms` and `kth_term` 2.3-2.6x), a convolution of two 2^19-term
sequences is **5x faster than ac-library-rs**. An experimental AVX-512 engine (cargo feature
`avx512`, off by default) adds another 1.05x overall, up to 1.4x on large transforms.

## 1. Where the time goes (callgrind, heaviest generated test, instruction counts)

| solution | NTT-related share of all instructions (scalar) | total instructions: scalar -> AVX2 |
|---|---:|---:|
| `exp_of_formal_power_series` | 83% (forward 52%, inverse + pointwise + copies inlined in `mul_trunc` 29%) | 6726M -> 1853M |
| `inv_of_formal_power_series` | 77% | 2476M -> 848M |
| `pow_of_formal_power_series` | 84% | 9492M -> 2496M |
| `multipoint_evaluation` | 72% (+ malloc / free ~14% after the speed-up) | 5448M -> 1871M |
| `bell_number` | 81% | 1543M -> 410M |
| `convolution_mod` | 40%; **the rest is I/O** (number formatting 20%, parsing 29%) | 1216M -> 804M |
| `convolution_mod_large` | ~25% of the 6.3 s; parsing + printing 2 x 16M numbers is ~4 s | |

After the speed-up the transforms still take ~45-55% of the instructions of the polynomial solutions,
`memset` / `memcpy` / `malloc` ~10-20% (zero-padding, temporary vectors in the divide-and-conquer
code), and for the plain convolution problems the shared input / output layer dominates.

## 2. Engines

All engines compute the same transform (decimation in frequency, bit-reversed-like output, so
`intt(pointwise(ntt(a), ntt(b)))` is the cyclic convolution) and the optimised ones are checked to be
word-for-word identical to the portable radix-4 code (`convolution/tests/simd.rs`, `simd512.rs`: all sizes up
to 2^20, five different primes).

| engine | description |
|---|---|
| `r2 naive` | textbook iterative NTT: bit reversal, `u64 %` reductions |
| `r2 table` | radix-2 DIF with per-level twiddle tables and Montgomery products |
| `r4 scalar` | the radix-4 transform of the crate before this change (ac-library-rs style, twiddles from `O(log n)` constants) |
| `r4 avx2` | the same arithmetic on 8 lanes: layers of span >= 16 are vectorised along the span with broadcast twiddles (precomputed `w * (-P^-1)` so the Montgomery quotient comes straight from the input), cache-blocked; the last four layers run inside two registers per 16 elements with lane shuffles and a per-thread table of block rotations |
| `r4 avx512` | sixteen lanes, last five layers in two registers per 32 elements via `vpermt2d`, shuffle tables are compile-time constants |
| ac-library-rs | `ac_library::convolution` (end to end only) |

Forward transform, ms (minimum of 7-25 runs, input copy subtracted; machine: Xeon 2.8 GHz, AVX-512,
4 MB L2 per core) and the inverse transform for the radix-4 engines:

| size | r2 naive | r2 table | r4 scalar | r4 AVX2 | r4 AVX-512 | inv scalar | inv AVX2 | inv AVX-512 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 2^10 | 0.026 | 0.007 | 0.012 | 0.003 | 0.003 | 0.016 | 0.004 | 0.003 |
| 2^11 | 0.066 | 0.017 | 0.025 | 0.007 | 0.005 | 0.024 | 0.006 | 0.005 |
| 2^12 | 0.143 | 0.036 | 0.053 | 0.015 | 0.010 | 0.050 | 0.014 | 0.010 |
| 2^13 | 0.290 | 0.076 | 0.107 | 0.032 | 0.024 | 0.103 | 0.032 | 0.021 |
| 2^14 | 0.576 | 0.163 | 0.227 | 0.068 | 0.046 | 0.218 | 0.071 | 0.046 |
| 2^15 | 1.284 | 0.342 | 0.464 | 0.147 | 0.095 | 0.448 | 0.151 | 0.096 |
| 2^16 | 2.695 | 0.688 | 1.001 | 0.312 | 0.306 | 1.485 | 0.422 | 0.309 |
| 2^17 | 7.318 | 1.593 | 2.065 | 0.673 | 0.434 | 1.983 | 0.694 | 0.426 |
| 2^18 | 12.601 | 3.459 | 4.499 | 1.482 | 1.021 | 4.378 | 1.554 | 1.044 |
| 2^19 | 34.887 | 6.832 | 9.123 | 3.224 | 2.226 | 9.172 | 3.265 | 2.210 |
| 2^20 | 57.091 | 26.070 | 32.227 | 8.288 | 5.264 | 19.274 | 7.174 | 5.192 |
| 2^21 | - | 35.753 | 41.391 | 15.077 | 11.280 | 40.737 | 15.936 | 15.128 |
| 2^22 | - | 91.325 | 86.967 | 32.905 | 30.059 | 83.717 | 30.572 | 31.690 |

Convolution of two sequences of length n/2 (result length n), mod 998244353, ms:

| n | ac-library-rs | ours, portable | ours, best SIMD (AVX-512 here) | vs ac |
|---|---:|---:|---:|---:|
| 2^12 | 0.240 | 0.165 | 0.035 | 6.9x |
| 2^13 | 0.514 | 0.334 | 0.084 | 6.1x |
| 2^14 | 1.118 | 0.730 | 0.207 | 5.4x |
| 2^15 | 2.444 | 1.484 | 0.330 | 7.4x |
| 2^16 | 5.243 | 3.170 | 0.700 | 7.5x |
| 2^17 | 11.221 | 6.726 | 1.508 | 7.4x |
| 2^18 | 23.765 | 14.061 | 3.574 | 6.6x |
| 2^19 | 50.144 | 29.064 | 9.888 | 5.1x |
| 2^20 | 113.744 | 73.809 | 22.012 | 5.2x |
| 2^21 | 263.570 | 137.469 | 53.483 | 4.9x |
| 2^22 | 515.622 | 312.081 | 136.985 | 3.8x |

* The scalar radix-4 code is **not** faster than a plain radix-2 with twiddle tables (`r2 table`);
  the portable fallback could use that, but every judge CPU has AVX2.
* AVX2 is 2.7-4.4x faster than the scalar radix-4 at every size (4x for small, 2.7x above L2), AVX-512
  another 1.3-1.6x for sizes 2^11..2^21. The compute-bound inner loop is the Montgomery product
  (six `vpmuludq` per eight lanes); cache blocking made no measurable difference on this 4 MB-L2
  machine (block sizes 2^11..2^30 all within noise) but is kept for CPUs with a smaller L2.
* AVX-512 needs Rust >= 1.89 for the intrinsics and an `avx512f` CPU, hence the cargo feature
  `avx512` (off by default, not inlined into the submission files). Its first version rebuilt the
  shuffle tables on every call and lost on the many small transforms of divide-and-conquer code
  (0.5-0.9x on `polynomial_root_finding`, `multipoint_evaluation`); with compile-time tables and a
  dispatch threshold of 2^12 it is never slower than AVX2 beyond noise.

## 3. End to end on the Library Checker tests

Worst time over the four largest generated inputs of each problem, ms, same binary with the SIMD
dispatch switched off / AVX2 only / AVX-512 (`NOSIMD`, `NOAVX512` environment switches of the
benchmark harness; outputs identical in all modes). Input parsing and output printing included.

| problem | portable | AVX2 | AVX-512 | AVX2 / portable |
|---|---:|---:|---:|---:|
| `addition_of_big_integers` | 68 | 66 | 65 | 1.03x |
| `addition_of_hex_big_integers` | 68 | 69 | 75 | 0.99x |
| `bell_number` | 849 | 392 | 313 | 2.17x |
| `bernoulli_number` | 310 | 166 | 167 | 1.87x |
| `binomial_coefficient` | 144 | 138 | 162 | 1.04x |
| `binomial_coefficient_prime_mod` | 473 | 459 | 478 | 1.03x |
| `characteristic_polynomial` | 315 | 305 | 297 | 1.03x |
| `composition_of_formal_power_series` | 59 | 25 | 23 | 2.37x |
| `composition_of_formal_power_series_large` | 1579 | 601 | 504 | 2.63x |
| `compositional_inverse_of_formal_power_series` | 128 | 56 | 39 | 2.30x |
| `compositional_inverse_of_formal_power_series_large` | 3017 | 1246 | 997 | 2.42x |
| `consecutive_terms_of_linear_recurrent_sequence` | 2682 | 1018 | 820 | 2.64x |
| `conversion_from_monomial_basis_to_newton_basis` | 789 | 431 | 393 | 1.83x |
| `convolution_mod` | 214 | 166 | 212 | 1.29x |
| `convolution_mod_1000000007` | 392 | 273 | 229 | 1.44x |
| `convolution_mod_2_64` | 990 | 645 | 556 | 1.54x |
| `convolution_mod_large` | 8273 | 6314 | 6340 | 1.31x |
| `division_of_big_integers` | 151 | 153 | 147 | 0.99x |
| `division_of_hex_big_integers` | 131 | 134 | 130 | 0.97x |
| `division_of_polynomials` | 202 | 129 | 132 | 1.57x |
| `exp_of_formal_power_series` | 869 | 395 | 363 | 2.20x |
| `exp_of_formal_power_series_sparse` | 94 | 86 | 89 | 1.10x |
| `factorial` | 238 | 125 | 119 | 1.90x |
| `find_linear_recurrence` | 98 | 102 | 95 | 0.96x |
| `frequency_table_of_tree_distance` | 2282 | 1830 | 1824 | 1.25x |
| `inv_of_formal_power_series` | 351 | 192 | 178 | 1.83x |
| `inv_of_formal_power_series_2d` | 848 | 385 | 419 | 2.20x |
| `inv_of_formal_power_series_sparse` | 72 | 66 | 69 | 1.09x |
| `inv_of_polynomials` | 891 | 435 | 425 | 2.05x |
| `kth_term_of_linearly_recurrent_sequence` | 990 | 424 | 334 | 2.34x |
| `log_of_formal_power_series` | 423 | 223 | 229 | 1.90x |
| `log_of_formal_power_series_sparse` | 120 | 94 | 92 | 1.28x |
| `many_factorials` | 359 | 254 | 254 | 1.41x |
| `montmort_number_mod` | 98 | 94 | 86 | 1.04x |
| `mul_mod2n_convolution` | 818 | 515 | 489 | 1.59x |
| `mul_modp_convolution` | 204 | 169 | 174 | 1.20x |
| `multiplication_of_big_integers` | 109 | 82 | 90 | 1.33x |
| `multiplication_of_hex_big_integers` | 91 | 78 | 87 | 1.16x |
| `multipoint_evaluation` | 750 | 402 | 342 | 1.87x |
| `multipoint_evaluation_on_geometric_sequence` | 232 | 166 | 169 | 1.40x |
| `multivariate_convolution` | 1052 | 692 | 700 | 1.52x |
| `multivariate_convolution_cyclic` | 935 | 613 | 601 | 1.52x |
| `number_of_increasing_sequences_between_two_sequences` | 1160 | 472 | 462 | 2.46x |
| `number_of_subsequences` | 111 | 97 | 93 | 1.15x |
| `partition_function` | 332 | 151 | 160 | 2.19x |
| `polynomial_interpolation` | 926 | 483 | 474 | 1.92x |
| `polynomial_interpolation_on_geometric_sequence` | 286 | 214 | 197 | 1.34x |
| `polynomial_root_finding` | 807 | 379 | 378 | 2.13x |
| `polynomial_taylor_shift` | 149 | 119 | 115 | 1.25x |
| `pow_of_formal_power_series` | 1193 | 501 | 436 | 2.38x |
| `pow_of_formal_power_series_sparse` | 103 | 108 | 105 | 0.96x |
| `prefix_sum_of_polynomial` | 798 | 396 | 325 | 2.02x |
| `product_of_polynomial_sequence` | 1231 | 681 | 582 | 1.81x |
| `q_binomial_coefficient_prime_mod` | 742 | 747 | 760 | 0.99x |
| `sharp_p_subset_sum` | 1016 | 494 | 434 | 2.06x |
| `shift_of_sampling_points_of_polynomial` | 269 | 199 | 207 | 1.35x |
| `sqrt_of_formal_power_series` | 800 | 354 | 290 | 2.26x |
| `sqrt_of_formal_power_series_sparse` | 138 | 114 | 104 | 1.21x |
| `stirling_number_of_the_first_kind` | 203 | 111 | 111 | 1.84x |
| `stirling_number_of_the_first_kind_fixed_k` | 614 | 233 | 227 | 2.64x |
| `stirling_number_of_the_first_kind_small_p_large_n` | 469 | 397 | 500 | 1.18x |
| `stirling_number_of_the_second_kind` | 135 | 123 | 126 | 1.09x |
| `stirling_number_of_the_second_kind_fixed_k` | 577 | 242 | 206 | 2.39x |
| `stirling_number_of_the_second_kind_small_p_large_n` | 126 | 138 | 126 | 0.91x |
| `sum_of_exponential_times_polynomial` | 1060 | 1097 | 1092 | 0.97x |
| `sum_of_exponential_times_polynomial_limit` | 827 | 819 | 817 | 1.01x |
| `wildcard_pattern_matching` | 847 | 377 | 336 | 2.24x |
| **sum** | **47680** | **29253** | **27972** | **1.63x** |

Problems near 1.0x are not transform-bound (sparse series, big-integer addition / division, the
Stirling `small_p_large_n` and `q_binomial` families, `sum_of_exponential_*`).

## 4. Other changes made on the way

* `convolution_mod` folds the pointwise product and the `1 / n` of the inverse transform into one
  vector pass (`pointwise_scaled`) and pads its inputs with a single write.
* `polynomial::fps::mul` passes the residues to the convolution without copying (`StaticModInt` is now
  `repr(transparent)`) and takes the result vector over: ~5-15% on the polynomial solutions.
* `ntt` / `intt_unscaled` / `intt` dispatch automatically; `ntt_scalar` / `intt_unscaled_scalar` keep the portable
  code, `simd::set_enabled(false)` switches the dispatch off for A/B runs.

## 5. Not done / observations

* **Shared input / output.** With the template reader (`str::parse` per token) and `write!`-based
  writer, `convolution_mod` spends about half of its time in I/O; a hand-written integer parser and
  `itoa`-style printer took the heaviest test from 138 ms to 103 ms (`convolution_mod`), and
  `convolution_mod_large` (6.3 s, 33M numbers in, 33M out) is mostly I/O. This would help every solution but
  changes the template, so it was left as a proposal.
* **Floating-point FFT** (complex doubles, 15-bit splitting, 4 transforms of size n instead of 9 NTTs
  for an arbitrary modulus) was not implemented: the AVX2 NTT already makes the 3-prime convolution
  (`convolution_arbitrary_mod`, `convolution_mod_1000000007`: 392 -> 273 ms) competitive and
  exact, and a good vectorised complex FFT is a large piece of code for one problem.
* **Lazy reductions** (values kept in `[0, 2P)`, valid for `P < 2^30`) would save two instructions of
  thirteen per Montgomery product (~10%); not implemented.
* The last four layers cost ~30% of a transform for 20% of the layers (shuffles and vector twiddles);
  a radix-4 in-register tail is the next candidate.
