# Library Checker `convolution`: dsa-templates vs ac-library-rs

All 16 problems of the [Library Checker](https://judge.yosupo.jp/) `convolution` category are
solved with the Rust crates of this repository and every solution passes the problem's own
checker on every generated test. This page compares run time with
[ac-library-rs](https://github.com/rust-lang-ja/ac-library-rs) (`ac`) wherever ac has something
to build on.

## Method

* Same I/O layer for both (`FastScan` / `Out` from `template/rust`), `cargo build --release`
  (no `target-cpu` flags; the only architecture-specific code is the run-time detected
  `PCLMULQDQ` path of `gf2_64`).
* For each problem every test is run once per implementation and checked with the problem's
  checker; the heaviest test is then re-run 5 times (3 for `convolution_mod_large`, whose inputs
  are 335 MB). The table reports the **minimum** wall-clock time in milliseconds, input parsing
  and output printing included.
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10%.
* ac-library-rs provides number-theoretic transforms for static moduli (`convolution`,
  `convolution_raw`, `convolution_i64`). Everything else here — bitwise and divisor convolutions,
  min-plus, multiplicative index structure, multivariate products, characteristic 2 — has no ac
  counterpart, so the "ac" column is **—** for those rows. Where a layer on top of ac's NTT is
  possible (arbitrary modulus), the ac column is that hand-written layer.

## Summary

Compared against ac on 2 problems: geometric mean of ours/ac = **0.97**; 2 within 10% or faster, 0 slower than that.

> **Note on the timings.** They were measured with the portable radix-4 transforms. The `convolution` crate now runs its NTT on AVX2 when the CPU has it (identical results): the NTT-bound solutions of this table are 1.4-2.6x faster, see [`../ntt_engines.md`](../ntt_engines.md) for the profile, the engines compared and per-problem before / after numbers.

New or reworked crates driven by this category: `convolution` (radix-4 Montgomery NTT without bit
reversal, arbitrary modulus by three primes, block decomposition for results beyond `2^23`,
convolution modulo `2^64`), `set_power_series` (zeta / Moebius / Walsh-Hadamard, AND / OR / XOR
and subset convolution), `divisor_transform` (GCD / LCM convolution), `min_plus` (convex and
concave min-plus convolution, monotone minima), `multiplicative_convolution` (products modulo a
prime and modulo `2^n`), `multivariate_convolution` (truncated and cyclic, Bluestein DFT),
`gf2_64` (field arithmetic and additive FFT) and `number_theory` (Miller-Rabin, Pollard rho,
primitive roots). Each has a randomized stress test against the definition.

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `bitwise_and_convolution` | 5 | `set_power_series` | 136 | — | — | Superset zeta transform, pointwise product, Moebius (`and_convolution`). ac: nothing. |
| `bitwise_xor_convolution` | 5 | `set_power_series` | 132 | — | — | Walsh-Hadamard transform (`xor_convolution`). ac: nothing. |
| `convolution_F_2_64` | 10 | `gf2_64` | 687 | — | — | Polynomial product over GF(2^64): additive FFT over a Cantor basis (`x^2 + x` halving, XOR-only Taylor split), carry-less multiplication with PCLMULQDQ. O(n log n) field multiplications. ac: nothing. |
| `convolution_mod` | 5 | `convolution` | 170 | 166 | 1.02× | Radix-4 Montgomery NTT. ac column: `convolution` on `ModInt998244353`. Both are dominated by parsing and printing 2^20 numbers. |
| `convolution_mod_1000000007` | 10 | `convolution` | 341 | 373 | 0.91× | `convolution_arbitrary_mod`: three NTT primes and Garner. ac column: `convolution_raw` on three hand-defined moduli and Garner written by hand. |
| `convolution_mod_2_64` | 10 | `convolution` | 880 | — | — | `convolution_mod_2_64`: 16-bit limbs, the four limb-position sums formed in the transform domain, two primes (24 transforms). ac has no such routine. |
| `convolution_mod_large` | 10 | `convolution` | 7055 | — | — | `convolution_large`: 2^22-blocks transformed once at size 2^23 (the largest the prime allows), spectra of block pairs accumulated per output block, overlap-add. 15 transforms of 2^23 instead of 9 of 2^25 over three primes. Input and output are ~340 MB, so parsing and printing are a large share. ac: a single transform is limited to 2^23 and nothing blocks it. |
| `gcd_convolution` | 5 | `divisor_transform` | 142 | — | — | Zeta transform over multiples (one pass per prime), pointwise product, Moebius. O(n log log n). ac: nothing. |
| `lcm_convolution` | 5 | `divisor_transform` | 137 | — | — | Zeta transform over divisors, pointwise product, Moebius. ac: nothing. |
| `min_plus_convolution_concave_arbitrary` | 5 | `min_plus` | 452 | — | — | b cut into blocks of the length of a; each block's triangles are halved recursively into fully valid rectangles solved by monotone minima (inverse Monge). O((n+m) log n). ac: nothing. |
| `min_plus_convolution_convex_arbitrary` | 5 | `min_plus` | 93 | — | — | Monotone minima (divide and conquer) on the band matrix, which is Monge for convex a. ac: nothing. |
| `min_plus_convolution_convex_convex` | 5 | `min_plus` | 75 | — | — | Merge of the two slope sequences, O(n+m). ac: nothing. |
| `mul_mod2n_convolution` | 5 | `multiplicative_convolution` | 645 | — | — | Odd parts as +-5^t (Z/2 x Z/2^(M-2)), folded per pair of 2-adic valuations; each pair is a 2-by-2^(M-2) cyclic convolution done with four forward and two inverse NTTs. ac: nothing. |
| `mul_modp_convolution` | 5 | `multiplicative_convolution` | 140 | — | — | Discrete logarithms via a primitive root turn products into a cyclic convolution of length p-1; index 0 from sums. ac: nothing. |
| `multivariate_convolution` | 10 | `multivariate_convolution` | 1019 | — | — | Mixed-radix Kronecker substitution with a carry-count grading chi mod K: K transforms per input, K^2 pointwise products. ac: nothing. |
| `multivariate_convolution_cyclic` | 10 | `multivariate_convolution` | 903 | — | — | Per-axis DFT over F_p (naive up to length 64, Bluestein with an arbitrary-modulus convolution beyond), pointwise product, inverse. ac: nothing. |

