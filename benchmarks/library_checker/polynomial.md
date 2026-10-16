# Library Checker `polynomial`: dsa-templates vs ac-library-rs

All 28 problems of the [Library Checker](https://judge.yosupo.jp/) `polynomial` category are
solved with the Rust crates of this repository and every solution passes the problem's own
checker on every generated test. This page compares run time with
[ac-library-rs](https://github.com/rust-lang-ja/ac-library-rs) (`ac`) wherever ac has something
to build on.

## Method

* Same I/O layer for both (`FastScan` / `Out` from `template/rust`), `cargo build --release`.
* For each problem every test is run once per implementation and checked with the problem's
  checker; the heaviest test is then re-run 5 times, alternating implementations. The table
  reports the **minimum** wall-clock time in milliseconds, input parsing and output printing
  included.
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10%.
* ac-library-rs offers `convolution` and modular integers but **no polynomial or power-series
  routines**. The "ac" column is therefore **—** for most rows. For the standard Newton-iteration
  problems (inverse, logarithm, exponential, power, division, Taylor shift) the ac column is the
  *same algorithm hand-written on top of ac's `convolution`* (`lcbench::acpoly`), which isolates the
  cost of the convolution layer: the NTT of this repository is radix-4 with Montgomery
  reductions.

## Summary

Compared against ac on 6 problems: geometric mean of ours/ac = **0.72**; 6 within 10% or faster, 0 slower than that.

New crates driven by this category: `polynomial` — `fps` (inverse, log, exp, pow, sqrt,
division with remainder, bivariate inverse), `sparse` (the `O(N K)` versions), `special` (Taylor
shift, shifting sampling points, geometric evaluation and interpolation, prefix sums, products of
sequences), `multipoint` (subproduct tree: evaluation, interpolation, Newton basis), `gcd`
(half-gcd, extended gcd, inverse modulo a polynomial), `roots` (root finding over `F_p`) and
`composition` (Kinoshita–Li composition and compositional inverse in `O(n log^2 n)`) — and
`poly_mod_p` (factorisation over `F_p` for arbitrary primes). Each has a randomized stress test
against the definition or a naive computation.

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `composition_of_formal_power_series` | 10 | `modint`, `polynomial` | 57 | — | — | `compose`: Kinoshita-Li, the transposed power projection (Bostan-Mori on 1 - y g(x), 2D products by Kronecker substitution). ac: nothing. |
| `composition_of_formal_power_series_large` | 10 | `modint`, `polynomial` | 1393 | — | — | Same `compose` at N = 131072: O(n log^2 n), about 17 rounds of bivariate products of total size O(n). ac: nothing. |
| `compositional_inverse_of_formal_power_series` | 10 | `modint`, `polynomial` | 109 | — | — | `compositional_inverse`: Newton's iteration with one composition per doubling (f'(g) = (f(g))' / g'). ac: nothing. |
| `compositional_inverse_of_formal_power_series_large` | 10 | `modint`, `polynomial` | 2551 | — | — | Same at N = 131072; about two compositions' worth of work. ac: nothing. |
| `conversion_from_monomial_basis_to_newton_basis` | 10 | `modint`, `polynomial` | 673 | — | — | `monomial_to_newton`: divide by the product of the first half of the points on the subproduct tree and recurse into both halves. ac: nothing. |
| `division_of_polynomials` | 10 | `modint`, `polynomial` | 331 | 448 | 0.74× | Quotient from the reversed polynomials and a Newton inverse, remainder by one product. ac column: the same algorithm hand-written with ac's `convolution`. |
| `exp_of_formal_power_series` | 10 | `modint`, `polynomial` | 766 | 1158 | 0.66× | Newton iteration `g <- g (1 - log g + f)`. ac column: the same, with ac's `convolution`. |
| `exp_of_formal_power_series_sparse` | 10 | `modint`, `polynomial` | 75 | — | — | O(N K) recurrence from g' = f' g. ac: nothing. |
| `factorization_of_polynomials` | 10 | `poly_mod_p` | 1139 | — | — | `poly_mod_p`: Yun's square-free split (p-th root when the derivative vanishes), distinct-degree factorisation, Cantor-Zassenhaus (trace map for p = 2). Schoolbook arithmetic, degree <= 100. ac: nothing. |
| `inv_of_formal_power_series` | 10 | `modint`, `polynomial` | 298 | 372 | 0.80× | Newton iteration. ac column: the same, with ac's `convolution` on `ModInt998244353`. |
| `inv_of_formal_power_series_2d` | 10 | `modint`, `polynomial` | 864 | — | — | Newton iteration in x; the bivariate products are Kronecker-substituted into one convolution. ac: nothing. |
| `inv_of_formal_power_series_sparse` | 10 | `modint`, `polynomial` | 58 | — | — | O(N K) recurrence. ac: nothing. |
| `inv_of_polynomials` | 10 | `modint`, `polynomial` | 871 | — | — | Extended Euclid with the recursive half-gcd (O(n log^2 n)); the cofactor comes from the accumulated 2x2 matrix. ac: nothing. |
| `log_of_formal_power_series` | 10 | `modint`, `polynomial` | 354 | 526 | 0.67× | f' * inv(f), integrated. ac column: the same, with ac's `convolution`. |
| `log_of_formal_power_series_sparse` | 10 | `modint`, `polynomial` | 73 | — | — | O(N K) recurrence from f = f g'. ac: nothing. |
| `multipoint_evaluation` | 10 | `modint`, `polynomial` | 630 | — | — | Subproduct tree and the transposed (middle product) descent, O(n log^2 n). ac: nothing. |
| `multipoint_evaluation_on_geometric_sequence` | 10 | `modint`, `polynomial` | 205 | — | — | Bluestein/chirp-z: one correlation. ac: nothing. |
| `polynomial_interpolation` | 10 | `modint`, `polynomial` | 837 | — | — | Subproduct tree, weights from the multipoint evaluation of P', bottom-up combination. ac: nothing. |
| `polynomial_interpolation_on_geometric_sequence` | 5 | `modint`, `polynomial` | 187 | — | — | Newton form on a geometric progression plus the q-binomial theorem: two convolutions, O(n log n). ac: nothing. |
| `polynomial_root_finding` | 10 | `modint`, `polynomial` | 834 | — | — | x^p mod f by repeated squaring (NTT products), gcd with x^p - x, then random splitting with (x + a)^((p-1)/2) - 1. ac: nothing. |
| `polynomial_taylor_shift` | 5 | `modint`, `polynomial` | 119 | 151 | 0.79× | One correlation with factorials. ac column: the same, with ac's `convolution`. |
| `pow_of_formal_power_series` | 10 | `modint`, `polynomial` | 1098 | 1696 | 0.65× | Leading zeros shifted out, exp(e log f) with e reduced mod p, constant^(e mod p-1). ac column: the same with the ac-based log/exp above. |
| `pow_of_formal_power_series_sparse` | 10 | `modint`, `polynomial` | 92 | — | — | O(N K) recurrence from h g' = e h' g, after removing the leading monomial. ac: nothing. |
| `prefix_sum_of_polynomial` | 5 | `modint`, `polynomial` | 674 | — | — | Faulhaber with Bernoulli numbers (inverse of (e^x - 1)/x): one inverse and one correlation, O(n log n). ac: nothing. |
| `product_of_polynomial_sequence` | 10 | `modint`, `polynomial` | 1041 | — | — | Merge the two shortest first with a heap: O(D log D log N). ac: nothing. |
| `shift_of_sampling_points_of_polynomial` | 5 | `modint`, `polynomial` | 198 | — | — | Lagrange weights and one convolution with the inverses of consecutive integers; arguments that land on the samples are copied, runs are split at wrap-arounds. ac: nothing. |
| `sqrt_of_formal_power_series` | 10 | `modint`, `polynomial` | 628 | — | — | Tonelli-Shanks for the leading term, Newton g <- (g + f/g)/2; an odd leading exponent or a non-residue gives -1. ac: nothing. |
| `sqrt_of_formal_power_series_sparse` | 10 | `modint`, `polynomial` | 107 | — | — | O(N K) from the pow recurrence with exponent 1/2. ac: nothing. |

