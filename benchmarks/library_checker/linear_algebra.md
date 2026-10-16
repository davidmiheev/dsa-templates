# Library Checker `linear_algebra`: dsa-templates vs ac-library-rs

All 18 problems of the [Library Checker](https://judge.yosupo.jp/) `linear_algebra` category are
solved with the Rust crates of this repository and every solution passes the problem's own
checker on every generated test. This page compares run time with
[ac-library-rs](https://github.com/rust-lang-ja/ac-library-rs) (`ac`) wherever ac has something to
build on, and notes the run time of the reference C++ solution of the Library Checker repository
for a few heavy problems.

## Method

* Same I/O layer for both (`FastScan` / `Out` from `template/rust`), `cargo build --release`.
* For each problem every test is run once per implementation and checked with the problem's
  checker; the heaviest test is then re-run 3 times, alternating implementations. The table
  reports the **minimum** wall-clock time in milliseconds, input parsing and output printing
  included.
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10%.
* ac-library-rs has **no matrix or linear-algebra routines**. For the three basic problems
  (product, determinant, inverse) the ac column is the *straightforward implementation written on
  top of ac's `ModInt998244353`*, i.e. what a user of ac would write; everywhere else it is **—**.

## Summary

Compared against ac on 3 problems: geometric mean of ours/ac = **0.62**; 3 within 10% or faster, 0 slower than that.

New crates driven by this category:

* `modmatrix` (extended): a flat `Matrix` with a tiled AVX2-selected product, power, determinant,
  rank, inverse, adjugate (also for singular matrices), linear systems with kernel bases,
  characteristic polynomial (Hessenberg), Pfaffian; `hafnian`, `sparse_determinant` (Wiedemann) and
  `determinant_mod` (any modulus).
* `bit_matrix`: matrices over `F_2` in 64-bit words — four-Russians product, rank, determinant,
  inverse, linear systems — and `intersect_spaces`.
* `linear_recurrence`: Berlekamp–Massey.

Each has a randomized stress test against the definition (Laplace expansion, Pfaffian / hafnian
expansion, cofactors, brute-force subspaces, ...).

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `adjugate_matrix` | 5 | `modmatrix` | 301 | — | — | `Matrix::adjugate`: `det * inverse` when invertible; zero when the rank is at most n-2; for rank n-1 the rank-one product `c x y^T` of the right and left kernel vectors, `c` from one cofactor. `O(n^3)` in all cases. ac: nothing. Reference C++: 849 ms. |
| `characteristic_polynomial` | 10 | `modmatrix` | 189 | — | — | `Matrix::characteristic_polynomial`: similarity reduction to upper Hessenberg form, then the `O(n^3)` three-term recurrence on the leading minors. ac: nothing. Reference C++: 213 ms. |
| `hafnian_of_matrix` | 10 | `modmatrix` | 1361 | — | — | `modmatrix::hafnian`: recursion over pairs of vertices with polynomial entries (include / exclude the pair, merge it into the rest), `O(2^{n/2} n^2 d^2)`; lazy reduction of the `u64` accumulators, no allocation below the recursion buffers. ac: nothing. Reference C++: 1.58 s. |
| `intersection_of_f2_vector_spaces` | 10 | `bit_matrix` | 619 | — | — | `bit_matrix::intersect_spaces`: Zassenhaus on `(u | u)`, `(v | 0)` packed in one `u64`; rows whose left half vanishes are the intersection. ac: nothing. Reference C++: 712 ms. |
| `inverse_matrix` | 5 | `modmatrix` | 204 | 244 | 0.84× | `Matrix::inverse`: Gauss-Jordan on `[A | I]`; rows updated with Shoup's division-free modular multiplication (vectorizes). Ac column: the same elimination on `ModInt998244353`. Reference C++: 245 ms. |
| `inverse_matrix_mod_2` | 10 | `bit_matrix` | 785 | — | — | `BitMatrix::inverse`: Gauss-Jordan on `[A | I]` with 64-bit word row operations, `O(n^3 / 64)`. ac: nothing. Reference C++: 1.78 s. |
| `matrix_det` | 5 | `modmatrix` | 54 | 59 | 0.90× | `Matrix::det`: elimination with Shoup multiplication. Ac column: the same on `ModInt998244353`. Reference C++: 85 ms. |
| `matrix_det_arbitrary_mod` | 5 | `modmatrix` | 313 | — | — | `determinant_mod`: Euclid on the rows instead of a division (no inverses needed for composite modulus), Barrett reduction. ac: nothing. Reference C++: 283 ms. |
| `matrix_det_mod_2` | 5 | `bit_matrix` | 144 | — | — | `BitMatrix::det` (full rank test by word-parallel elimination). ac: nothing. |
| `matrix_product` | 10 | `modmatrix` | 464 | 1464 | 0.32× | `Matrix::mul`: i-k-j loop tiled by 256 columns, `u64` accumulators reduced every 18 terms, AVX2 build selected at run time. Ac column: the same loop order on `ModInt998244353` (a reduction per multiplication). Reference C++: 1.55 s. |
| `matrix_product_mod_2` | 10 | `bit_matrix` | 174 | — | — | `BitMatrix::mul`: method of four Russians (tables of the 256 xor-combinations of 8 rows). ac: nothing. Reference C++: 1.22 s. |
| `matrix_rank` | 5 | `modmatrix` | 55 | — | — | `Matrix::rank`: forward elimination only (no back substitution). ac: nothing. |
| `matrix_rank_mod_2` | 10 | `bit_matrix` | 213 | — | — | `BitMatrix::rank` after transposing so that the matrix has the fewer rows; word-parallel elimination. ac: nothing. Reference C++: 191 ms. |
| `pfaffian_of_matrix` | 10 | `modmatrix` | 95 | — | — | `Matrix::pfaffian`: congruence elimination (Parlett-Reid style) clearing row `k` with a symmetric row / column operation, `O(n^3)`. ac: nothing. Reference C++: 181 ms. |
| `pow_of_matrix` | 5 | `modmatrix` | 191 | — | — | `Matrix::pow`: square-and-multiply on the tiled product, ~120 products of order 200. ac: nothing. Reference C++: 608 ms. |
| `sparse_matrix_det` | 5 | `modmatrix` | 356 | — | — | `modmatrix::sparse_determinant`: Wiedemann with a random diagonal preconditioner; `2n` sparse matrix-vector products and `linear_recurrence::berlekamp_massey`; `det(A D) = (-1)^{n+1} c_n`. ac: nothing. Reference C++: 172 ms. |
| `system_of_linear_equations` | 5 | `modmatrix` | 76 | — | — | `Matrix::solve`: Gauss-Jordan on `[A | b]`, one particular solution and a kernel basis from the free columns. ac: nothing. |
| `system_of_linear_equations_mod_2` | 5 | `bit_matrix` | 292 | — | — | `BitMatrix::solve`: the same over `F_2` with packed rows. ac: nothing. |

