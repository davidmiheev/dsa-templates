# Library Checker `number_theory`: dsa-templates vs ac-library-rs

All 26 problems of the [Library Checker](https://judge.yosupo.jp/) `number_theory` category are
solved with the Rust crates of this repository and every solution passes the problem's own
checker on every generated test. This page compares run time with
[ac-library-rs](https://github.com/rust-lang-ja/ac-library-rs) (`ac`) wherever ac has something
to build on, and notes the run time of the reference C++ solution of the Library Checker
repository on the heaviest test where that is informative.

## Method

* Same I/O layer for both (`FastScan` / `Out` from `template/rust`), `cargo build --release`.
* For each problem every test is run once per implementation and checked with the problem's
  checker; the heaviest test is then re-run 3 times (the three largest problems once),
  alternating implementations. The table reports the **minimum** wall-clock time in
  milliseconds, input parsing and output printing included.
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10%.
* ac-library-rs only has `floor_sum`, `crt`, `inv_mod`, `pow_mod`, an `i32` `primitive_root` and
  `convolution`. The "ac" column is therefore **—** for most rows; where it is filled the ac side
  is `floor_sum` itself, or (Bernoulli numbers) the same Newton-iteration algorithm hand-written
  on top of ac's `convolution`.
* "Reference C++" is the `sol/correct.cpp` shipped with the problem, compiled by the Library
  Checker scripts and timed on the same test and machine.

## Summary

Compared against ac on 2 problems: geometric mean of ours/ac = **0.82**; 2 within 10% or faster, 0 slower than that.

New crates driven by this category:

* `number_theory` (extended): Montgomery-form Miller–Rabin and Pollard rho (3–4x faster than the
  previous `u128 %` version), `isqrt`, `kth_root`, `kth_root_mod`, `discrete_log`,
  `discrete_log_prime` / `PrimeLog` (Pohlig–Hellman, batched baby-step giant-step with prefetch),
  `floor_sum`, `min_mod_linear`, `euler_phi`, `tetration_mod`, `NimMultiplier`, `gaussian_gcd`
  and `two_squares`.
* `stern_brocot`: paths, `lca`, `ancestor`, subtree ranges and best rational approximations.
* `prime_sieve`: segmented prime sieve, linear Moebius / totient sieves.
* `sublinear`: sums over floor quotients — Lucy prime counting, Min_25's sieve, Dirichlet products
  and inverses of prefix sums, the `n^{2/3}` division trick (`totient_sum`, `mertens`) and
  `squarefree_count`.
* `polynomial::special::bernoulli_numbers`.

Every new function has a randomized stress test against the definition or a brute-force
computation.

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `bernoulli_number` | 10 | `modint`, `polynomial` | 285 | 431 | 0.66× | `polynomial::special::bernoulli_numbers`: `B_j / j!` are the coefficients of `x / (e^x - 1)`, one Newton inverse. ac column: the same with ac's `convolution` (`lcbench::acpoly`). |
| `counting_primes` | 5 | `sublinear` | 253 | — | — | `sublinear::prime_pi`: Lucy_Hedgehog, `O(N^{3/4})`; unchecked indexing and float reciprocals instead of 64-bit divisions. ac: nothing. Reference C++: 485 ms. |
| `counting_squarefrees` | 10 | `sublinear` | 1033 | — | — | `sublinear::squarefree_count`: `sum_{d <= L} mu(d) floor(N / d^2)` from a linear sieve with `L = (4e5 N)^{1/3}`, the tail regrouped as `sum_k (M(sqrt(N/k)) - M(L))` with Mertens values by the division trick. ac: nothing. Reference C++: 1645 ms. |
| `dirichlet_convolution_and_prefix_sums` | 30 | `sublinear` | 10673 | — | — | `sublinear::dirichlet_product`: hyperbola method for every quotient, `O(N^{3/4})`; F and G interleaved so one cache line serves both lookups, `floor(x / i)` as the large entry `k i` or a float reciprocal, lazy `u64` accumulation. ac: nothing. Reference C++: 10.8 s (an `O(N^{2/3} log N)` bucket method). |
| `dirichlet_inverse_and_prefix_sums` | 30 | `sublinear` | 11089 | — | — | `sublinear::dirichlet_inverse`: the same hyperbola identity solved for `G(x)`, quotients in increasing order. ac: nothing. Reference C++: 4.6 s (the `N^{2/3}` method); ours is 2.4x slower but a third of the limit. |
| `discrete_logarithm_fixed_mod` | 10 | `number_theory` | 2300 | — | — | `number_theory::PrimeLog::log_many`: `2^24`-entry baby-step table (key and value in one cache line, filled with prefetching), giant steps for a window of 32 queries at a time with software prefetch. ac: nothing. Reference C++: 641 ms (Farey-sequence technique). |
| `discrete_logarithm_mod` | 10 | `number_theory` | 65 | — | — | `number_theory::discrete_log`: extended BSGS (peels the common factors of `x` and `m`), Barrett reduction, open-addressing table. ac: nothing. |
| `enumerate_primes` | 10 | `prime_sieve` | 905 | — | — | `prime_sieve::primes_up_to`: odd-only segmented sieve in 128 KiB segments. ac: nothing. Reference C++: 2.3 s. |
| `enumerate_quotients` | 5 | `number_theory` | 75 | — | — | `isqrt`; the set is `{1..sqrt(N)}` plus `floor(N / i)`. ac: nothing. |
| `factorize` | 10 | `number_theory` | 127 | — | — | `number_theory::factorize`: trial division by 15 primes, Miller-Rabin and Brent's rho, both in 64-bit Montgomery arithmetic. ac: nothing. Reference C++ (plain rho, one gcd per step): 5x faster on one hand-made test (23 vs 134 ms), 10-15x slower on the big semiprime tests (0.5-1.2 s vs 35-70 ms). |
| `gcd_of_gaussian_integers` | 5 | `number_theory` | 116 | — | — | `number_theory::gaussian_gcd`: Euclid with rounded Gaussian division in `i128`. ac: nothing. |
| `kth_root_integer` | 10 | `number_theory` | 96 | — | — | `number_theory::kth_root`: float estimate corrected with exact `u128` powers; `isqrt` for `k = 2`. ac: nothing. |
| `kth_root_mod` | 5 | `number_theory` | 37 | — | — | `number_theory::kth_root_mod`: the part of `p - 1` coprime to `k` is inverted by one exponentiation; Pohlig-Hellman discrete logarithms only for the primes shared with `k`. ac: nothing. Reference C++: 72 ms. |
| `min_of_mod_of_linear` | 10 | `number_theory` | 45 | — | — | `number_theory::min_mod_linear`: Euclid-like recursion in `O(log m)`: the minima of the blocks between wraps form a linear sequence modulo a smaller modulus, mirrored when `a > m / 2`. ac: nothing. |
| `nim_product_64` | 10 | `number_theory` | 433 | — | — | `number_theory::NimMultiplier`: three recursive products per halving (Karatsuba-like) and a 64 KiB table for the 8-bit base case. ac: nothing. |
| `primality_test` | 5 | `number_theory` | 206 | — | — | `number_theory::is_prime`: deterministic Miller-Rabin, 7 bases, 64-bit Montgomery arithmetic. ac: nothing. Reference C++ (`__int128 %`, 12 bases): 618 ms. |
| `primitive_root` | 10 | `number_theory` | 24 | — | — | `number_theory::primitive_root`: Pollard-rho factorisation of `p - 1`, then candidates in increasing order. ac: `primitive_root` takes an `i32` modulus, so the `10^18` range is out of reach. |
| `rational_approximation` | 5 | `stern_brocot` | 160 | — | — | `stern_brocot::best_approximations`: Stern-Brocot descent with run jumps, `O(log N)`. ac: nothing. |
| `sqrt_mod` | 10 | `number_theory` | 121 | — | — | `number_theory::mod_sqrt` (Tonelli-Shanks, `NONE` for non-residues). ac: nothing. |
| `stern_brocot_tree` | 5 | `stern_brocot` | 160 | — | — | `stern_brocot`: run-length paths from the continued fraction (`encode_path`, `decode_path`), `lca` by comparing runs, `ancestor` by truncating the path, `descendant_range` from the Farey neighbours. ac: nothing. |
| `sum_of_floor_of_linear` | 5 | `number_theory` | 37 | 36 | 1.03× | `number_theory::floor_sum` (exact `u128` arithmetic). ac column: `ac_library::floor_sum`. |
| `sum_of_multiplicative_function` | 10 | `prime_sieve`, `sublinear` | 639 | — | — | `sublinear`: Lucy for `pi` and the prime sum modulo 469762049, then Min_25's DFS with `f(p) = a + b p`. ac: nothing. Reference C++: 641 ms. |
| `sum_of_multiplicative_function_large` | 30 | `prime_sieve`, `sublinear` | 20234 | — | — | Same code at `N = 10^13`: Lucy about 11 s, DFS about 8 s (leaf primes `p^2 <= v < p^3` summed in a tight loop with prefix sums of `f`). ac: nothing. Reference C++: 24.0 s. |
| `sum_of_totient_function` | 5 | `sublinear` | 277 | — | — | `sublinear::totient_sum`: `phi * 1 = id`, so `Phi(v) = v(v+1)/2 - sum_{d >= 2} Phi(v / d)` with a linear-sieve table up to `2 N^{2/3}`, `O(N^{2/3})`. ac: nothing. Reference C++: 418 ms. |
| `tetration_mod` | 10 | `number_theory` | 12 | — | — | `number_theory::tetration_mod`: generalised Euler theorem with capped values along the chain `m, phi(m), phi(phi(m)), ...` (at most 30 levels). ac: nothing. |
| `two_square_sum` | 5 | `number_theory` | 57 | — | — | `number_theory::two_squares`: factorise, a Gaussian prime above each `p = 1 (mod 4)` from `sqrt(-1)` and a Gaussian gcd, all splits `pi^j conj(pi)^(e-j)`, the four unit rotations, sort and dedupe. ac: nothing. |

