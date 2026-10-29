# Library Checker `set_power_series`: dsa-templates vs ac-library-rs

All 5 problems of the [Library Checker](https://judge.yosupo.jp/) `set_power_series` category are
solved with the `set_power_series` crate of this repository and every solution passes the
problem's own checker on every generated test.

## Method

* Same I/O layer for both (`FastScan` / `Out` from `template/rust`), `cargo build --release`.
* For each problem every test is run once per implementation and checked with the problem's
  checker; the heaviest test is then re-run 5 times, alternating implementations. The table
  reports the **minimum** wall-clock time in milliseconds, input parsing and output printing
  included.
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10%.
* ac-library-rs has nothing for set power series (it has no subset convolution, exp or log on the
  subset ring), so the ac column is **—** everywhere except `subset_convolution`, where it is a
  straightforward layered implementation on ac's `ModInt998244353`.

## Summary

Compared against ac on 1 problems: geometric mean of ours/ac = **0.56**; 1 within 10% or faster, 0 slower than that.

> **Note on the timings.** They were measured with the portable radix-4 transforms. The `convolution` crate now runs its NTT on AVX2 when the CPU has it (identical results): the NTT-bound solutions of this table are 1.4-2.6x faster, see [`../ntt_engines.md`](../ntt_engines.md) for the profile, the engines compared and per-problem before / after numbers.

New code driven by this category (all in `set_power_series`): `subset_convolution` rewritten on a
`[mask][rank]` layout, `subset_exp`, `subset_log`, `subset_compose` (composition with a univariate
polynomial) and `subset_power_projection` (its transpose). The ranked representation turns the
products, `exp` and `log` into per-mask polynomial operations, and composition / power projection
need only `O(n^2 2^n + M n)`. Each has a randomized stress test against the definition (powers by
repeated subset convolution).

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `exp_of_set_power_series` | 10 | `set_power_series` | 973 | — | — | Ranked transform, then per mask the O(n * popcount) recurrence `k g_k = sum j f_j g_{k-j}` in the rank variable, inverse transform. ac: nothing. |
| `log_of_set_power_series` | 10 | `set_power_series` | 955 | — | — | Same ranked domain, recurrence `g_k = f_k - (1/k) sum_{j<k} j g_j f_{k-j}`. ac: nothing. |
| `polynomial_composite_set_power_series` | 10 | `set_power_series` | 1148 | — | — | Taylor values f^(k)(s_0), k <= n, by synthetic division (O(M n)); then peel one variable at a time with `h_k(a + x b) = h_k(a) + x b h_{k+1}(a)`: n-j+1 subset products on level j, O(n^2 2^n) in total. ac: nothing. |
| `power_projection_of_set_power_series` | 10 | `set_power_series` | 1187 | — | — | The transposed composition: the same levels run backwards, each product replaced by its transpose (reverse, multiply, reverse), then `sum_k W_k i!/(i-k)! s_0^(i-k)` for every m. ac: nothing. |
| `subset_convolution` | 10 | `set_power_series` | 1142 | 2025 | 0.56× | Ranked zeta transform stored as `[mask][rank]` (stride n+1, so each butterfly is a contiguous vector add), per-mask truncated polynomial product, Moebius, read the rank equal to the popcount. ac column: ac has no subset convolution; plain layered `[rank][mask]` version on `ModInt998244353`. |

