# Library Checker `big_integer`: dsa-templates vs ac-library-rs

All 6 problems of the [Library Checker](https://judge.yosupo.jp/) `big_integer` category are solved
with the new `bigint` crate and every solution passes the problem's own checker on every generated
test. ac-library-rs has **no big integers**, so the comparison is against the reference C++
solution that ships with each problem (`sol/correct.cpp`, timed on the same test and machine),
quoted in the comments.

## Method

* Same I/O layer (`FastScan` / `Out` from `template/rust`), `cargo build --release`.
* For each problem every test is run once and checked with the problem's checker; the heaviest
  test is then re-run 5 times. The table reports the **minimum** wall-clock time in milliseconds,
  input parsing and output printing included.
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10% (the reference
  solutions vary by more than that between runs).

## Summary

ac-library-rs has nothing for this category. Against the reference C++ solutions: within 10% or faster on 6 of 6 problems; geometric mean of ours / reference = **0.39**.

New crate `bigint`: `BigUint<BASE>` / `BigInt<BASE>` for any `BASE <= 2^20` with `Decimal` (`10^6`)
and `Hex` (`16^5`) aliases, parsing and printing, addition, subtraction, NTT multiplication and
Newton division with remainder. Randomized stress tests against `i128`, a schoolbook product and
`q b + r` reconstruction, with lengths straddling the algorithm thresholds.

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `addition_of_big_integers` | 5 | `bigint` | 38 | — | — | `bigint::Decimal::add` (signed, base `10^6` limbs; magnitude add or subtract after a comparison). ac: nothing. Reference C++: 70 ms. |
| `addition_of_hex_big_integers` | 5 | `bigint` | 47 | — | — | `bigint::Hex::add` (base `16^5` limbs). ac: nothing. Reference C++: 97 ms. |
| `division_of_big_integers` | 10 | `bigint` | 418 | — | — | `bigint::Decimal::divrem`: Knuth's algorithm D for short operands or short quotients, otherwise a Newton reciprocal (precision doubling), Barrett-style quotient estimate and an exact correction by a product; long quotients in chunks of the divisor length sharing one reciprocal. `O(M(n))`. ac: nothing. Reference C++: 424 ms. |
| `division_of_hex_big_integers` | 10 | `bigint` | 375 | — | — | `bigint::Hex::divrem`, the same code in base `16^5`. ac: nothing. Reference C++: 1.71 s. |
| `multiplication_of_big_integers` | 5 | `bigint` | 114 | — | — | `bigint::Decimal::mul`: schoolbook below 32 limbs, otherwise an exact convolution modulo two NTT primes (`998244353`, `469762049`) and a CRT; the base `10^6` keeps every coefficient below `P1 P2` for operands up to about 2.8 million digits (a three-prime fallback beyond). ac: nothing (its convolution is modulo one prime). Reference C++: 126 ms; most of our time is the six transforms of size `2^20`. |
| `multiplication_of_hex_big_integers` | 5 | `bigint` | 142 | — | — | `bigint::Hex::mul`, the same code in base `16^5`. ac: nothing. Reference C++: 2.15 s. |

