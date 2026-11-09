# Library Checker `sample`: dsa-templates vs ac-library-rs

The 3 problems of the [Library Checker](https://judge.yosupo.jp/) `sample` category (`A + B` with
one case, many cases, and signed 128-bit operands) are solved with the `bigint` crate and pass
the problem's own checker on every generated test. ac-library-rs has nothing for these, so the
comparison is against the reference C++ solution that ships with each problem (`sol/correct.cpp`,
timed on the same test and machine), quoted in the comments. They mostly measure input and output.

## Method

* Same I/O layer (`FastScan` / `Out` from `template/rust`), `cargo build --release`.
* For each problem every test is run once and checked with the problem's checker; the heaviest
  test is then re-run 5 times. The table reports the **minimum** wall-clock time in milliseconds,
  input parsing and output printing included.
* Machine: shared 4-core cloud container, so run-to-run noise is around ±10%.

## Summary

ac-library-rs has nothing for this category. Against the reference C++ solutions: within 10% or faster on 3 of 3 problems; geometric mean of ours / reference = **0.52**.

## Results

| Problem | TL (s) | Our lib | ours (ms) | ac (ms) | ours / ac | Comment |
|---|---:|---|---:|---:|---:|---|
| `aplusb` | 2 | `bigint` | 1 | — | — | `bigint::Decimal::add`: nothing here needs a library, the crate is used only to follow the rule that every solution goes through one. ac: nothing. Reference C++: 3 ms. |
| `many_aplusb` | 5 | `bigint` | 195 | — | — | `bigint::Decimal::add` per case (sums up to `2 * 10^18`). ac: nothing. Reference C++: 260 ms (`scanf` / `printf`). |
| `many_aplusb_128bit` | 5 | `bigint` | 130 | — | — | `bigint::Decimal::add`: signed operands up to `10^37` in magnitude, which is exactly what the signed wrapper does. ac: nothing. Reference C++: 325 ms. |

