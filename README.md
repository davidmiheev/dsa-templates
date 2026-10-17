# dsa-templates
Collection of data structures and algos in Python3 and Rust

## Stress tests

Every Rust crate has a randomized stress test in `<crate>/tests/stress.rs` that checks it against a naive reference model (deterministic seeds, so a failure names the case that reproduces it) and prints a one-line summary:

```console
$ cargo test --release --workspace --test stress -- --nocapture --test-threads=1
$ STRESS_SCALE=10 cargo test --release --workspace --test stress   # ten times as many cases
```

Shared helpers (RNG, summary printer) live in `stress_common.rs`. The largest cases run on smaller inputs in debug builds so plain `cargo test` stays quick.
