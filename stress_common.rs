//! Shared helpers for the per-crate stress tests (`<crate>/tests/stress.rs`).
//!
//! Each test pulls this file in with
//! `#[path = "../../stress_common.rs"] mod common;`
//! so no crate needs an extra dependency.
//!
//! Every stress test runs many small random cases against a naive reference
//! model (deterministic seeds, so a failure names the case that reproduces it),
//! and prints a one-line summary. Run them with
//!
//! ```text
//! cargo test --release --workspace -- --nocapture --test-threads=1
//! ```
//!
//! Set `STRESS_SCALE=10` to run ten times as many cases.
#![allow(dead_code)]

use std::time::Instant;

/// xorshift64* generator: deterministic and dependency-free.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform in `0..n` (`n > 0`).
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    /// Uniform in `lo..=hi`.
    pub fn range(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next_u64() % (hi - lo + 1) as u64) as i64
    }

    /// True with probability `num / den`.
    pub fn chance(&mut self, num: u64, den: u64) -> bool {
        self.next_u64() % den < num
    }

    pub fn pick<T: Copy>(&mut self, xs: &[T]) -> T {
        xs[self.below(xs.len())]
    }

    /// Random lowercase string over the first `alphabet` letters.
    pub fn string(&mut self, len: usize, alphabet: usize) -> String {
        (0..len).map(|_| (b'a' + self.below(alphabet) as u8) as char).collect()
    }

    /// Random labelled tree on `n` vertices as an adjacency list.
    /// `shape`: 0 = uniform random, 1 = path, 2 = star, 3 = caterpillar-ish.
    pub fn tree(&mut self, n: usize, shape: usize) -> Vec<Vec<usize>> {
        let mut adj = vec![Vec::new(); n];
        for v in 1..n {
            let p = match shape {
                1 => v - 1,
                2 => 0,
                3 => if v % 2 == 1 { v.saturating_sub(2).max(0) } else { v - 1 },
                _ => self.below(v),
            };
            adj[p].push(v);
            adj[v].push(p);
        }
        adj
    }
}

/// Number of random cases is multiplied by this (`STRESS_SCALE` env var, default 1).
pub fn scale() -> usize {
    std::env::var("STRESS_SCALE").ok().and_then(|s| s.parse().ok()).unwrap_or(1)
}

/// Case/operation counter that prints a summary line when finished.
pub struct Stress {
    name: &'static str,
    cases: u64,
    ops: u64,
    start: Instant,
}

impl Stress {
    pub fn new(name: &'static str) -> Self {
        Stress { name, cases: 0, ops: 0, start: Instant::now() }
    }

    pub fn case(&mut self) {
        self.cases += 1;
    }

    pub fn ops(&mut self, k: usize) {
        self.ops += k as u64;
    }

    pub fn done(self) {
        let ms = self.start.elapsed().as_secs_f64() * 1000.0;
        println!(
            "[stress] {:<24} cases={:>7} checks={:>10} time={:>9.1} ms",
            self.name, self.cases, self.ops, ms
        );
    }
}
