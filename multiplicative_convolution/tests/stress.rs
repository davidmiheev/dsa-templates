//! Stress test: multiplicative convolutions against the double loop.

use multiplicative_convolution::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const M: u64 = 998_244_353;

fn naive(m: usize, a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut c = vec![0u64; m];
    for i in 0..m {
        for j in 0..m {
            let k = i * j % m;
            c[k] = (c[k] + a[i] as u64 * b[j] as u64) % M;
        }
    }
    c.into_iter().map(|x| x as u32).collect()
}

#[test]
fn prime_modulus_matches_double_loop() {
    let mut st = Stress::new("multiplicative/prime");
    for &p in &[2usize, 3, 5, 7, 11, 13, 17, 31, 61, 101, 257, 997, 1009] {
        let mut rng = Rng::new(p as u64);
        for _ in 0..4 * scale() {
            let a: Vec<u32> = (0..p).map(|_| rng.below(M as usize) as u32).collect();
            let b: Vec<u32> = (0..p).map(|_| rng.below(M as usize) as u32).collect();
            assert_eq!(mul_mod_prime_convolution(p, &a, &b), naive(p, &a, &b), "p={p}");
            st.case();
        }
    }
    st.done();
}

#[test]
fn power_of_two_modulus_matches_double_loop() {
    let mut st = Stress::new("multiplicative/pow2");
    for n in 0..=10usize {
        let mut rng = Rng::new(n as u64 + 100);
        for _ in 0..4 * scale() {
            let m = 1usize << n;
            let a: Vec<u32> = (0..m).map(|_| rng.below(M as usize) as u32).collect();
            let b: Vec<u32> = (0..m).map(|_| rng.below(M as usize) as u32).collect();
            assert_eq!(mul_mod_pow2_convolution(n, &a, &b), naive(m, &a, &b), "n={n}");
            st.case();
        }
    }
    st.done();
}
