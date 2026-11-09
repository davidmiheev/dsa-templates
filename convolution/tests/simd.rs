//! The AVX2 transforms must reproduce the scalar ones word for word.
#![cfg(target_arch = "x86_64")]

use convolution::simd::{available, intt_unscaled_avx2, ntt_avx2};
use convolution::{intt_unscaled, ntt, P1, P2, P3};

fn lcg(seed: &mut u64) -> u64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *seed >> 33
}

fn check<const P: u32>() {
    let mut seed = 12345u64;
    for h in 0..=20usize {
        let n = 1usize << h;
        let a: Vec<u32> = (0..n).map(|_| (lcg(&mut seed) % P as u64) as u32).collect();
        let (mut x, mut y) = (a.clone(), a.clone());
        ntt::<P>(&mut x);
        ntt_avx2::<P>(&mut y);
        assert_eq!(x, y, "forward, P={P}, h={h}");
        intt_unscaled::<P>(&mut x);
        intt_unscaled_avx2::<P>(&mut y);
        assert_eq!(x, y, "inverse, P={P}, h={h}");
        // also from an arbitrary (not transformed) input
        let (mut u, mut v) = (a.clone(), a);
        intt_unscaled::<P>(&mut u);
        intt_unscaled_avx2::<P>(&mut v);
        assert_eq!(u, v, "inverse of raw data, P={P}, h={h}");
    }
}

#[test]
fn avx2_equals_scalar() {
    if !available() {
        eprintln!("AVX2 not available, skipping");
        return;
    }
    check::<P1>();
    check::<P2>();
    check::<P3>();
    check::<1_004_535_809>();
    check::<2_013_265_921>();
}
