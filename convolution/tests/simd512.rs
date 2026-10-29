//! The AVX-512 transforms must reproduce the scalar ones word for word.
#![cfg(all(target_arch = "x86_64", feature = "avx512"))]

use convolution::simd512::{available, intt_unscaled_avx512, ntt_avx512};
use convolution::{intt_unscaled_scalar, ntt_scalar, P1, P3};

fn lcg(seed: &mut u64) -> u64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *seed >> 33
}

fn check<const P: u32>() {
    let mut seed = 777u64;
    for h in 0..=20usize {
        let a: Vec<u32> = (0..1usize << h).map(|_| (lcg(&mut seed) % P as u64) as u32).collect();
        let (mut x, mut y) = (a.clone(), a.clone());
        ntt_scalar::<P>(&mut x);
        ntt_avx512::<P>(&mut y);
        assert_eq!(x, y, "forward, P={P}, h={h}");
        intt_unscaled_scalar::<P>(&mut x);
        intt_unscaled_avx512::<P>(&mut y);
        assert_eq!(x, y, "inverse, P={P}, h={h}");
        let (mut u, mut v) = (a.clone(), a);
        intt_unscaled_scalar::<P>(&mut u);
        intt_unscaled_avx512::<P>(&mut v);
        assert_eq!(u, v, "inverse of raw data, P={P}, h={h}");
    }
}

#[test]
fn avx512_equals_scalar() {
    if !available() {
        eprintln!("AVX-512 not available, skipping");
        return;
    }
    check::<P1>();
    check::<P3>();
    check::<2_013_265_921>();
}
