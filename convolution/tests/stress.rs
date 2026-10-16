//! Stress test: NTT convolutions against the schoolbook product (mod p and exact u64).

use convolution::{
    convolution, convolution_arbitrary_mod, convolution_large, convolution_mod, convolution_mod_2_64, convolution_three_primes,
    convolution_u64, intt, ntt, P1, P2, P3, P4,
};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn schoolbook_mod(a: &[u32], b: &[u32], m: u32) -> Vec<u32> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut c = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            c[i + j] = (c[i + j] + x as u64 * y as u64) % m as u64;
        }
    }
    c.into_iter().map(|x| x as u32).collect()
}

#[test]
fn matches_schoolbook_for_all_primes() {
    let mut st = Stress::new("convolution/mod");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        let (la, lb) = (rng.range(0, 150) as usize, rng.range(0, 150) as usize);
        for (name, m) in [("p1", P1), ("p2", P2), ("p3", P3)] {
            let a: Vec<u32> = (0..la).map(|_| rng.below(m as usize) as u32).collect();
            let b: Vec<u32> = (0..lb).map(|_| rng.below(m as usize) as u32).collect();
            let got = match m {
                P1 => convolution_mod::<P1>(&a, &b),
                P2 => convolution_mod::<P2>(&a, &b),
                _ => convolution_mod::<P3>(&a, &b),
            };
            assert_eq!(got, schoolbook_mod(&a, &b, m), "{name} la={la} lb={lb}");
        }
        st.case();
        st.ops(1);
    }
    assert_eq!(convolution(&[], &[1]), Vec::<u32>::new());
    st.done();
}

#[test]
fn exact_u64_convolution() {
    let mut st = Stress::new("convolution/u64");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(case as u64 + 5);
        let (la, lb) = (rng.range(1, 120) as usize, rng.range(1, 120) as usize);
        let bound = *[2u64, 1000, 1 << 30].get(case % 3).unwrap();
        let a: Vec<u64> = (0..la).map(|_| rng.next_u64() % bound).collect();
        let b: Vec<u64> = (0..lb).map(|_| rng.next_u64() % bound).collect();
        let mut want = vec![0u128; la + lb - 1];
        for i in 0..la {
            for j in 0..lb {
                want[i + j] += a[i] as u128 * b[j] as u128;
            }
        }
        if want.iter().any(|&x| x >= 1 << 64) {
            continue;
        }
        let want: Vec<u64> = want.into_iter().map(|x| x as u64).collect();
        assert_eq!(convolution_u64(&a, &b), want, "case {case}");
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn large_sizes() {
    let mut st = Stress::new("convolution/large");
    let mut rng = Rng::new(1);
    let n = 1 << 19;
    let a: Vec<u32> = (0..n).map(|_| rng.below(P1 as usize) as u32).collect();
    let b: Vec<u32> = (0..n).map(|_| rng.below(P1 as usize) as u32).collect();
    st.case();
    let c = convolution(&a, &b);
    assert_eq!(c.len(), 2 * n - 1);
    // Spot-check coefficients against the definition.
    for &k in &[0usize, 1, 12345, n - 1, n, 2 * n - 2] {
        let lo = k.saturating_sub(n - 1);
        let hi = k.min(n - 1);
        let want = (lo..=hi).fold(0u64, |acc, i| (acc + a[i] as u64 * b[k - i] as u64) % P1 as u64);
        assert_eq!(c[k] as u64, want, "k={k}");
    }
    // Exact counting: squares of 0/1 vectors.
    let ones: Vec<u64> = (0..200_000).map(|i| (i % 3 == 0) as u64).collect();
    let sq = convolution_u64(&ones, &ones);
    assert_eq!(sq[0], 1);
    assert_eq!(sq.iter().sum::<u64>(), ones.iter().sum::<u64>().pow(2));
    st.ops(n);
    st.done();
}

#[test]
#[should_panic(expected = "2-adic order")]
fn too_long_for_p2_panics() {
    let a = vec![1u32; 1 << 25];
    convolution_mod::<P2>(&a, &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);
}

#[test]
fn every_transform_size_against_schoolbook() {
    let mut st = Stress::new("convolution/sizes");
    let mut rng = Rng::new(77);
    // all power-of-two transform sizes 2^0 .. 2^11 exercise radix-2 and radix-4 stages
    for h in 0..=11 {
        let n = 1usize << h;
        for _ in 0..3 {
            let a: Vec<u32> = (0..n).map(|_| rng.below(P1 as usize) as u32).collect();
            let mut f = a.clone();
            ntt::<P1>(&mut f);
            intt::<P1>(&mut f);
            assert_eq!(f, a, "ntt/intt roundtrip at 2^{h}");
            // transform of a delta is constant, of a shifted delta a character
            let mut delta = vec![0u32; n];
            delta[0] = 1;
            ntt::<P1>(&mut delta);
            assert!(delta.iter().all(|&x| x == 1));
        }
        st.case();
    }
    for case in 0..300 * scale() {
        let (la, lb) = (rng.range(25, 700) as usize, rng.range(25, 700) as usize);
        for (name, m) in [("p4", P4), ("p1", P1)] {
            let a: Vec<u32> = (0..la).map(|_| rng.below(m as usize) as u32).collect();
            let b: Vec<u32> = (0..lb).map(|_| rng.below(m as usize) as u32).collect();
            let got = if m == P4 { convolution_mod::<P4>(&a, &b) } else { convolution(&a, &b) };
            assert_eq!(got, schoolbook_mod(&a, &b, m), "{name} case {case} la={la} lb={lb}");
        }
        st.ops(1);
    }
    st.done();
}

#[test]
fn arbitrary_modulus_and_large() {
    let mut st = Stress::new("convolution/arbitrary");
    let mut rng = Rng::new(5);
    for case in 0..200 * scale() {
        let m = *[1_000_000_007u32, 2, 998_244_353, 1_999_999_973, 65_537].get(case % 5).unwrap();
        let (la, lb) = (rng.range(1, 300) as usize, rng.range(1, 300) as usize);
        let a: Vec<u32> = (0..la).map(|_| rng.below(m as usize) as u32).collect();
        let b: Vec<u32> = (0..lb).map(|_| rng.below(m as usize) as u32).collect();
        assert_eq!(convolution_arbitrary_mod(&a, &b, m), schoolbook_mod(&a, &b, m), "case {case} m={m}");
        let a2: Vec<u32> = a.iter().map(|&x| x % P1).collect();
        let b2: Vec<u32> = b.iter().map(|&x| x % P1).collect();
        // the three-prime route on the large primes
        assert_eq!(convolution_three_primes::<P2, P3, P4>(&a2, &b2, P1), convolution(&a2, &b2), "case {case}");
        st.case();
    }
    st.done();
}

/// Result longer than 2^23: all-ones inputs give a known triangle. About 25 s, so it only runs
/// with `STRESS_SCALE >= 2`.
#[test]
fn convolution_beyond_two_to_the_23() {
    if scale() < 2 {
        return;
    }
    let mut st = Stress::new("convolution/2^24");
    let n = (1usize << 23) + 5;
    let ones = vec![1u32; n];
    let c = convolution_large(&ones, &ones);
    assert_eq!(c.len(), 2 * n - 1);
    for &k in &[0usize, 1, 100, n - 1, n, n + 77, 2 * n - 2] {
        let want = (if k < n { k + 1 } else { 2 * n - 1 - k }) as u32 % P1;
        assert_eq!(c[k], want, "k={k}");
    }
    // random data, coefficients spot-checked against the definition
    let mut rng = Rng::new(3);
    let a: Vec<u32> = (0..(1usize << 23) + 5).map(|_| rng.below(P1 as usize) as u32).collect();
    let b: Vec<u32> = (0..(1usize << 22) + 9).map(|_| rng.below(P1 as usize) as u32).collect();
    let c = convolution_large(&a, &b);
    assert_eq!(c.len(), a.len() + b.len() - 1);
    for &k in &[0usize, 17, (1 << 22) - 1, 1 << 22, (1 << 23) + 4, (1 << 23) + (1 << 22), c.len() - 1] {
        let lo = k.saturating_sub(b.len() - 1);
        let hi = k.min(a.len() - 1);
        let want = (lo..=hi).fold(0u64, |acc, i| (acc + a[i] as u64 * b[k - i] as u64) % P1 as u64);
        assert_eq!(c[k] as u64, want, "k={k}");
    }
    st.done();
}

#[test]
fn mod_2_64_matches_wrapping_schoolbook() {
    let mut st = Stress::new("convolution/2^64");
    let mut rng = Rng::new(9);
    for case in 0..200 * scale() {
        let (la, lb) = (rng.range(1, 200) as usize, rng.range(1, 200) as usize);
        let a: Vec<u64> = (0..la).map(|_| rng.next_u64()).collect();
        let b: Vec<u64> = (0..lb).map(|_| rng.next_u64()).collect();
        let mut want = vec![0u64; la + lb - 1];
        for i in 0..la {
            for j in 0..lb {
                want[i + j] = want[i + j].wrapping_add(a[i].wrapping_mul(b[j]));
            }
        }
        assert_eq!(convolution_mod_2_64(&a, &b), want, "case {case}");
        st.case();
    }
    // all-ones limbs maximise the per-coefficient sums
    let a = vec![u64::MAX; 100_000];
    let c = convolution_mod_2_64(&a, &a);
    for &k in &[0usize, 5, 99_999, 150_000, 199_998] {
        let terms = (if k < 100_000 { k + 1 } else { 199_999 - k }) as u64;
        assert_eq!(c[k], terms, "k={k}"); // (-1) * (-1) = 1 per term
    }
    st.done();
}
