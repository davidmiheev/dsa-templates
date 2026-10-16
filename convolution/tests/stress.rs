//! Stress test: NTT convolutions against the schoolbook product (mod p and exact u64).

use convolution::{convolution, convolution_mod, convolution_u64, P1, P2, P3};

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
#[should_panic(expected = "result too long")]
fn too_long_for_p2_panics() {
    let a = vec![1u32; 1 << 25];
    convolution_mod::<P2>(&a, &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);
}
