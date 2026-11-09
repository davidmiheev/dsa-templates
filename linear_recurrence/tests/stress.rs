//! Stress test: BM recovers the length of random recurrences and predicts the sequence.

use linear_recurrence::berlekamp_massey;
use modint::ModInt998244353 as Mint;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn recovers_random_recurrences() {
    let mut st = Stress::new("linear_recurrence");
    for case in 0..2000 * scale() {
        let mut rng = Rng::new(case as u64);
        let l = rng.below(8);
        let c: Vec<Mint> = (0..l).map(|_| Mint::new(rng.next_u64())).collect();
        let n = 2 * l + rng.below(10);
        let mut s: Vec<Mint> = (0..l).map(|_| Mint::new(rng.next_u64())).collect();
        while s.len() < n {
            let k = s.len();
            s.push((0..l).fold(Mint::zero(), |a, i| a + c[i] * s[k - 1 - i]));
        }
        let r = berlekamp_massey(&s);
        assert!(r.len() <= l, "found a longer recurrence: {} > {}", r.len(), l);
        // it reproduces the sequence
        for k in r.len()..s.len() {
            let v = (0..r.len()).fold(Mint::zero(), |a, i| a + r[i] * s[k - 1 - i]);
            assert_eq!(v, s[k], "case {case}");
        }
        st.case();
        st.ops(n * l.max(1));
    }
    // the zero sequence and a constant one
    assert!(berlekamp_massey::<998_244_353>(&[Mint::zero(); 5]).is_empty());
    assert_eq!(berlekamp_massey(&[Mint::new(7); 5]), vec![Mint::new(1)]);
    st.done();
}

#[test]
fn kth_and_consecutive_terms() {
    use linear_recurrence::{consecutive_terms, kth_term};
    let mut st = Stress::new("linear_recurrence/terms");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(5_000 + case as u64);
        let d = 1 + rng.below(9);
        let a: Vec<Mint> = (0..d).map(|_| Mint::new(rng.next_u64())).collect();
        let c: Vec<Mint> = (0..d).map(|_| Mint::new(rng.next_u64())).collect();
        let total = 400;
        let mut s = a.clone();
        while s.len() < total {
            let k = s.len();
            s.push((0..d).fold(Mint::zero(), |acc, j| acc + c[j] * s[k - 1 - j]));
        }
        let k = rng.below(300) as u64;
        assert_eq!(kth_term(&a, &c, k), s[k as usize], "k={k} d={d}");
        let m = 1 + rng.below(60);
        let got = consecutive_terms(&a, &c, k, m);
        assert_eq!(got, s[k as usize..k as usize + m].to_vec(), "consecutive k={k} m={m} d={d}");
        st.case();
        st.ops(d * 300);
    }
    // larger orders exercise the spectral loop
    for case in 0..30 * scale() {
        let mut rng = Rng::new(6_000 + case as u64);
        let d = 60 + rng.below(150);
        let a: Vec<Mint> = (0..d).map(|_| Mint::new(rng.next_u64())).collect();
        let c: Vec<Mint> = (0..d).map(|_| Mint::new(rng.next_u64())).collect();
        let mut s = a.clone();
        while s.len() < 3000 {
            let k = s.len();
            s.push((0..d).fold(Mint::zero(), |acc, j| acc + c[j] * s[k - 1 - j]));
        }
        for _ in 0..10 {
            let k = rng.below(3000);
            assert_eq!(kth_term(&a, &c, k as u64), s[k], "k={k} d={d}");
        }
        st.case();
    }
    // huge indices: the two functions and shifts of the start agree
    let mut rng = Rng::new(1);
    for _ in 0..20 {
        let d = 1 + rng.below(12);
        let a: Vec<Mint> = (0..d).map(|_| Mint::new(rng.next_u64())).collect();
        let c: Vec<Mint> = (0..d).map(|_| Mint::new(rng.next_u64())).collect();
        let k = rng.next_u64() >> 2;
        let block = consecutive_terms(&a, &c, k, 20);
        for (i, &v) in block.iter().enumerate() {
            assert_eq!(kth_term(&a, &c, k + i as u64), v);
        }
    }
    st.done();
}
