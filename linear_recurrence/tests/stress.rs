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
