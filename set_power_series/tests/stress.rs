//! Stress test: transforms and convolutions against the definitions.

use set_power_series::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const M: u32 = 998_244_353;

fn naive(a: &[u32], b: &[u32], ok: impl Fn(usize, usize) -> Option<usize>) -> Vec<u32> {
    let mut c = vec![0u64; a.len()];
    for i in 0..a.len() {
        for j in 0..b.len() {
            if let Some(k) = ok(i, j) {
                c[k] = (c[k] + a[i] as u64 * b[j] as u64) % M as u64;
            }
        }
    }
    c.into_iter().map(|x| x as u32).collect()
}

#[test]
fn convolutions_match_definitions() {
    let mut st = Stress::new("set_power_series");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(0, 7) as usize;
        let len = 1usize << n;
        let a: Vec<u32> = (0..len).map(|_| rng.below(M as usize) as u32).collect();
        let b: Vec<u32> = (0..len).map(|_| rng.below(M as usize) as u32).collect();
        assert_eq!(and_convolution::<M>(&a, &b), naive(&a, &b, |i, j| Some(i & j)));
        assert_eq!(or_convolution::<M>(&a, &b), naive(&a, &b, |i, j| Some(i | j)));
        assert_eq!(xor_convolution::<M>(&a, &b), naive(&a, &b, |i, j| Some(i ^ j)));
        assert_eq!(subset_convolution::<M>(&a, &b), naive(&a, &b, |i, j| if i & j == 0 { Some(i | j) } else { None }));
        // transforms invert each other
        let mut z = a.clone();
        zeta_subset::<M>(&mut z);
        mobius_subset::<M>(&mut z);
        assert_eq!(z, a);
        zeta_superset::<M>(&mut z);
        mobius_superset::<M>(&mut z);
        assert_eq!(z, a);
        st.case();
        st.ops(len);
    }
    st.done();
}
