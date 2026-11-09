//! Stress test: GCD/LCM convolutions and the transforms against the definitions.

use divisor_transform::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const M: u32 = 998_244_353;

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 { a } else { gcd(b, a % b) }
}

#[test]
fn matches_definitions() {
    let mut st = Stress::new("divisor_transform");
    for case in 0..300 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 60) as usize;
        let mut a: Vec<u32> = (0..=n).map(|_| rng.below(M as usize) as u32).collect();
        let mut b: Vec<u32> = (0..=n).map(|_| rng.below(M as usize) as u32).collect();
        a[0] = 0;
        b[0] = 0;
        let (mut g, mut l) = (vec![0u64; n + 1], vec![0u64; n + 1]);
        for i in 1..=n {
            for j in 1..=n {
                let prod = a[i] as u64 * b[j] as u64 % M as u64;
                let d = gcd(i, j);
                g[d] = (g[d] + prod) % M as u64;
                let m = i / d * j;
                if m <= n {
                    l[m] = (l[m] + prod) % M as u64;
                }
            }
        }
        assert_eq!(gcd_convolution::<M>(&a, &b), g.iter().map(|&x| x as u32).collect::<Vec<_>>(), "gcd n={n}");
        assert_eq!(lcm_convolution::<M>(&a, &b), l.iter().map(|&x| x as u32).collect::<Vec<_>>(), "lcm n={n}");
        let mut z = a.clone();
        zeta_divisors::<M>(&mut z);
        for k in 1..=n {
            let want = (1..=k).filter(|d| k % d == 0).fold(0u64, |s, d| (s + a[d] as u64) % M as u64);
            assert_eq!(z[k] as u64, want);
        }
        mobius_divisors::<M>(&mut z);
        assert_eq!(z, a);
        zeta_multiples::<M>(&mut z);
        mobius_multiples::<M>(&mut z);
        assert_eq!(z, a);
        st.case();
        st.ops(n);
    }
    st.done();
}
