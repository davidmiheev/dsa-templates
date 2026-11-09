//! Stress test: chromatic number and polynomial against enumeration of all colourings.

use graph_coloring::{chromatic_number, chromatic_polynomial};
use modint::ModInt998244353 as Mint;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn colourings(n: usize, edges: &[(usize, usize)], k: usize) -> u64 {
    let mut count = 0;
    let mut color = vec![0usize; n];
    let total = (k as u64).pow(n as u32);
    for code in 0..total {
        let mut c = code;
        for x in color.iter_mut() {
            *x = (c % k as u64) as usize;
            c /= k as u64;
        }
        if edges.iter().all(|&(u, v)| color[u] != color[v]) {
            count += 1;
        }
    }
    count
}

#[test]
fn polynomial_matches_colouring_counts() {
    let mut st = Stress::new("graph_coloring/poly");
    for case in 0..250 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 6) as usize;
        let m = rng.range(0, 9) as usize;
        let edges: Vec<(usize, usize)> = (0..m)
            .map(|_| {
                let u = rng.below(n);
                // occasional loops, frequent parallel edges
                if rng.chance(1, 12) { (u, u) } else { (u, rng.below(n)) }
            })
            .collect();
        let p = chromatic_polynomial::<998_244_353>(n, &edges);
        assert_eq!(p.len(), n + 1);
        for k in 0..=n + 1 {
            let mut value = Mint::zero();
            for &c in p.iter().rev() {
                value = value * Mint::from(k) + c;
            }
            let expect = if k == 0 { 0 } else { colourings(n, &edges, k) };
            assert_eq!(value, Mint::new(expect), "case {case}: n={n} {edges:?} k={k}");
        }
        st.case();
    }
    st.done();
}

#[test]
fn number_matches_smallest_colourable_k() {
    let mut st = Stress::new("graph_coloring/number");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64 + 31);
        let n = rng.range(1, 7) as usize;
        let edges: Vec<(usize, usize)> = (0..rng.range(0, 14))
            .map(|_| {
                let u = rng.below(n);
                (u, (u + 1 + rng.below(n.max(2) - 1)) % n)
            })
            .filter(|&(u, v)| u != v)
            .collect();
        let expect = (1..=n).find(|&k| colourings(n, &edges, k) > 0).unwrap();
        assert_eq!(chromatic_number(n, &edges), expect, "case {case}: n={n} {edges:?}");
        st.case();
    }
    st.done();
}
