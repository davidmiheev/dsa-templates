//! Stress test: factorisations multiply back to the input, the factors are irreducible (checked by
//! exhaustive trial division on small fields) and pairwise distinct.

use poly_mod_p::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn mul_poly(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = (r[i + j] + x * y) % p;
        }
    }
    r
}

/// Exhaustive irreducibility: no monic polynomial of degree 1..=deg/2 divides it.
fn irreducible(f: &[u64], p: u64) -> bool {
    let n = f.len() - 1;
    for d in 1..=n / 2 {
        let count = p.pow(d as u32);
        for code in 0..count {
            let mut g = Vec::with_capacity(d + 1);
            let mut c = code;
            for _ in 0..d {
                g.push(c % p);
                c /= p;
            }
            g.push(1);
            if div_rem(f, &g, p).1.is_empty() {
                return false;
            }
        }
    }
    true
}

#[test]
fn factorisation_is_correct() {
    let mut st = Stress::new("poly_mod_p");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(case as u64);
        let p = *[2u64, 3, 5, 7, 11, 13, 998_244_353].get(rng.below(7)).unwrap();
        // product of a few random monic polynomials, with repeated factors and p-th powers
        let mut f = vec![1u64];
        for _ in 0..rng.range(0, 4) {
            let d = rng.range(1, 4) as usize;
            let mut g: Vec<u64> = (0..d).map(|_| rng.next_u64() % p).collect();
            g.push(1);
            for _ in 0..rng.range(1, 3) {
                f = mul_poly(&f, &g, p);
            }
        }
        if rng.chance(1, 6) && p < 8 {
            // substitute x -> x^p: forces a vanishing derivative
            let mut h = vec![0u64; (f.len() - 1) * p as usize + 1];
            for (i, &c) in f.iter().enumerate() {
                h[i * p as usize] = c;
            }
            f = h;
        }
        let factors = factorize(&f, p);
        let mut back = vec![1u64];
        for (g, e) in &factors {
            assert_eq!(*g.last().unwrap(), 1);
            if g.len() <= 5 && p.pow((g.len() - 1) as u32) <= 20_000 {
                assert!(irreducible(g, p), "case {case}: factor {g:?} over F_{p} is reducible");
            }
            for _ in 0..*e {
                back = mul_poly(&back, g, p);
            }
        }
        assert_eq!(back, f, "case {case} p={p}");
        let mut keys: Vec<&Vec<u64>> = factors.iter().map(|(g, _)| g).collect();
        keys.sort();
        keys.dedup();
        assert_eq!(keys.len(), factors.len(), "factors must be distinct");
        st.case();
    }
    st.done();
}
