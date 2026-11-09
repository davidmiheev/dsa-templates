//! Stress test: field axioms, the Cantor basis, and the additive-FFT convolution against schoolbook.

use gf2_64::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn field_arithmetic() {
    let mut st = Stress::new("gf2_64/field");
    let mut rng = Rng::new(1);
    for _ in 0..2000 * scale() {
        let (a, b, c) = (rng.next_u64(), rng.next_u64(), rng.next_u64());
        assert_eq!(mul(a, b), mul_portable(a, b));
        assert_eq!(mul(a, b), mul(b, a));
        assert_eq!(mul(a, mul(b, c)), mul(mul(a, b), c));
        assert_eq!(mul(a, b ^ c), mul(a, b) ^ mul(a, c));
        assert_eq!(mul(a, 1), a);
        st.ops(1);
    }
    assert_eq!(mul(1 << 63, 2), 0b11011);
    st.case();
    st.done();
}

#[test]
fn cantor_basis_is_a_chain() {
    let b = cantor_basis(40);
    assert_eq!(b[0], 1);
    for i in 0..39 {
        assert_eq!(square(b[i + 1]) ^ b[i + 1], b[i], "i={i}");
    }
    // independence: no non-empty subset XORs to zero (check via Gaussian elimination)
    let mut basis: Vec<u64> = Vec::new();
    for &v in &b {
        let mut x = v;
        for &r in &basis {
            if x & (r & r.wrapping_neg()) != 0 {
                x ^= r;
            }
        }
        assert_ne!(x, 0, "dependent basis");
        basis.push(x);
    }
}

#[test]
fn convolution_matches_schoolbook() {
    let mut st = Stress::new("gf2_64/convolution");
    let mut rng = Rng::new(2);
    for case in 0..150 * scale() {
        let (la, lb) = (rng.range(1, 70) as usize, rng.range(1, 70) as usize);
        let a: Vec<u64> = (0..la).map(|_| rng.next_u64()).collect();
        let b: Vec<u64> = (0..lb).map(|_| rng.next_u64()).collect();
        let mut want = vec![0u64; la + lb - 1];
        for i in 0..la {
            for j in 0..lb {
                want[i + j] ^= mul_portable(a[i], b[j]);
            }
        }
        assert_eq!(convolution(&a, &b), want, "case {case} la={la} lb={lb}");
        st.case();
    }
    // power-of-two-ish sizes around the transform boundaries
    for &(la, lb) in &[(16usize, 17usize), (64, 65), (128, 129), (300, 213)] {
        let a: Vec<u64> = (0..la).map(|_| rng.next_u64()).collect();
        let b: Vec<u64> = (0..lb).map(|_| rng.next_u64()).collect();
        let mut want = vec![0u64; la + lb - 1];
        for i in 0..la {
            for j in 0..lb {
                want[i + j] ^= mul(a[i], b[j]);
            }
        }
        assert_eq!(convolution(&a, &b), want, "la={la} lb={lb}");
    }
    st.done();
}
