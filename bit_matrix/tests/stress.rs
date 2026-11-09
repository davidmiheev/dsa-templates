//! Stress test: every operation against the definition on small matrices.

use bit_matrix::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn random(rng: &mut Rng, r: usize, c: usize, density: u64) -> BitMatrix {
    let mut m = BitMatrix::new(r, c);
    for i in 0..r {
        for j in 0..c {
            m.set(i, j, rng.chance(density, 8));
        }
    }
    m
}

fn naive_mul(a: &BitMatrix, b: &BitMatrix) -> BitMatrix {
    let mut c = BitMatrix::new(a.rows(), b.cols());
    for i in 0..a.rows() {
        for j in 0..b.cols() {
            let v = (0..a.cols()).fold(false, |s, k| s ^ (a.get(i, k) & b.get(k, j)));
            c.set(i, j, v);
        }
    }
    c
}

/// Rank by the most naive elimination on a `Vec<Vec<bool>>`.
fn naive_rank(m: &BitMatrix) -> usize {
    let mut a: Vec<Vec<bool>> = (0..m.rows()).map(|i| (0..m.cols()).map(|j| m.get(i, j)).collect()).collect();
    let mut r = 0;
    for c in 0..m.cols() {
        if let Some(p) = (r..a.len()).find(|&i| a[i][c]) {
            a.swap(r, p);
            for i in 0..a.len() {
                if i != r && a[i][c] {
                    let pivot_row = a[r].clone();
                    for (x, y) in a[i].iter_mut().zip(&pivot_row) {
                        *x ^= y;
                    }
                }
            }
            r += 1;
        }
    }
    r
}

#[test]
fn matrix_operations_match_definitions() {
    let mut st = Stress::new("bit_matrix");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let dim = |rng: &mut Rng| [1usize, 2, 3, 7, 8, 9, 63, 64, 65, 70, 130][rng.below(11)];
        let (n, m, k) = (dim(&mut rng), dim(&mut rng), dim(&mut rng));
        let (n, m, k) = if case % 3 == 0 { (n.min(20), m.min(20), k.min(20)) } else { (n, m, k) };
        let density = 1 + rng.below(7) as u64;
        let a = random(&mut rng, n, m, density);
        let b = random(&mut rng, m, k, density);
        assert_eq!(a.mul(&b), naive_mul(&a, &b));
        assert_eq!(a.transpose().transpose(), a);
        assert_eq!(a.rank(), naive_rank(&a));
        assert_eq!(a.rank(), a.transpose().rank());
        // square: determinant, inverse
        let s = random(&mut rng, n, n, density);
        let inv = s.inverse();
        if s.det() {
            assert_eq!(s.mul(&inv), BitMatrix::identity(n));
            assert_eq!(inv.mul(&s), BitMatrix::identity(n));
        } else {
            assert_eq!(inv.rows(), 0);
        }
        // linear system: A x = b with b from a random x (always solvable) or random b
        let x0: Vec<u64> = (0..m.div_ceil(64)).map(|_| rng.next_u64()).collect();
        let mut rhs = vec![0u64; n.div_ceil(64)];
        for i in 0..n {
            let v = (0..m).fold(false, |s, j| s ^ (a.get(i, j) & (x0[j / 64] >> (j % 64) & 1 == 1)));
            if v {
                rhs[i / 64] |= 1 << (i % 64);
            }
        }
        if case % 2 == 1 {
            rhs.iter_mut().for_each(|w| *w = rng.next_u64());
        }
        let (part, basis) = a.solve(&rhs);
        let apply = |x: &[u64]| -> Vec<u64> {
            let mut y = vec![0u64; n.div_ceil(64)];
            for i in 0..n {
                let v = (0..m).fold(false, |s, j| s ^ (a.get(i, j) & (x[j / 64] >> (j % 64) & 1 == 1)));
                if v {
                    y[i / 64] |= 1 << (i % 64);
                }
            }
            y
        };
        if part.is_empty() {
            // inconsistent: rank of [A | b] exceeds rank of A
            let mut aug = BitMatrix::new(n, m + 1);
            for i in 0..n {
                for j in 0..m {
                    aug.set(i, j, a.get(i, j));
                }
                aug.set(i, m, rhs[i / 64] >> (i % 64) & 1 == 1);
            }
            assert!(aug.rank() > a.rank());
        } else {
            let mut y = apply(&part);
            // mask padding bits of rhs beyond n
            if n % 64 != 0 {
                *y.last_mut().unwrap() &= (1u64 << (n % 64)) - 1;
            }
            let mut want = rhs.clone();
            if n % 64 != 0 {
                *want.last_mut().unwrap() &= (1u64 << (n % 64)) - 1;
            }
            assert_eq!(y, want);
            assert_eq!(basis.len(), m - a.rank());
            for v in &basis {
                assert!(apply(v).iter().all(|&w| w == 0));
            }
            assert_eq!(BitMatrix::new(0, 0).rows(), 0);
            // basis is independent: its rank equals its size
            let mut bm = BitMatrix::new(basis.len(), m);
            for (i, v) in basis.iter().enumerate() {
                for j in 0..m {
                    bm.set(i, j, v[j / 64] >> (j % 64) & 1 == 1);
                }
            }
            assert_eq!(bm.rank(), basis.len());
        }
        st.case();
        st.ops(n * m * k / 64 + 1);
    }
    st.done();
}

fn span(basis: &[u32]) -> std::collections::HashSet<u32> {
    let mut s = std::collections::HashSet::from([0u32]);
    for &b in basis {
        let add: Vec<u32> = s.iter().map(|&x| x ^ b).collect();
        s.extend(add);
    }
    s
}

#[test]
fn space_intersection_matches_enumeration() {
    let mut st = Stress::new("bit_matrix/intersect");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(900 + case as u64);
        let bits = 1 + rng.below(7);
        let pick = |rng: &mut Rng| -> Vec<u32> {
            let mut basis: Vec<u32> = Vec::new();
            for _ in 0..rng.below(bits + 1) {
                let x = rng.below(1 << bits) as u32;
                if !span(&basis).contains(&x) {
                    basis.push(x);
                }
            }
            basis
        };
        let (u, v) = (pick(&mut rng), pick(&mut rng));
        let w = intersect_spaces(&u, &v);
        let (su, sv) = (span(&u), span(&v));
        let inter: std::collections::HashSet<u32> = su.intersection(&sv).copied().collect();
        let sw = span(&w);
        assert_eq!(sw, inter, "u={u:?} v={v:?} w={w:?}");
        assert_eq!(sw.len(), 1 << w.len(), "w must be independent");
        st.case();
    }
    // 30-bit values do not overflow
    assert_eq!(intersect_spaces(&[(1 << 30) - 1], &[(1 << 30) - 1]), vec![(1 << 30) - 1]);
    st.done();
}
