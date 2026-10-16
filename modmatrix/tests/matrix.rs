//! Stress tests for `Matrix`, `hafnian`, `sparse_determinant` and `determinant_mod` against
//! naive definitions.

use modmatrix::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const M: u32 = 998_244_353;
type Mat = Matrix<M>;

fn mulm(a: u32, b: u32) -> u32 {
    (a as u64 * b as u64 % M as u64) as u32
}
fn addm(a: u32, b: u32) -> u32 {
    ((a as u64 + b as u64) % M as u64) as u32
}
fn subm(a: u32, b: u32) -> u32 {
    ((a as u64 + M as u64 - b as u64) % M as u64) as u32
}
fn powm(mut b: u32, mut e: u64) -> u32 {
    let mut r = 1;
    while e > 0 {
        if e & 1 == 1 {
            r = mulm(r, b);
        }
        b = mulm(b, b);
        e >>= 1;
    }
    r
}

fn random(rng: &mut Rng, r: usize, c: usize) -> Mat {
    let small = rng.chance(1, 3);
    let v = (0..r * c).map(|_| if small { rng.below(3) as u32 } else { rng.below(M as usize) as u32 }).collect();
    Mat::from_vec(r, c, v)
}

fn naive_mul(a: &Mat, b: &Mat) -> Mat {
    let mut c = Mat::zeros(a.rows(), b.cols());
    for i in 0..a.rows() {
        for j in 0..b.cols() {
            let v = (0..a.cols()).fold(0, |s, k| addm(s, mulm(a.get(i, k), b.get(k, j))));
            c.set(i, j, v);
        }
    }
    c
}

/// Determinant by Laplace expansion along the first row.
fn laplace(a: &Mat) -> u32 {
    let n = a.rows();
    if n == 0 {
        return 1;
    }
    if n == 1 {
        return a.get(0, 0);
    }
    let mut total = 0;
    for j in 0..n {
        let mut minor = Vec::new();
        for i in 1..n {
            for c in (0..n).filter(|&c| c != j) {
                minor.push(a.get(i, c));
            }
        }
        let d = mulm(a.get(0, j), laplace(&Mat::from_vec(n - 1, n - 1, minor)));
        total = if j % 2 == 0 { addm(total, d) } else { subm(total, d) };
    }
    total
}

/// A random matrix of the given rank (product of two thin random matrices).
fn with_rank(rng: &mut Rng, n: usize, m: usize, rank: usize) -> Mat {
    if rank == 0 {
        return Mat::zeros(n, m);
    }
    let full = |rng: &mut Rng, r: usize, c: usize| Mat::from_vec(r, c, (0..r * c).map(|_| 1 + rng.below(M as usize - 1) as u32).collect());
    naive_mul(&full(rng, n, rank), &full(rng, rank, m))
}

#[test]
fn product_power_det_rank_inverse() {
    let mut st = Stress::new("modmatrix/basic");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let dim = |rng: &mut Rng| 1 + rng.below(if case % 50 == 0 { 40 } else { 8 });
        let (n, m, k) = (dim(&mut rng), dim(&mut rng), dim(&mut rng));
        let a = random(&mut rng, n, m);
        let b = random(&mut rng, m, k);
        assert_eq!(a.mul(&b), naive_mul(&a, &b));
        assert_eq!(a.transpose().transpose(), a);
        // power
        let s = random(&mut rng, n, n);
        let e = rng.below(20) as u64;
        let mut want = Mat::identity(n);
        for _ in 0..e {
            want = naive_mul(&want, &s);
        }
        assert_eq!(s.pow(e), want);
        // determinant, rank
        if n <= 6 {
            assert_eq!(s.det(), laplace(&s));
        }
        let r = rng.below(n.min(m) + 1);
        let lr = with_rank(&mut rng, n, m, r);
        assert_eq!(lr.rank(), r, "n={n} m={m} r={r}");
        assert_eq!(lr.transpose().rank(), r);
        // inverse
        let inv = s.inverse();
        if s.det() != 0 {
            assert_eq!(s.mul(&inv), Mat::identity(n));
        } else {
            assert_eq!(inv.rows(), 0);
        }
        let sr = n.saturating_sub(1 + rng.below(2));
        let sing = with_rank(&mut rng, n, n, sr);
        if n >= 1 && sing.rank() < n {
            assert_eq!(sing.inverse().rows(), 0);
            assert_eq!(sing.det(), 0);
        }
        st.case();
        st.ops(n * m * k);
    }
    // a product wider than the column tile
    let mut rng = Rng::new(5);
    let (a, b) = (random(&mut rng, 7, 9), random(&mut rng, 9, 600));
    assert_eq!(a.mul(&b), naive_mul(&a, &b));
    st.done();
}

#[test]
fn adjugate_solve_charpoly() {
    let mut st = Stress::new("modmatrix/adj");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(100_000 + case as u64);
        let n = 1 + rng.below(6);
        let rank = if rng.chance(1, 3) { n } else { rng.below(n + 1) };
        let a = with_rank(&mut rng, n, n, rank);
        // adjugate against the cofactor definition
        let adj = a.adjugate();
        for i in 0..n {
            for j in 0..n {
                let mut minor = Vec::new();
                for r in (0..n).filter(|&r| r != j) {
                    for c in (0..n).filter(|&c| c != i) {
                        minor.push(a.get(r, c));
                    }
                }
                let d = laplace(&Mat::from_vec(n - 1, n - 1, minor));
                let want = if (i + j) % 2 == 0 { d } else { subm(0, d) };
                assert_eq!(adj.get(i, j), want, "adj n={n} rank={rank} ({i},{j})");
            }
        }
        // linear systems
        let m = 1 + rng.below(6);
        let sysr = rng.below(n.min(m) + 1);
        let sys = with_rank(&mut rng, n, m, sysr);
        let x0: Vec<u32> = (0..m).map(|_| rng.below(M as usize) as u32).collect();
        let apply = |x: &[u32]| -> Vec<u32> { (0..n).map(|i| (0..m).fold(0, |s, j| addm(s, mulm(sys.get(i, j), x[j])))).collect() };
        let mut b = apply(&x0);
        if case % 3 == 1 {
            b[rng.below(n)] = rng.below(M as usize) as u32;
        }
        let (part, basis) = sys.solve(&b);
        if part.is_empty() {
            let mut aug = Vec::new();
            for (i, &bi) in b.iter().enumerate() {
                aug.extend_from_slice(sys.row(i));
                aug.push(bi);
            }
            assert!(Mat::from_vec(n, m + 1, aug).rank() > sys.rank());
        } else {
            assert_eq!(apply(&part), b);
            assert_eq!(basis.len(), m - sys.rank());
            for v in &basis {
                assert!(apply(v).iter().all(|&x| x == 0));
            }
            let flat: Vec<u32> = basis.iter().flatten().copied().collect();
            assert_eq!(Mat::from_vec(basis.len(), m, flat).rank(), basis.len());
        }
        // characteristic polynomial at random points
        let cp = a.characteristic_polynomial();
        assert_eq!(cp.len(), n + 1);
        assert_eq!(cp[n], 1);
        for _ in 0..3 {
            let x = rng.below(M as usize) as u32;
            let mut shifted = a.data().iter().map(|&v| subm(0, v)).collect::<Vec<_>>();
            for i in 0..n {
                shifted[i * n + i] = addm(shifted[i * n + i], x);
            }
            let want = Mat::from_vec(n, n, shifted).det();
            let got = cp.iter().rev().fold(0, |s, &c| addm(mulm(s, x), c));
            assert_eq!(got, want, "charpoly n={n}");
        }
        st.case();
    }
    // 0 x 0 and a matrix that needs pivoting in the Hessenberg reduction
    assert_eq!(Mat::zeros(0, 0).characteristic_polynomial(), vec![1]);
    let a = Mat::from_vec(3, 3, vec![0, 0, 1, 1, 0, 0, 0, 1, 0]);
    assert_eq!(a.characteristic_polynomial(), vec![M - 1, 0, 0, 1]);
    st.done();
}

fn pfaffian_expansion(a: &Mat, rows: &[usize]) -> u32 {
    if rows.is_empty() {
        return 1;
    }
    let mut total = 0;
    for t in 1..rows.len() {
        let rest: Vec<usize> = rows.iter().enumerate().filter(|&(k, _)| k != 0 && k != t).map(|(_, &r)| r).collect();
        let d = mulm(a.get(rows[0], rows[t]), pfaffian_expansion(a, &rest));
        total = if t % 2 == 1 { addm(total, d) } else { subm(total, d) };
    }
    total
}

fn hafnian_expansion(a: &Mat, rows: &[usize]) -> u32 {
    if rows.is_empty() {
        return 1;
    }
    let mut total = 0;
    for t in 1..rows.len() {
        let rest: Vec<usize> = rows.iter().enumerate().filter(|&(k, _)| k != 0 && k != t).map(|(_, &r)| r).collect();
        total = addm(total, mulm(a.get(rows[0], rows[t]), hafnian_expansion(a, &rest)));
    }
    total
}

#[test]
fn pfaffian_hafnian_sparse_and_modular_determinants() {
    let mut st = Stress::new("modmatrix/special");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(200_000 + case as u64);
        let n = 2 * rng.below(5);
        let mut skew = Mat::zeros(n, n);
        let mut sym = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..i {
                let v = if rng.chance(1, 4) { 0 } else { rng.below(M as usize) as u32 };
                skew.set(i, j, subm(0, v));
                skew.set(j, i, v);
                sym.set(i, j, v);
                sym.set(j, i, v);
            }
        }
        let rows: Vec<usize> = (0..n).collect();
        let pf = skew.pfaffian();
        assert_eq!(pf, pfaffian_expansion(&skew, &rows), "pfaffian n={n}");
        assert_eq!(mulm(pf, pf), skew.det());
        assert_eq!(hafnian(&sym), hafnian_expansion(&sym, &rows), "hafnian n={n}");
        st.case();
    }
    // larger hafnians against the recursion over matchings (n = 12)
    for seed in 0..5u64 {
        let mut rng = Rng::new(seed);
        let n = 12;
        let mut sym = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..i {
                let v = rng.below(M as usize) as u32;
                sym.set(i, j, v);
                sym.set(j, i, v);
            }
        }
        assert_eq!(hafnian(&sym), hafnian_expansion(&sym, &(0..n).collect::<Vec<_>>()));
    }
    assert_eq!(hafnian(&Mat::zeros(0, 0)), 1);
    // sparse determinant against the dense one, including singular matrices
    for case in 0..400 * scale() {
        let mut rng = Rng::new(300_000 + case as u64);
        let n = 1 + rng.below(25);
        let k = 1 + rng.below(3 * n);
        let mut cells = std::collections::BTreeMap::new();
        for _ in 0..k {
            cells.insert((rng.below(n), rng.below(n)), 1 + rng.below(M as usize - 1) as u32);
        }
        if case % 4 == 0 {
            for i in 0..n {
                cells.insert((i, i), 1 + rng.below(M as usize - 1) as u32);
            }
        }
        let entries: Vec<(usize, usize, u32)> = cells.iter().map(|(&(i, j), &c)| (i, j, c)).collect();
        let mut dense = Mat::zeros(n, n);
        for &(i, j, c) in &entries {
            dense.set(i, j, c);
        }
        assert_eq!(sparse_determinant::<M>(n, &entries), dense.det(), "sparse n={n}");
        st.case();
    }
    // determinant modulo composite and prime moduli against Laplace over the integers
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(400_000 + case as u64);
        let n = rng.below(6);
        let m = [2u32, 3, 4, 6, 8, 9, 12, 100, 1000, 65_536, 999_999_999, 1_000_000_000, 4_294_967_295, 998_244_353][rng.below(14)];
        let a: Vec<u32> = (0..n * n).map(|_| rng.below(m as usize) as u32).collect();
        let want = {
            fn lap(a: &[i128], n: usize, m: i128) -> i128 {
                if n == 0 {
                    return 1;
                }
                let mut t = 0;
                for j in 0..n {
                    let minor: Vec<i128> = (1..n).flat_map(|i| (0..n).filter(move |&c| c != j).map(move |c| (i, c))).map(|(i, c)| a[i * n + c]).collect();
                    let d = a[j] * lap(&minor, n - 1, m) % m;
                    t += if j % 2 == 0 { d } else { -d };
                }
                t.rem_euclid(m)
            }
            let ai: Vec<i128> = a.iter().map(|&x| x as i128).collect();
            lap(&ai, n, m as i128).rem_euclid(m as i128) as u32
        };
        assert_eq!(determinant_mod(n, &a, m), want, "n={n} m={m} a={a:?}");
        st.case();
    }
    assert_eq!(powm(2, 10), 1024);
    st.done();
}
