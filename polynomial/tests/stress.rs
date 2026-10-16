//! Stress tests for the series core: every operation against its defining identity or a naive
//! O(n^2) computation.

use modint::ModInt998244353 as Mint;
use polynomial::fps;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn random(rng: &mut Rng, n: usize) -> Vec<Mint> {
    (0..n).map(|_| Mint::new(rng.next_u64())).collect()
}

fn naive_mul(a: &[Mint], b: &[Mint]) -> Vec<Mint> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut c = vec![Mint::zero(); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            c[i + j] += a[i] * b[j];
        }
    }
    c
}

#[test]
fn series_identities() {
    let mut st = Stress::new("polynomial/fps");
    for case in 0..300 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 70) as usize;
        let len = rng.range(1, 80) as usize;
        let mut f = random(&mut rng, len);
        // inverse
        if f[0].val() == 0 {
            f[0] = Mint::one();
        }
        let g = fps::inv(&f, n);
        let mut prod = naive_mul(&f, &g);
        prod.resize(n.max(prod.len()), Mint::zero());
        assert_eq!(prod[0], Mint::one());
        assert!(prod[1..n].iter().all(|x| x.val() == 0), "inv case {case}");
        // log / exp
        let mut h = f.clone();
        h[0] = Mint::one();
        let lg = fps::log(&h, n);
        let mut z = lg.clone();
        z[0] = Mint::zero();
        let back = fps::exp(&lg, n);
        assert_eq!(back, fps::resized(&h, n), "exp(log) case {case}");
        // derivative of log h = h' / h
        let d = fps::derivative(&lg);
        let want = fps::mul_trunc(&fps::derivative(&h), &fps::inv(&h, n), n - 1);
        assert_eq!(d, want[..n.saturating_sub(1)].to_vec());
        // pow against repeated multiplication
        let e = rng.range(0, 9) as u64;
        let mut p = vec![Mint::one()];
        for _ in 0..e {
            p = fps::mul_trunc(&p, &f, n);
        }
        p.resize(n, Mint::zero());
        assert_eq!(fps::pow(&f, e, n), p, "pow case {case} e={e}");
        // pow with leading zeros
        let shift = rng.range(0, 4) as usize;
        let mut fz = vec![Mint::zero(); shift];
        fz.extend_from_slice(&f);
        let mut pz = vec![Mint::one()];
        for _ in 0..e {
            pz = fps::mul_trunc(&pz, &fz, n);
        }
        pz.resize(n, Mint::zero());
        assert_eq!(fps::pow(&fz, e, n), pz, "pow with zeros case {case}");
        st.case();
        st.ops(n);
    }
    st.done();
}

#[test]
fn square_roots() {
    let mut rng = Rng::new(99);
    for case in 0..200 * scale() {
        let n = rng.range(1, 60) as usize;
        let shift = 2 * rng.range(0, 3) as usize;
        let mut g = vec![Mint::zero(); shift];
        let glen = rng.range(1, 40) as usize;
        g.extend(random(&mut rng, glen));
        let f = naive_mul(&g, &g);
        let r = fps::sqrt(&f, n);
        assert_eq!(r.len(), n, "a square must have a root (case {case})");
        let mut sq = naive_mul(&r, &r);
        sq.resize(n.max(sq.len()), Mint::zero());
        assert_eq!(&sq[..n], &fps::resized(&f, n)[..], "case {case}");
        // an odd leading exponent below n has no root
        let mut odd = vec![Mint::zero()];
        odd.extend(f.iter().copied());
        let k = odd.iter().position(|x| x.val() != 0).unwrap();
        if k % 2 == 1 && k < n {
            assert!(fps::sqrt(&odd, n).is_empty(), "case {case}");
        }
    }
}

#[test]
fn division_with_remainder() {
    let mut rng = Rng::new(7);
    for case in 0..300 * scale() {
        let (n, m) = (rng.range(1, 60) as usize, rng.range(1, 40) as usize);
        let f = random(&mut rng, n);
        let mut g = random(&mut rng, m);
        if g[m - 1].val() == 0 {
            g[m - 1] = Mint::one();
        }
        let (q, r) = fps::div_rem(&f, &g);
        let mut recomposed = naive_mul(&q, &g);
        recomposed.resize(recomposed.len().max(r.len()).max(f.len()), Mint::zero());
        for (i, x) in r.iter().enumerate() {
            recomposed[i] += *x;
        }
        let mut ft = f.clone();
        ft.resize(recomposed.len(), Mint::zero());
        assert_eq!(recomposed, ft, "case {case}");
        assert!(r.len() < m || n < m);
    }
}

#[test]
fn large_sizes() {
    let mut rng = Rng::new(3);
    let n = 1 << 17;
    let mut f = random(&mut rng, n);
    f[0] = Mint::zero();
    let e = fps::exp(&f, n);
    let mut h = e.clone();
    assert_eq!(h[0], Mint::one());
    let l = fps::log(&h, n);
    assert_eq!(l, f);
    h[0] = Mint::new(5);
    let g = fps::inv(&h, n);
    let p = fps::mul_trunc(&g, &h, n);
    assert_eq!(p[0], Mint::one());
    assert!(p[1..].iter().all(|x| x.val() == 0));
}

#[test]
fn sparse_versions_match_dense() {
    use polynomial::sparse;
    let mut rng = Rng::new(21);
    for case in 0..300 * scale() {
        let n = rng.range(1, 80) as usize;
        let k = rng.range(0, 5) as usize;
        let mut idx: Vec<usize> = (0..k).map(|_| rng.below(n)).collect();
        idx.sort();
        idx.dedup();
        let terms_with = |first: Option<usize>, rng: &mut Rng| -> Vec<(usize, Mint)> {
            let mut t: Vec<(usize, Mint)> = Vec::new();
            if let Some(i) = first {
                t.push((i, Mint::new(1 + rng.next_u64() % 998_244_352)));
            }
            for &i in &idx {
                if first.map_or(true, |f| i > f) {
                    t.push((i, Mint::new(1 + rng.next_u64() % 998_244_352)));
                }
            }
            t
        };
        let dense = |t: &[(usize, Mint)]| -> Vec<Mint> {
            let mut v = vec![Mint::zero(); n];
            for &(i, x) in t {
                v[i] = x;
            }
            v
        };
        let t0 = terms_with(Some(0), &mut rng);
        assert_eq!(sparse::inv_sparse(n, &t0), fps::inv(&dense(&t0), n), "inv case {case}");
        let mut t1 = t0.clone();
        t1[0].1 = Mint::one();
        assert_eq!(sparse::log_sparse(n, &t1), fps::log(&dense(&t1), n), "log case {case}");
        let tz = terms_with(None, &mut rng).into_iter().filter(|&(i, _)| i > 0).collect::<Vec<_>>();
        assert_eq!(sparse::exp_sparse(n, &tz), fps::exp(&dense(&tz), n), "exp case {case}");
        let e = if rng.chance(1, 4) { rng.next_u64() % 1_000_000_000_000_000_000 } else { rng.range(0, 12) as u64 };
        let tp = terms_with(if rng.chance(1, 2) { Some(rng.below(n)) } else { None }, &mut rng);
        assert_eq!(sparse::pow_sparse(n, &tp, e), fps::pow(&dense(&tp), e, n), "pow case {case} e={e}");
        assert_eq!(sparse::sqrt_sparse(n, &tp).is_empty(), fps::sqrt(&dense(&tp), n).is_empty(), "sqrt existence case {case}");
        let r = sparse::sqrt_sparse(n, &tp);
        if !r.is_empty() {
            let sq = fps::mul_trunc(&r, &r, n);
            assert_eq!(sq, dense(&tp), "sqrt case {case}");
        }
    }
}

fn eval(f: &[Mint], x: Mint) -> Mint {
    f.iter().rev().fold(Mint::zero(), |acc, &c| acc * x + c)
}

#[test]
fn special_transforms() {
    use polynomial::special::*;
    let mut rng = Rng::new(55);
    for case in 0..250 * scale() {
        let n = rng.range(1, 40) as usize;
        let f = random(&mut rng, n);
        // Taylor shift
        let c = Mint::new(rng.next_u64());
        let g = taylor_shift(&f, c);
        for _ in 0..3 {
            let x = Mint::new(rng.next_u64());
            assert_eq!(eval(&g, x), eval(&f, x + c), "taylor case {case}");
        }
        // sampling shift, including shifts that wrap around the modulus and hit the sample range
        let samples: Vec<Mint> = (0..n).map(|i| eval(&f, Mint::from(i))).collect();
        let m = rng.range(1, 40) as usize;
        let shifts = [Mint::new(rng.next_u64()), Mint::from(rng.below(n + 3)), Mint::from(-(rng.below(8) as i64)), Mint::new(998_244_353 - 3)];
        for &c in &shifts {
            let got = shift_sampling_points(&samples, c, m);
            for (i, &v) in got.iter().enumerate() {
                assert_eq!(v, eval(&f, c + Mint::from(i)), "shift case {case} c={c:?} i={i} n={n}");
            }
        }
        // geometric evaluation
        let (a, r) = (Mint::new(rng.next_u64()), if rng.chance(1, 8) { Mint::zero() } else { Mint::new(rng.next_u64()) });
        let ev = eval_geometric(&f, a, r, m);
        let mut x = a;
        for (i, &v) in ev.iter().enumerate() {
            assert_eq!(v, eval(&f, x), "geometric case {case} i={i}");
            x *= r;
        }
        // geometric interpolation (random ratio is of huge order, so the points are distinct)
        let q = if n == 2 && rng.chance(1, 3) { Mint::zero() } else { Mint::new(2 + rng.next_u64() % 1000) };
        let a2 = if n >= 2 { Mint::new(1 + rng.next_u64() % 998_244_352) } else { a };
        let ys: Vec<Mint> = {
            let mut x = a2;
            (0..n).map(|_| { let v = eval(&f, x); x *= q; v }).collect()
        };
        assert_eq!(interpolate_geometric(&ys, a2, q), f, "interpolation case {case}");
        // prefix sums
        let pre = prefix_sum_polynomial(&f);
        assert_eq!(pre.len(), n + 1);
        let mut acc = Mint::zero();
        for t in 0..n + 3 {
            assert_eq!(eval(&pre, Mint::from(t)), acc, "prefix case {case} t={t}");
            acc += eval(&f, Mint::from(t));
        }
        // product of a sequence
        let k = rng.range(0, 6) as usize;
        let polys: Vec<Vec<Mint>> = (0..k)
            .map(|_| {
                let len = rng.range(1, 6) as usize;
                random(&mut rng, len)
            })
            .collect();
        let mut want = vec![Mint::one()];
        for p in &polys {
            want = naive_mul(&want, p);
        }
        assert_eq!(product_all(polys), want);
    }
}

#[test]
fn multipoint_and_interpolation() {
    use polynomial::multipoint::*;
    let mut rng = Rng::new(77);
    for case in 0..250 * scale() {
        let n = rng.range(1, 50) as usize;
        let m = rng.range(1, 50) as usize;
        let f = random(&mut rng, n);
        let pts = random(&mut rng, m);
        let got = multipoint_evaluate(&f, &pts);
        for (i, &x) in pts.iter().enumerate() {
            assert_eq!(got[i], eval(&f, x), "eval case {case} n={n} m={m}");
        }
        // interpolation through distinct points recovers the polynomial
        let k = rng.range(1, 50) as usize;
        let mut xs: Vec<Mint> = Vec::new();
        while xs.len() < k {
            let x = Mint::new(rng.next_u64());
            if !xs.contains(&x) {
                xs.push(x);
            }
        }
        let g = random(&mut rng, k);
        let ys = multipoint_evaluate(&g, &xs);
        assert_eq!(interpolate(&xs, &ys), g, "interpolation case {case}");
        // Newton basis: expand and compare
        let b = monomial_to_newton(&g, &xs);
        let mut total = vec![Mint::zero(); k];
        let mut basis = vec![Mint::one()];
        for i in 0..k {
            for (j, &c) in basis.iter().enumerate() {
                total[j] += b[i] * c;
            }
            basis = naive_mul(&basis, &[-xs[i], Mint::one()]);
        }
        assert_eq!(total, g, "newton case {case}");
    }
}

#[test]
fn gcd_and_modular_inverse() {
    use polynomial::gcd::*;
    let mut rng = Rng::new(123);
    for case in 0..400 * scale() {
        // sizes straddle the half-gcd threshold; a shared factor half of the time
        let big = case % 8 == 0;
        let (na, nb) = if big { (rng.range(100, 400) as usize, rng.range(100, 400) as usize) } else { (rng.range(1, 30) as usize, rng.range(1, 30) as usize) };
        let mut a = random(&mut rng, na);
        let mut b = random(&mut rng, nb);
        if rng.chance(1, 2) {
            let cl = rng.range(1, 6) as usize;
            let common = random(&mut rng, cl);
            a = naive_mul(&a, &common);
            b = naive_mul(&b, &common);
        }
        let want = gcd_naive(&a, &b);
        let (g, s, t) = ext_gcd(&a, &b);
        assert_eq!(g, want, "gcd case {case}");
        let mut comb = {
            let x = naive_mul(&s, &a);
            let y = naive_mul(&t, &b);
            let mut r = vec![Mint::zero(); x.len().max(y.len())];
            for (i, &c) in x.iter().enumerate() { r[i] += c; }
            for (i, &c) in y.iter().enumerate() { r[i] += c; }
            r
        };
        fps::trim(&mut comb);
        assert_eq!(comb, g, "bezout case {case}");
        // inverse modulo a polynomial
        let h = inverse_mod(&b, &a);
        if want.len() == 1 || a.len() == 1 {
            if !h.is_empty() && a.len() > 1 {
                let (_, mut r) = fps::div_rem(&naive_mul(&h, &b), &a);
                fps::trim(&mut r);
                assert_eq!(r, vec![Mint::one()], "inverse case {case}");
            }
        } else {
            assert!(h.is_empty(), "inverse must not exist when the gcd is non-trivial (case {case})");
        }
    }
}

#[test]
fn root_finding() {
    use polynomial::roots::find_roots;
    let mut rng = Rng::new(31);
    // a quadratic non-residue gives the irreducible x^2 - nr
    let nr = (2..).map(Mint::from).find(|&z: &Mint| z.pow(998_244_352 / 2) == -Mint::one()).unwrap();
    for case in 0..60 * scale() {
        let k = rng.range(0, 12) as usize;
        let mut want: Vec<Mint> = (0..k).map(|_| Mint::new(rng.next_u64())).collect();
        if rng.chance(1, 4) {
            want.push(Mint::zero());
        }
        let mut f = vec![Mint::new(1 + rng.next_u64() % 998_244_352)];
        for &r in &want {
            f = naive_mul(&f, &[-r, Mint::one()]);
        }
        for _ in 0..rng.range(0, 3) {
            f = naive_mul(&f, &[-nr, Mint::zero(), Mint::one()]);
        }
        if rng.chance(1, 3) && !want.is_empty() {
            f = naive_mul(&f, &[-want[0], Mint::one()]); // a repeated root
        }
        let mut expect = want.clone();
        expect.sort_by_key(|r| r.val());
        expect.dedup();
        assert_eq!(find_roots(&f), expect, "case {case}");
    }
}

#[test]
fn bivariate_inverse() {
    let mut rng = Rng::new(8);
    for case in 0..150 * scale() {
        let (n, m) = (rng.range(1, 12) as usize, rng.range(1, 12) as usize);
        let mut f = random(&mut rng, n * m);
        if f[0].val() == 0 {
            f[0] = Mint::one();
        }
        let g = fps::inv_bivariate(&f, n, m);
        // multiply back naively
        let mut prod = vec![Mint::zero(); n * m];
        for i in 0..n {
            for j in 0..m {
                for k in 0..=i {
                    for l in 0..=j {
                        prod[i * m + j] += f[k * m + l] * g[(i - k) * m + (j - l)];
                    }
                }
            }
        }
        for (idx, &v) in prod.iter().enumerate() {
            assert_eq!(v, if idx == 0 { Mint::one() } else { Mint::zero() }, "case {case} {n}x{m} idx {idx}");
        }
    }
}
