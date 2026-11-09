//! Stress test: transforms and convolutions against the definitions.

use set_power_series::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const M: u32 = 998_244_353;

fn naive(a: &[u32], b: &[u32], ok: impl Fn(usize, usize) -> Option<usize>) -> Vec<u32> {
    let mut c = vec![0u64; a.len()];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            if let Some(k) = ok(i, j) {
                c[k] = (c[k] + x as u64 * y as u64) % M as u64;
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

fn mulmod(a: u64, b: u64) -> u64 {
    a * b % M as u64
}

fn inv(a: u64) -> u64 {
    let (mut r, mut b, mut e) = (1u64, a % M as u64, M as u64 - 2);
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, b);
        }
        b = mulmod(b, b);
        e >>= 1;
    }
    r
}

/// `powers[m] = s^m` for `m = 0..=count` by repeated subset convolution.
fn powers(s: &[u32], count: usize) -> Vec<Vec<u32>> {
    let mut one = vec![0u32; s.len()];
    one[0] = 1;
    let mut out = vec![one];
    for _ in 0..count {
        let next = naive(out.last().unwrap(), s, |i, j| if i & j == 0 { Some(i | j) } else { None });
        out.push(next);
    }
    out
}

#[test]
fn exp_log_compose_projection_match_definitions() {
    let mut st = Stress::new("set_power_series_exp_compose");
    for case in 0..250 * scale() {
        let mut rng = Rng::new(7_000 + case as u64);
        let n = rng.range(0, 6) as usize;
        let len = 1usize << n;
        let mut s: Vec<u32> = (0..len).map(|_| rng.below(M as usize) as u32).collect();
        let pw = powers(&s, n + 1);
        // exp of a series without constant term
        let mut t = s.clone();
        t[0] = 0;
        let pt = powers(&t, n);
        let mut want = vec![0u64; len];
        let mut fact_inv = 1u64;
        for (m, p) in pt.iter().enumerate() {
            if m > 0 {
                fact_inv = mulmod(fact_inv, inv(m as u64));
            }
            for (w, &x) in want.iter_mut().zip(p) {
                *w = (*w + mulmod(x as u64, fact_inv)) % M as u64;
            }
        }
        let e = subset_exp::<M>(&t);
        assert_eq!(e, want.iter().map(|&x| x as u32).collect::<Vec<_>>());
        // log inverts exp, and exp inverts log
        assert_eq!(subset_log::<M>(&e), t);
        s[0] = 1;
        assert_eq!(subset_exp::<M>(&subset_log::<M>(&s)), s);
        // composition and power projection with a random polynomial
        let m_count = rng.below(9);
        let f: Vec<u32> = (0..m_count).map(|_| rng.below(M as usize) as u32).collect();
        s[0] = rng.below(M as usize) as u32;
        let pw2 = if s[0] == 1 { pw.clone() } else { powers(&s, m_count.max(n + 1)) };
        let mut want = vec![0u64; len];
        for (i, &c) in f.iter().enumerate() {
            for (w, &x) in want.iter_mut().zip(&pw2[i]) {
                *w = (*w + mulmod(x as u64, c as u64)) % M as u64;
            }
        }
        assert_eq!(subset_compose::<M>(&f, &s), want.iter().map(|&x| x as u32).collect::<Vec<_>>());
        let w: Vec<u32> = (0..len).map(|_| rng.below(M as usize) as u32).collect();
        let cnt = rng.below(2 * n + 6);
        let pw3 = powers(&s, cnt);
        let want: Vec<u32> = (0..cnt)
            .map(|m| (0..len).fold(0u64, |a, i| (a + mulmod(w[i] as u64, pw3[m][i] as u64)) % M as u64) as u32)
            .collect();
        assert_eq!(subset_power_projection::<M>(&s, &w, cnt), want);
        st.case();
        st.ops(len * (n + 1));
    }
    st.done();
}
