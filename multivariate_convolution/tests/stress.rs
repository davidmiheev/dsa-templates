//! Stress test: both multivariate products against the definition on random shapes.

use multivariate_convolution::*;
use number_theory::is_prime;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn digits(mut x: usize, shape: &[usize]) -> Vec<usize> {
    shape.iter().map(|&n| { let d = x % n; x /= n; d }).collect()
}

fn index(d: &[usize], shape: &[usize]) -> usize {
    let (mut x, mut base) = (0, 1);
    for (i, &n) in shape.iter().enumerate() {
        x += d[i] * base;
        base *= n;
    }
    x
}

fn naive(shape: &[usize], f: &[u32], g: &[u32], p: u64, cyclic: bool) -> Vec<u32> {
    let total: usize = shape.iter().product();
    let mut c = vec![0u64; total];
    for i in 0..total {
        let di = digits(i, shape);
        for j in 0..total {
            let dj = digits(j, shape);
            let mut dk = vec![0; shape.len()];
            let mut ok = true;
            for t in 0..shape.len() {
                let s = di[t] + dj[t];
                if s >= shape[t] {
                    if cyclic { dk[t] = s - shape[t]; } else { ok = false; break; }
                } else {
                    dk[t] = s;
                }
            }
            if ok {
                let k = index(&dk, shape);
                c[k] = (c[k] + f[i] as u64 * g[j] as u64) % p;
            }
        }
    }
    c.into_iter().map(|x| x as u32).collect()
}

#[test]
fn truncated_product_matches_definition() {
    let mut st = Stress::new("multivariate/truncated");
    const P: u64 = 998_244_353;
    for case in 0..300 * scale() {
        let mut rng = Rng::new(case as u64);
        let k = rng.range(0, 4) as usize;
        let shape: Vec<usize> = (0..k).map(|_| rng.range(2, 5) as usize).collect();
        let total: usize = shape.iter().product();
        let f: Vec<u32> = (0..total).map(|_| rng.below(P as usize) as u32).collect();
        let g: Vec<u32> = (0..total).map(|_| rng.below(P as usize) as u32).collect();
        assert_eq!(multivariate_convolution(&shape, &f, &g), naive(&shape, &f, &g, P, false), "shape {shape:?}");
        st.case();
    }
    // many binary variables (the subset-convolution regime)
    let shape = vec![2usize; 9];
    let mut rng = Rng::new(5);
    let f: Vec<u32> = (0..512).map(|_| rng.below(P as usize) as u32).collect();
    let g: Vec<u32> = (0..512).map(|_| rng.below(P as usize) as u32).collect();
    assert_eq!(multivariate_convolution(&shape, &f, &g), naive(&shape, &f, &g, P, false));
    st.done();
}

#[test]
fn cyclic_product_matches_definition() {
    let mut st = Stress::new("multivariate/cyclic");
    // primes p and divisors of p - 1 as axis lengths, including axes longer than 64
    let cases: Vec<(u64, Vec<usize>)> = vec![
        (3, vec![2]), (3, vec![2, 2]), (7, vec![2, 3]), (7, vec![6]), (13, vec![3, 4, 2]),
        (97, vec![96]), (97, vec![8, 12]), (257, vec![256]), (1_000_000_007, vec![2, 2]),
        (998_244_353, vec![4, 7, 17]), (65_537, vec![128]), (65_537, vec![256, 4]), (12_289, vec![192]),
    ];
    for (p, shape) in cases {
        assert!(is_prime(p));
        let total: usize = shape.iter().product();
        let mut rng = Rng::new(p);
        for _ in 0..2 * scale() {
            let f: Vec<u32> = (0..total).map(|_| rng.below(p as usize) as u32).collect();
            let g: Vec<u32> = (0..total).map(|_| rng.below(p as usize) as u32).collect();
            assert_eq!(multivariate_convolution_cyclic(p, &shape, &f, &g), naive(&shape, &f, &g, p, true), "p={p} shape={shape:?}");
            st.case();
        }
    }
    assert_eq!(multivariate_convolution_cyclic(7, &[], &[3], &[4]), vec![5]);
    st.done();
}
