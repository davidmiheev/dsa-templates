//! Stress test: `FenwickTree` against plain prefix sums and a sqrt-decomposition model.

use fenwick::FenwickTree;
use std::panic::catch_unwind;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_naive_prefix_sums() {
    let mut st = Stress::new("fenwick/naive");
    for case in 0..300 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 70) as usize;
        let mut model: Vec<i64> = (0..n).map(|_| rng.range(-1000, 1000)).collect();
        let mut ft = FenwickTree::new(&model);
        st.case();
        for step in 0..150 {
            match rng.below(3) {
                0 => {
                    let (i, d) = (rng.below(n), rng.range(-1000, 1000));
                    model[i] += d;
                    ft.update(i, d);
                }
                1 => {
                    let i = rng.below(n);
                    assert_eq!(ft.pref(i), model[..=i].iter().sum::<i64>(), "pref({i}) case {case} step {step}");
                }
                _ => {
                    let l = rng.below(n);
                    let r = l + rng.below(n - l);
                    assert_eq!(
                        ft.sum_range(l, r),
                        model[l..=r].iter().sum::<i64>(),
                        "sum_range({l},{r}) case {case} step {step}"
                    );
                }
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn works_for_other_numeric_types() {
    let mut st = Stress::new("fenwick/u64+f64");
    for case in 0..100 * scale() {
        let mut rng = Rng::new(500 + case as u64);
        let n = rng.range(1, 50) as usize;
        let mut mu: Vec<u64> = (0..n).map(|_| rng.range(0, 1000) as u64).collect();
        let mut mf: Vec<f64> = mu.iter().map(|&x| x as f64).collect();
        let mut fu = FenwickTree::new(&mu);
        let mut ff = FenwickTree::new(&mf);
        st.case();
        for _ in 0..100 {
            let i = rng.below(n);
            if rng.chance(1, 2) {
                let d = rng.range(0, 1000) as u64;
                mu[i] += d;
                mf[i] += d as f64;
                fu.update(i, d);
                ff.update(i, d as f64);
            } else {
                let l = rng.below(n);
                let r = l + rng.below(n - l);
                assert_eq!(fu.sum_range(l, r), mu[l..=r].iter().sum::<u64>());
                assert_eq!(ff.sum_range(l, r), mf[l..=r].iter().sum::<f64>());
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn out_of_range_arguments_panic() {
    let ft = FenwickTree::new(&[1i64, 2, 3]);
    assert!(catch_unwind(|| ft.pref(3)).is_err());
    assert!(catch_unwind(|| ft.sum_range(2, 1)).is_err());
    assert!(catch_unwind(|| ft.sum_range(0, 3)).is_err());
    let mut ft2 = FenwickTree::new(&[1i64, 2, 3]);
    assert!(catch_unwind(move || ft2.update(3, 1)).is_err());
    let empty: FenwickTree<i64> = FenwickTree::new(&[]);
    assert!(catch_unwind(|| empty.pref(0)).is_err());
}

/// Independent model at scale: blocks of 512 with cached sums.
struct Blocks {
    v: Vec<i64>,
    sums: Vec<i64>,
}

impl Blocks {
    fn prefix(&self, i: usize) -> i64 {
        let b = i / 512;
        self.sums[..b].iter().sum::<i64>() + self.v[b * 512..=i].iter().sum::<i64>()
    }
}

#[test]
fn large_matches_block_model() {
    let mut st = Stress::new("fenwick/large");
    let n = 200_000;
    let mut rng = Rng::new(77);
    let init: Vec<i64> = (0..n).map(|_| rng.range(-1_000_000, 1_000_000)).collect();
    let mut ft = FenwickTree::new(&init);
    let mut m = Blocks { sums: init.chunks(512).map(|c| c.iter().sum()).collect(), v: init };
    st.case();
    for _ in 0..200_000 * scale() {
        let i = rng.below(n);
        if rng.chance(1, 2) {
            let d = rng.range(-1_000_000, 1_000_000);
            m.v[i] += d;
            m.sums[i / 512] += d;
            ft.update(i, d);
        } else {
            assert_eq!(ft.pref(i), m.prefix(i), "pref({i})");
        }
        st.ops(1);
    }
    st.done();
}
