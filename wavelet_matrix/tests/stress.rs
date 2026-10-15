//! Stress test: `WaveletMatrix` / `WeightedWaveletMatrix` against naive scans of the range.

use wavelet_matrix::{WaveletMatrix, WeightedWaveletMatrix};
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn range(rng: &mut Rng, n: usize) -> (usize, usize) {
    let l = rng.below(n + 1);
    (l, l + rng.below(n + 1 - l)) // empty ranges included
}

#[test]
fn matches_naive_scans() {
    let mut st = Stress::new("wavelet/naive");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(0, 70) as usize;
        let spread = *[1i64, 3, 10, 1_000_000_000].get(case % 4).unwrap(); // many/few distinct values
        let a: Vec<i64> = (0..n).map(|_| rng.range(-spread, spread)).collect();
        let w = WaveletMatrix::new(&a);
        st.case();
        assert_eq!((w.len(), w.is_empty()), (n, n == 0));
        for _ in 0..60 {
            let (l, r) = range(&mut rng, n);
            let x = if rng.chance(1, 2) && n > 0 { a[rng.below(n)] } else { rng.range(-spread - 1, spread + 1) };
            let ctx = format!("[{l},{r}) x={x} case {case} n={n}");
            let slice = &a[l..r];
            assert_eq!(w.count_less(l, r, x), slice.iter().filter(|&&v| v < x).count(), "count_less {ctx}");
            assert_eq!(w.count_less_eq(l, r, x), slice.iter().filter(|&&v| v <= x).count(), "count_less_eq {ctx}");
            assert_eq!(w.count_equal(l, r, x), slice.iter().filter(|&&v| v == x).count(), "count_equal {ctx}");
            let y = rng.range(-spread - 1, spread + 1);
            assert_eq!(w.range_freq(l, r, x, y), slice.iter().filter(|&&v| x <= v && v < y).count(), "range_freq [{x},{y}) {ctx}");
            if l < r {
                let k = rng.below(r - l);
                let mut sorted = slice.to_vec();
                sorted.sort_unstable();
                assert_eq!(w.kth_smallest(l, r, k), sorted[k], "kth_smallest k={k} {ctx}");
            }
            st.ops(5);
        }
    }
    st.done();
}

#[test]
fn weighted_sums_match_naive() {
    let mut st = Stress::new("wavelet/weighted");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(50_000 + case as u64);
        let n = rng.range(0, 70) as usize;
        let spread = *[2i64, 20, 1_000_000_000].get(case % 3).unwrap();
        let a: Vec<i64> = (0..n).map(|_| rng.range(0, spread)).collect();
        let weight: Vec<i64> = (0..n).map(|_| rng.range(-1000, 1000)).collect(); // independent of the value
        let w = WeightedWaveletMatrix::new(&a, &weight);
        st.case();
        for _ in 0..60 {
            let (l, r) = range(&mut rng, n);
            let x = rng.range(-1, spread + 1);
            let want = |keep: &dyn Fn(i64) -> bool| {
                (l..r).filter(|&i| keep(a[i])).fold((0usize, 0i64), |(c, s), i| (c + 1, s + weight[i]))
            };
            assert_eq!(w.count_sum_less_eq(l, r, x), want(&|v| v <= x), "<= {x} [{l},{r}) case {case}");
            assert_eq!(w.count_sum_less(l, r, x), want(&|v| v < x), "< {x} [{l},{r}) case {case}");
            assert_eq!(w.matrix().count_less_eq(l, r, x), want(&|v| v <= x).0);
            st.ops(3);
        }
    }
    st.done();
}

#[test]
fn works_for_other_ordered_types() {
    let words = ["pear", "apple", "fig", "apple", "kiwi"];
    let w = WaveletMatrix::new(&words);
    assert_eq!((w.kth_smallest(0, 5, 0), w.kth_smallest(0, 5, 4)), ("apple", "pear"));
    assert_eq!((w.count_equal(0, 5, "apple"), w.count_less(0, 5, "g")), (2, 3));
}

#[test]
fn invalid_arguments_panic() {
    let w = WaveletMatrix::new(&[3, 1, 2]);
    assert!(catch_unwind(AssertUnwindSafe(|| w.kth_smallest(1, 1, 0))).is_err()); // empty range
    assert!(catch_unwind(AssertUnwindSafe(|| w.kth_smallest(0, 3, 3))).is_err()); // k too large
    assert!(catch_unwind(AssertUnwindSafe(|| w.kth_smallest(0, 4, 0))).is_err()); // r > n
    assert!(catch_unwind(AssertUnwindSafe(|| w.count_less(2, 1, 0))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| WeightedWaveletMatrix::new(&[1, 2], &[1i64]))).is_err());
}

#[test]
fn large_random_ranges() {
    let mut st = Stress::new("wavelet/large");
    let n = 200_000;
    let mut rng = Rng::new(9);
    let a: Vec<u32> = (0..n).map(|_| rng.below(1_000_000_000) as u32).collect();
    let w = WaveletMatrix::new(&a);
    st.case();
    for _ in 0..2_000 * scale() {
        let l = rng.below(n);
        let r = (l + 1 + rng.below(3000)).min(n); // moderate ranges keep the naive scan cheap
        let k = rng.below(r - l);
        let mut sorted = a[l..r].to_vec();
        sorted.sort_unstable();
        assert_eq!(w.kth_smallest(l, r, k), sorted[k]);
        let x = rng.below(1_000_000_000) as u32;
        assert_eq!(w.count_less(l, r, x), sorted.partition_point(|&v| v < x));
        st.ops(2);
    }
    // whole-array order statistics
    let mut all = a.clone();
    all.sort_unstable();
    for k in [0, 1, n / 2, n - 1] {
        assert_eq!(w.kth_smallest(0, n, k), all[k]);
    }
    st.done();
}
