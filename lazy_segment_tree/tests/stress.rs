//! Stress test: `LazySegmentTree` under three different actions against naive arrays, for sizes
//! that are not powers of two, including `max_right` / `min_left`.

use algebra::{MapMonoid, Monoid};
use lazy_segment_tree::LazySegmentTree;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

// ---- 1. range add / range min --------------------------------------------------------------
enum Min {}
impl Monoid for Min {
    type S = i64;
    fn identity() -> i64 {
        i64::MAX
    }
    fn op(a: &i64, b: &i64) -> i64 {
        *a.min(b)
    }
}
enum AddMin {}
impl MapMonoid for AddMin {
    type M = Min;
    type F = i64;
    fn id_map() -> i64 {
        0
    }
    fn apply(f: &i64, x: &i64) -> i64 {
        if *x == i64::MAX { *x } else { x + f }
    }
    fn compose(f: &i64, g: &i64) -> i64 {
        f + g
    }
}

// ---- 2. range affine / range sum (mod p), elements carry their length ------------------------
const P: u64 = 998_244_353;
enum SumLen {}
impl Monoid for SumLen {
    type S = (u64, u64); // (sum, length)
    fn identity() -> (u64, u64) {
        (0, 0)
    }
    fn op(a: &(u64, u64), b: &(u64, u64)) -> (u64, u64) {
        ((a.0 + b.0) % P, a.1 + b.1)
    }
}
enum AffineSum {}
impl MapMonoid for AffineSum {
    type M = SumLen;
    type F = (u64, u64); // x -> a * x + b
    fn id_map() -> (u64, u64) {
        (1, 0)
    }
    fn apply(f: &(u64, u64), x: &(u64, u64)) -> (u64, u64) {
        ((f.0 * x.0 + f.1 * (x.1 % P)) % P, x.1)
    }
    fn compose(f: &(u64, u64), g: &(u64, u64)) -> (u64, u64) {
        (f.0 * g.0 % P, (f.0 * g.1 + f.1) % P) // f after g
    }
}

// ---- 3. range assign / range sum (non-commutative composition) -------------------------------
enum SumLenI {}
impl Monoid for SumLenI {
    type S = (i64, i64);
    fn identity() -> (i64, i64) {
        (0, 0)
    }
    fn op(a: &(i64, i64), b: &(i64, i64)) -> (i64, i64) {
        (a.0 + b.0, a.1 + b.1)
    }
}
enum AssignSum {}
impl MapMonoid for AssignSum {
    type M = SumLenI;
    type F = Option<i64>;
    fn id_map() -> Option<i64> {
        None
    }
    fn apply(f: &Option<i64>, x: &(i64, i64)) -> (i64, i64) {
        match f {
            Some(c) => (c * x.1, x.1),
            None => *x,
        }
    }
    fn compose(f: &Option<i64>, g: &Option<i64>) -> Option<i64> {
        f.or(*g) // the later action wins
    }
}

fn range(rng: &mut Rng, n: usize) -> (usize, usize) {
    let l = rng.below(n + 1);
    (l, l + rng.below(n + 1 - l)) // empty ranges included
}

#[test]
fn add_min_matches_naive() {
    let mut st = Stress::new("lazy_segtree/add_min");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 70) as usize;
        let mut model: Vec<i64> = (0..n).map(|_| rng.range(-50, 50)).collect();
        let mut t = LazySegmentTree::<AddMin>::from_slice(&model);
        st.case();
        for step in 0..120 {
            let ctx = format!("case {case} step {step} n={n}");
            let (l, r) = range(&mut rng, n);
            match rng.below(7) {
                0 | 1 => {
                    let d = rng.range(-20, 20);
                    t.apply_range(l, r, d);
                    model[l..r].iter_mut().for_each(|x| *x += d);
                }
                2 => {
                    let (i, v) = (rng.below(n), rng.range(-50, 50));
                    t.set(i, v);
                    model[i] = v;
                }
                3 => {
                    let (i, d) = (rng.below(n), rng.range(-9, 9));
                    t.apply(i, d);
                    model[i] += d;
                }
                4 => {
                    let i = rng.below(n);
                    assert_eq!(t.get(i), model[i], "get({i}) {ctx}");
                }
                5 => {
                    let thr = rng.range(-60, 60);
                    // largest r' with min(a[l..r']) > thr  (monotone: true on a prefix)
                    let want = (l..=n).take_while(|&e| model[l..e].iter().min().map_or(true, |&m| m > thr)).last().unwrap();
                    assert_eq!(t.max_right(l, |&m| m > thr), want, "max_right({l}, >{thr}) {ctx}");
                    let want = (0..=r).rev().take_while(|&s| model[s..r].iter().min().map_or(true, |&m| m > thr)).last().unwrap();
                    assert_eq!(t.min_left(r, |&m| m > thr), want, "min_left({r}, >{thr}) {ctx}");
                }
                _ => {
                    let want = model[l..r].iter().copied().min().unwrap_or(i64::MAX);
                    assert_eq!(t.prod(l, r), want, "prod({l},{r}) {ctx}");
                }
            }
            assert_eq!(t.all_prod(), model.iter().copied().min().unwrap(), "all_prod {ctx}");
            st.ops(1);
        }
        for i in 0..n {
            assert_eq!(t.get(i), model[i], "final get({i}) case {case}");
        }
        st.ops(n);
    }
    st.done();
}

#[test]
fn affine_sum_matches_naive() {
    let mut st = Stress::new("lazy_segtree/affine_sum");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(10_000 + case as u64);
        let n = rng.range(1, 70) as usize;
        let mut model: Vec<u64> = (0..n).map(|_| rng.below(P as usize) as u64).collect();
        let init: Vec<(u64, u64)> = model.iter().map(|&x| (x, 1)).collect();
        let mut t = LazySegmentTree::<AffineSum>::from_slice(&init);
        st.case();
        for step in 0..120 {
            let (l, r) = range(&mut rng, n);
            let ctx = format!("case {case} step {step} n={n}");
            if rng.chance(1, 2) {
                let (a, b) = (rng.below(P as usize) as u64, rng.below(P as usize) as u64);
                t.apply_range(l, r, (a, b));
                model[l..r].iter_mut().for_each(|x| *x = (a * *x + b) % P);
            } else {
                let sum = model[l..r].iter().fold(0u64, |s, &x| (s + x) % P);
                assert_eq!(t.prod(l, r), (sum, (r - l) as u64), "prod({l},{r}) {ctx}");
            }
            st.ops(1);
        }
        for i in 0..n {
            assert_eq!(t.get(i), (model[i], 1), "final get({i}) case {case}");
        }
        st.ops(n);
    }
    st.done();
}

#[test]
fn assign_sum_matches_naive_including_search() {
    let mut st = Stress::new("lazy_segtree/assign_sum");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(20_000 + case as u64);
        let n = rng.range(1, 70) as usize;
        let mut model: Vec<i64> = (0..n).map(|_| rng.range(0, 9)).collect(); // non-negative: sums monotone
        let init: Vec<(i64, i64)> = model.iter().map(|&x| (x, 1)).collect();
        let mut t = LazySegmentTree::<AssignSum>::from_slice(&init);
        st.case();
        for step in 0..120 {
            let (l, r) = range(&mut rng, n);
            let ctx = format!("case {case} step {step} n={n}");
            match rng.below(4) {
                0 | 1 => {
                    let c = rng.range(0, 9);
                    t.apply_range(l, r, Some(c));
                    model[l..r].iter_mut().for_each(|x| *x = c);
                }
                2 => {
                    let limit = rng.range(0, 60);
                    let want = (l..=n).take_while(|&e| model[l..e].iter().sum::<i64>() <= limit).last().unwrap();
                    assert_eq!(t.max_right(l, |s| s.0 <= limit), want, "max_right({l}, <={limit}) {ctx}");
                    let want = (0..=r).rev().take_while(|&s| model[s..r].iter().sum::<i64>() <= limit).last().unwrap();
                    assert_eq!(t.min_left(r, |s| s.0 <= limit), want, "min_left({r}, <={limit}) {ctx}");
                }
                _ => assert_eq!(t.prod(l, r), (model[l..r].iter().sum(), (r - l) as i64), "prod({l},{r}) {ctx}"),
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn empty_tree_and_invalid_arguments_panic() {
    let mut t = LazySegmentTree::<AddMin>::new(0);
    assert!(t.is_empty());
    assert_eq!((t.len(), t.prod(0, 0), t.all_prod()), (0, i64::MAX, i64::MAX));
    assert_eq!(t.max_right(0, |_| true), 0);
    let mut t = LazySegmentTree::<AddMin>::from_slice(&[1, 2, 3]);
    assert!(catch_unwind(AssertUnwindSafe(|| t.prod(2, 4))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.prod(2, 1))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.get(3))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.set(3, 0))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.apply_range(0, 4, 1))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.max_right(0, |&m| m < 0))).is_err()); // pred(identity) is false
    assert!(catch_unwind(AssertUnwindSafe(|| t.min_left(3, |&m| m < 0))).is_err());
}

#[test]
fn large_random_matches_block_model() {
    // n = 200k, range add / range min against a sqrt-decomposition model.
    let mut st = Stress::new("lazy_segtree/large");
    let n = 200_000;
    let mut rng = Rng::new(3);
    let init: Vec<i64> = (0..n).map(|_| rng.range(-1_000_000, 1_000_000)).collect();
    let mut t = LazySegmentTree::<AddMin>::from_slice(&init);
    const B: usize = 512;
    let mut v = init;
    let mut add = vec![0i64; (n + B - 1) / B];
    let mut bmin: Vec<i64> = v.chunks(B).map(|c| *c.iter().min().unwrap()).collect();
    st.case();
    for _ in 0..30_000 * scale() {
        let l = rng.below(n);
        let r = l + 1 + rng.below(n - l);
        let (bl, br) = (l / B, (r - 1) / B);
        let rebuild = |v: &Vec<i64>, add: &Vec<i64>, bmin: &mut Vec<i64>, b: usize| {
            let (s, e) = (b * B, ((b + 1) * B).min(n));
            bmin[b] = v[s..e].iter().min().unwrap() + add[b];
        };
        if rng.chance(1, 2) {
            let d = rng.range(-1000, 1000);
            t.apply_range(l, r, d);
            if bl == br {
                v[l..r].iter_mut().for_each(|x| *x += d);
                rebuild(&v, &add, &mut bmin, bl);
            } else {
                v[l..(bl + 1) * B].iter_mut().for_each(|x| *x += d);
                v[br * B..r].iter_mut().for_each(|x| *x += d);
                rebuild(&v, &add, &mut bmin, bl);
                rebuild(&v, &add, &mut bmin, br);
                for b in bl + 1..br {
                    add[b] += d;
                    bmin[b] += d;
                }
            }
        } else {
            let want = if bl == br {
                v[l..r].iter().min().unwrap() + add[bl]
            } else {
                let a = v[l..(bl + 1) * B].iter().min().unwrap() + add[bl];
                let b = v[br * B..r].iter().min().unwrap() + add[br];
                (bl + 1..br).map(|k| bmin[k]).chain([a, b]).min().unwrap()
            };
            assert_eq!(t.prod(l, r), want, "prod({l},{r})");
        }
        st.ops(1);
    }
    st.done();
}
