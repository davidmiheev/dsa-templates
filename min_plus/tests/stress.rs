//! Stress test: structured min-plus convolutions against the double loop.

use min_plus::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn naive(a: &[i64], b: &[i64]) -> Vec<i64> {
    let mut c = vec![i64::MAX; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            c[i + j] = c[i + j].min(x + y);
        }
    }
    c
}

/// A convex sequence from sorted slopes (ties and negative slopes included).
fn convex(rng: &mut Rng, len: usize) -> Vec<i64> {
    let mut slopes: Vec<i64> = (0..len - 1).map(|_| rng.range(-8, 8)).collect();
    slopes.sort();
    let mut v = vec![rng.range(-20, 20)];
    for s in slopes {
        let last = *v.last().unwrap();
        v.push(last + s);
    }
    v
}

fn concave(rng: &mut Rng, len: usize) -> Vec<i64> {
    convex(rng, len).into_iter().map(|x| -x).collect()
}

#[test]
fn matches_double_loop() {
    let mut st = Stress::new("min_plus");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(case as u64);
        let big = case % 10 == 0;
        let hi = if big { 80 } else { 9 };
        let (n, m) = (rng.range(1, hi) as usize, rng.range(1, hi) as usize);
        let b: Vec<i64> = (0..m).map(|_| rng.range(-30, 30)).collect();
        let ca = convex(&mut rng, n);
        assert_eq!(min_plus_convex_arbitrary(&ca, &b), naive(&ca, &b), "convex/arbitrary {ca:?} {b:?}");
        let cb = convex(&mut rng, m);
        assert_eq!(min_plus_convex_convex(&ca, &cb), naive(&ca, &cb), "convex/convex {ca:?} {cb:?}");
        let co = concave(&mut rng, n);
        assert_eq!(min_plus_concave_arbitrary(&co, &b), naive(&co, &b), "concave/arbitrary {co:?} {b:?}");
        st.case();
        st.ops(n + m);
    }
    st.done();
}
