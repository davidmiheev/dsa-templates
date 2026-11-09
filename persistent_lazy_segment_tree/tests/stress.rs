//! Stress test: `PersistentLazySegmentTree` against one full `Vec` per version, with
//! non-commutative (affine) actions, point sets and range copies between versions.

use algebra::{MapMonoid, Monoid};
use persistent_lazy_segment_tree::PersistentLazySegmentTree;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const P: u64 = 998_244_353;

/// (sum, length) under range-affine actions `x -> a*x + b`.
enum SumLen {}
impl Monoid for SumLen {
    type S = (u64, u64);
    fn identity() -> (u64, u64) {
        (0, 0)
    }
    fn op(a: &(u64, u64), b: &(u64, u64)) -> (u64, u64) {
        ((a.0 + b.0) % P, a.1 + b.1)
    }
}
enum Affine {}
impl MapMonoid for Affine {
    type M = SumLen;
    type F = (u64, u64);
    fn id_map() -> (u64, u64) {
        (1, 0)
    }
    fn apply(f: &(u64, u64), x: &(u64, u64)) -> (u64, u64) {
        ((f.0 * x.0 + f.1 * x.1) % P, x.1)
    }
    fn compose(f: &(u64, u64), g: &(u64, u64)) -> (u64, u64) {
        (f.0 * g.0 % P, (f.0 * g.1 + f.1) % P)
    }
}

fn range(rng: &mut Rng, n: usize) -> (usize, usize) {
    let l = rng.below(n + 1);
    (l, l + rng.below(n + 1 - l))
}

fn check(t: &PersistentLazySegmentTree<Affine>, v: usize, model: &[u64], rng: &mut Rng, ctx: &str) {
    let (l, r) = range(rng, model.len());
    let want = (model[l..r].iter().sum::<u64>() % P, (r - l) as u64);
    assert_eq!(t.prod(v, l, r), want, "prod(v{v}, {l}, {r}) {ctx}");
    let i = rng.below(model.len());
    assert_eq!(t.get(v, i), (model[i], 1), "get(v{v}, {i}) {ctx}");
}

#[test]
fn matches_vec_per_version() {
    let mut st = Stress::new("persistent_lazy_segment_tree/naive");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 40) as usize;
        let init: Vec<u64> = (0..n).map(|_| rng.below(P as usize) as u64).collect();
        let mut t = PersistentLazySegmentTree::<Affine>::from_slice(&init.iter().map(|&x| (x, 1)).collect::<Vec<_>>());
        let mut model = vec![init];
        st.case();
        for step in 0..100 {
            let v = rng.below(model.len());
            let ctx = format!("case {case} step {step} on version {v}");
            let mut m = model[v].clone();
            let nv = match rng.below(5) {
                0 | 1 => {
                    let (l, r) = range(&mut rng, n);
                    let (a, b) = (rng.below(P as usize) as u64, rng.below(P as usize) as u64);
                    for x in &mut m[l..r] {
                        *x = (a * *x + b) % P;
                    }
                    t.apply_range(v, l, r, (a, b))
                }
                2 => {
                    let (i, x) = (rng.below(n), rng.below(P as usize) as u64);
                    m[i] = x;
                    t.set(v, i, (x, 1))
                }
                _ => {
                    let s = rng.below(model.len());
                    let (l, r) = range(&mut rng, n);
                    m[l..r].copy_from_slice(&model[s][l..r]);
                    t.copy_range(v, s, l, r)
                }
            };
            assert_eq!(nv, model.len(), "ids are consecutive {ctx}");
            model.push(m);
            // The new version, a random old one (updates must not leak into it), and the source.
            check(&t, nv, &model[nv], &mut rng, &ctx);
            let w = rng.below(model.len());
            check(&t, w, &model[w], &mut rng, &ctx);
            check(&t, v, &model[v], &mut rng, &ctx);
            st.ops(1);
        }
        for (w, m) in model.iter().enumerate() {
            for i in 0..n {
                assert_eq!(t.get(w, i), (m[i], 1), "final get(v{w}, {i}) case {case}");
            }
        }
    }
    st.done();
}

#[test]
fn panics_on_bad_arguments() {
    let mut t = PersistentLazySegmentTree::<Affine>::from_slice(&[(1, 1), (2, 1), (3, 1)]);
    assert!(catch_unwind(AssertUnwindSafe(|| t.prod(7, 0, 1))).is_err(), "unknown version");
    assert!(catch_unwind(AssertUnwindSafe(|| t.prod(0, 2, 4))).is_err(), "range end");
    assert!(catch_unwind(AssertUnwindSafe(|| t.prod(0, 2, 1))).is_err(), "reversed range");
    assert!(catch_unwind(AssertUnwindSafe(|| t.get(0, 3))).is_err(), "index");
    assert!(catch_unwind(AssertUnwindSafe(|| t.apply_range(0, 0, 4, (1, 1)))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.set(0, 3, (1, 1)))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| t.copy_range(0, 1, 0, 1))).is_err(), "unknown source");
    assert!(catch_unwind(AssertUnwindSafe(|| PersistentLazySegmentTree::<Affine>::from_slice(&[]))).is_err());
    // Empty ranges are fine and just alias the version's content.
    let v = t.apply_range(0, 1, 1, (5, 5));
    assert_eq!(t.prod(v, 0, 3), (6, 3));
    assert_eq!(t.prod(v, 2, 2), (0, 0));
}

#[test]
fn large_scale_memory_is_log_per_update() {
    let mut st = Stress::new("persistent_lazy_segment_tree/large");
    let n = 100_000;
    let q = 100_000 * scale();
    let mut rng = Rng::new(1);
    let init: Vec<u64> = (0..n).map(|i| i as u64 + 1).collect();
    let mut t = PersistentLazySegmentTree::<Affine>::from_slice(&init.iter().map(|&x| (x, 1)).collect::<Vec<_>>());
    st.case();
    for _ in 0..q {
        let v = rng.below(t.num_versions());
        let (l, r) = range(&mut rng, n);
        match rng.below(3) {
            0 => {
                t.apply_range(v, l, r, (rng.below(P as usize) as u64, rng.below(P as usize) as u64));
            }
            1 => {
                let s = rng.below(t.num_versions());
                t.copy_range(v, s, l, r);
            }
            _ => {
                let _ = t.prod(v, l, r);
            }
        }
        st.ops(1);
    }
    // Each update copies a bounded number of nodes (a few times 2*log2(n)), so the total node
    // count stays far below q * n.
    assert!(t.node_count() < 2 * n + q * 160, "{} nodes", t.node_count());
    st.done();
}
