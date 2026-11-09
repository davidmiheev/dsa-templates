//! Stress test: StaticModInt against plain i128 arithmetic.

use modint::StaticModInt;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn check<const M: u32>(st: &mut Stress, rng: &mut Rng, rounds: usize) {
    let m = M as i128;
    for _ in 0..rounds {
        let (a, b) = (rng.next_u64() as i64 as i128, rng.next_u64() as i64 as i128);
        let (x, y) = (StaticModInt::<M>::from(a as i64), StaticModInt::<M>::from(b as i64));
        assert_eq!((x + y).val() as i128, (a + b).rem_euclid(m));
        assert_eq!((x - y).val() as i128, (a - b).rem_euclid(m));
        assert_eq!((x * y).val() as i128, (a * b).rem_euclid(m));
        assert_eq!((-x).val() as i128, (-a).rem_euclid(m));
        if y.val() != 0 {
            assert_eq!(((x / y) * y).val(), x.val());
            assert_eq!((y * y.inv()).val(), 1);
        }
        let e = rng.below(40) as u64;
        let mut p = 1i128;
        for _ in 0..e {
            p = p * a.rem_euclid(m) % m;
        }
        assert_eq!(x.pow(e).val() as i128, p);
        st.ops(1);
    }
}

#[test]
fn matches_integer_arithmetic() {
    let mut st = Stress::new("modint");
    let mut rng = Rng::new(7);
    check::<998_244_353>(&mut st, &mut rng, 20000 * scale());
    check::<1_000_000_007>(&mut st, &mut rng, 20000 * scale());
    check::<7>(&mut st, &mut rng, 5000 * scale());
    check::<3_221_225_473>(&mut st, &mut rng, 5000 * scale()); // > 2^31, prime
    st.case();
    st.done();
}

#[test]
#[should_panic]
fn zero_has_no_inverse() {
    let _ = StaticModInt::<998_244_353>::new(0).inv();
}
