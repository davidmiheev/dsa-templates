//! Stress test: determinant against Leibniz expansion (all permutations).

use modint::ModInt998244353 as Mint;
use modmatrix::determinant;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn leibniz(a: &[Vec<Mint>]) -> Mint {
    let n = a.len();
    let mut perm: Vec<usize> = (0..n).collect();
    let mut total = Mint::zero();
    fn rec(k: usize, perm: &mut Vec<usize>, a: &[Vec<Mint>], total: &mut Mint) {
        if k == perm.len() {
            let mut sign = 1i32;
            let mut seen = vec![false; perm.len()];
            for s in 0..perm.len() {
                if !seen[s] {
                    let (mut x, mut len) = (s, 0);
                    while !seen[x] {
                        seen[x] = true;
                        x = perm[x];
                        len += 1;
                    }
                    if len % 2 == 0 {
                        sign = -sign;
                    }
                }
            }
            let mut prod = Mint::one();
            for (i, &p) in perm.iter().enumerate() {
                prod *= a[i][p];
            }
            *total += if sign > 0 { prod } else { -prod };
            return;
        }
        for i in k..perm.len() {
            perm.swap(k, i);
            rec(k + 1, perm, a, total);
            perm.swap(k, i);
        }
    }
    rec(0, &mut perm, a, &mut total);
    total
}

#[test]
fn matches_leibniz() {
    let mut st = Stress::new("modmatrix/det");
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(0, 6) as usize;
        let small = rng.chance(1, 2);
        let a: Vec<Vec<Mint>> = (0..n)
            .map(|_| (0..n).map(|_| Mint::from(if small { rng.range(-2, 2) } else { rng.next_u64() as i64 })).collect())
            .collect();
        assert_eq!(determinant(a.clone()), leibniz(&a), "case {case}: {a:?}");
        st.case();
    }
    st.done();
}
