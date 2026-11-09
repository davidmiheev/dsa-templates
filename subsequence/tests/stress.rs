//! Stress test: LIS length against the quadratic DP; the returned indices form a valid subsequence.

use subsequence::longest_increasing_subsequence;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_quadratic_dp() {
    let mut st = Stress::new("subsequence");
    for case in 0..8000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.below(30);
        let range = 1 + rng.below(12);
        let a: Vec<usize> = (0..n).map(|_| rng.below(range)).collect();
        for strict in [true, false] {
            let mut dp = vec![1usize; n];
            for i in 0..n {
                for j in 0..i {
                    if (strict && a[j] < a[i]) || (!strict && a[j] <= a[i]) {
                        dp[i] = dp[i].max(dp[j] + 1);
                    }
                }
            }
            let best = dp.iter().copied().max().unwrap_or(0);
            let idx = longest_increasing_subsequence(&a, strict);
            assert_eq!(idx.len(), best, "a={a:?} strict={strict}");
            assert!(idx.windows(2).all(|w| w[0] < w[1] && if strict { a[w[0]] < a[w[1]] } else { a[w[0]] <= a[w[1]] }));
        }
        st.case();
    }
    st.done();
}
