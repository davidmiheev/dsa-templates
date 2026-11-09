//! Bell, partition and subset-sum counts, derangements, distinct subsequences.

use hash_map::FastMap;
use modint::StaticModInt;
use polynomial::fps;

type Z<const M: u32> = StaticModInt<M>;

fn inverse_factorials<const M: u32>(n: usize) -> (Vec<Z<M>>, Vec<Z<M>>) {
    polynomial::special::factorials::<M>(n)
}

/// The Bell numbers `B_0, ..., B_n` (`B_k = k! [x^k] exp(e^x - 1)`), NTT prime `M`.
pub fn bell_numbers<const M: u32>(n: usize) -> Vec<Z<M>> {
    let (fact, inv_fact) = inverse_factorials::<M>(n + 1);
    let mut f: Vec<Z<M>> = inv_fact.clone();
    f[0] = Z::zero();
    let g = fps::exp(&f, n + 1);
    (0..=n).map(|k| g[k] * fact[k]).collect()
}

/// The partition numbers `p(0), ..., p(n)` from Euler's pentagonal number theorem
/// (`1 / prod (1 - x^k)`), NTT prime `M`.
pub fn partition_numbers<const M: u32>(n: usize) -> Vec<Z<M>> {
    let mut denom = vec![Z::<M>::zero(); n + 1];
    let mut j = 0i64;
    loop {
        let mut any = false;
        for &g in &[j * (3 * j - 1) / 2, j * (3 * j + 1) / 2] {
            if g <= n as i64 {
                any = true;
                denom[g as usize] = if j % 2 == 0 { Z::one() } else { -Z::one() };
            }
        }
        if !any && j > 0 {
            break;
        }
        j += 1;
    }
    fps::inv(&denom, n + 1)
}

/// `p_t` for `t = 1..=target`: the number of subsets of `s` (a multiset of positive integers)
/// with sum `t`, i.e. the coefficients of `prod (1 + x^{s_i})`, via `exp(sum_s c_s log(1 + x^s))`.
pub fn subset_sum_counts<const M: u32>(s: &[usize], target: usize) -> Vec<Z<M>> {
    let mut count = vec![0u64; target + 1];
    for &x in s {
        assert!(x >= 1 && x <= target, "elements must be in 1..=target");
        count[x] += 1;
    }
    let inv = fps::inverses::<M>(target + 1);
    let mut log = vec![Z::<M>::zero(); target + 1];
    for x in 1..=target {
        if count[x] == 0 {
            continue;
        }
        let c = Z::<M>::new(count[x]);
        let mut j = 1;
        while x * j <= target {
            let term = c * inv[j];
            log[x * j] += if j % 2 == 1 { term } else { -term };
            j += 1;
        }
    }
    let g = fps::exp(&log, target + 1);
    g[1..].to_vec()
}

/// The derangement (Montmort) numbers `a_1, ..., a_n` modulo any `m >= 1`:
/// `a_k = k a_{k-1} + (-1)^k`.
pub fn derangements_mod(n: usize, m: u64) -> Vec<u64> {
    assert!(m >= 1, "the modulus must be positive");
    let mut out = Vec::with_capacity(n);
    let mut a = 1 % m;
    for k in 1..=n as u64 {
        a = (k % m * a + if k % 2 == 0 { 1 } else { m - 1 % m }) % m;
        out.push(a);
    }
    out
}

/// The number of distinct non-empty subsequences of `a`, modulo the prime `M`.
pub fn distinct_subsequences<const M: u32>(a: &[u32]) -> u64 {
    let m = M as u64;
    // dp[i]: distinct subsequences (including the empty one) of the first i elements
    let mut dp = vec![0u64; a.len() + 1];
    dp[0] = 1;
    let mut last: FastMap<usize> = FastMap::with_capacity(a.len());
    for (i, &x) in a.iter().enumerate() {
        let mut v = 2 * dp[i] % m;
        if last.contains(x as u64) {
            let j = *last.get(x as u64);
            v = (v + m - dp[j]) % m;
        }
        dp[i + 1] = v;
        last.insert(x as u64, i);
    }
    (dp[a.len()] + m - 1) % m
}
