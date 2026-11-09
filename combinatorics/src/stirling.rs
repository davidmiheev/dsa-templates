//! Stirling numbers of both kinds.

use modint::StaticModInt;
use polynomial::{fps, special};

type Z<const M: u32> = StaticModInt<M>;

/// The signed Stirling numbers of the first kind `s(n, 0), ..., s(n, n)`, i.e. the coefficients
/// of `x (x - 1) ... (x - n + 1)`, in `O(n log n)` by doubling (`P_{2m}(x) = P_m(x) P_m(x - m)`).
/// NTT prime `M`.
pub fn stirling1_row<const M: u32>(n: usize) -> Vec<Z<M>> {
    if n == 0 {
        return vec![Z::one()];
    }
    if n % 2 == 1 {
        let p = stirling1_row::<M>(n - 1);
        // times (x - (n - 1))
        let c = Z::<M>::new(n as u64 - 1);
        let mut r = vec![Z::zero(); n + 1];
        for (i, &v) in p.iter().enumerate() {
            r[i + 1] += v;
            r[i] -= c * v;
        }
        return r;
    }
    let half = n / 2;
    let a = stirling1_row::<M>(half);
    let b = special::taylor_shift(&a, -Z::<M>::new(half as u64));
    fps::mul(&a, &b)
}

/// `s(m, k)` for `m = k..=n` (signed first kind), from `log(1 + x)^k / k!`. NTT prime `M`.
pub fn stirling1_fixed_k<const M: u32>(n: usize, k: usize) -> Vec<Z<M>> {
    assert!(k <= n, "k must not exceed n");
    let len = n - k + 1;
    let (fact, inv_fact) = special::factorials::<M>(n + 1);
    let inv = fps::inverses::<M>(len + 1);
    // log(1 + x) / x = sum (-1)^i x^i / (i + 1)
    let t: Vec<Z<M>> = (0..len).map(|i| if i % 2 == 0 { inv[i + 1] } else { -inv[i + 1] }).collect();
    let tk = fps::pow(&t, k as u64, len);
    (0..len).map(|i| tk[i] * fact[k + i] * inv_fact[k]).collect()
}

/// `j^e` for all `j <= n`, by a linear sieve (one modular power per prime).
fn power_table<const M: u32>(n: usize, e: u64) -> Vec<Z<M>> {
    let mut pw = vec![Z::<M>::zero(); n + 1];
    let mut primes: Vec<usize> = Vec::new();
    let mut composite = vec![false; n + 1];
    if n >= 1 {
        pw[1] = Z::one();
    }
    if e == 0 {
        pw[0] = Z::one();
    }
    for i in 2..=n {
        if !composite[i] {
            primes.push(i);
            pw[i] = Z::<M>::new(i as u64).pow(e);
        }
        for &p in &primes {
            if i * p > n {
                break;
            }
            composite[i * p] = true;
            pw[i * p] = pw[i] * pw[p];
            if i % p == 0 {
                break;
            }
        }
    }
    pw
}

/// The Stirling numbers of the second kind `S(n, 0), ..., S(n, n)`:
/// `S(n, k) = sum_j (-1)^{k-j} j^n / (j! (k-j)!)`, one convolution. NTT prime `M`.
pub fn stirling2_row<const M: u32>(n: usize) -> Vec<Z<M>> {
    let (_, inv_fact) = special::factorials::<M>(n + 1);
    let pw = power_table::<M>(n, n as u64);
    let a: Vec<Z<M>> = (0..=n).map(|j| pw[j] * inv_fact[j]).collect();
    let b: Vec<Z<M>> = (0..=n).map(|i| if i % 2 == 0 { inv_fact[i] } else { -inv_fact[i] }).collect();
    let mut c = fps::mul(&a, &b);
    c.truncate(n + 1);
    c
}

/// `S(m, k)` for `m = k..=n`, from `(e^x - 1)^k / k!`. NTT prime `M`.
pub fn stirling2_fixed_k<const M: u32>(n: usize, k: usize) -> Vec<Z<M>> {
    assert!(k <= n, "k must not exceed n");
    let len = n - k + 1;
    let (fact, inv_fact) = special::factorials::<M>(n + 2);
    // (e^x - 1) / x = sum x^i / (i + 1)!
    let t: Vec<Z<M>> = (0..len).map(|i| inv_fact[i + 1]).collect();
    let tk = fps::pow(&t, k as u64, len);
    (0..len).map(|i| tk[i] * fact[k + i] * inv_fact[k]).collect()
}

/// Binomials modulo the small prime `p` by Lucas' theorem.
struct SmallBinom {
    p: u64,
    fact: Vec<u32>,
    inv_fact: Vec<u32>,
}

impl SmallBinom {
    fn new(p: u64) -> Self {
        let mut fact = vec![1u32; p as usize];
        for i in 1..p as usize {
            fact[i] = (fact[i - 1] as u64 * i as u64 % p) as u32;
        }
        let mut inv_fact = vec![1u32; p as usize];
        let mut inv = |a: u64| -> u64 {
            let (mut r, mut b, mut e) = (1u64, a % p, p - 2);
            while e > 0 {
                if e & 1 == 1 {
                    r = r * b % p;
                }
                b = b * b % p;
                e >>= 1;
            }
            r
        };
        inv_fact[p as usize - 1] = inv(fact[p as usize - 1] as u64) as u32;
        for i in (1..p as usize).rev() {
            inv_fact[i - 1] = (inv_fact[i] as u64 * i as u64 % p) as u32;
        }
        let _ = &mut inv;
        SmallBinom { p, fact, inv_fact }
    }

    /// `C(n, k) mod p`.
    fn lucas(&self, mut n: u128, mut k: u128) -> u64 {
        let p = self.p as u128;
        let mut r = 1u64;
        while n > 0 || k > 0 {
            let (a, b) = ((n % p) as usize, (k % p) as usize);
            if b > a {
                return 0;
            }
            r = r * self.fact[a] as u64 % self.p * self.inv_fact[b] as u64 % self.p * self.inv_fact[a - b] as u64 % self.p;
            n /= p;
            k /= p;
        }
        r
    }
}

/// `S(n, k) mod p` for the prime `p` and arbitrary `k <= n < 2^63`, for every query `(n, k)`.
///
/// With `k = c p + r`: `sum_n S(n, k) x^{n-k} = 1 / prod_{j<=k} (1 - j x)` and modulo `p` the
/// product is `(1 - x^{p-1})^c prod_{j=1}^{r} (1 - j x)`; for `r > 0` this gives
/// `S(n, k) = [x^{(n-k) mod (p-1)}] prod_{j=r+1}^{p-1} (1 - j x) * C(c + I, I)` with
/// `I = (n-k) / (p-1)`, and for `r = 0` a single binomial. `O(p^2 + T log_p n)`.
#[allow(clippy::explicit_counter_loop)]
pub fn stirling2_small_prime(p: u64, queries: &[(u64, u64)]) -> Vec<u64> {
    assert!(p >= 2, "p must be prime");
    let sb = SmallBinom::new(p);
    let mut answers = vec![0u64; queries.len()];
    let mut by_r: Vec<Vec<usize>> = vec![Vec::new(); p as usize];
    for (idx, &(n, k)) in queries.iter().enumerate() {
        assert!(k <= n, "k must not exceed n");
        let (c, r, d) = (k / p, k % p, n - k);
        if k == 0 {
            answers[idx] = (n == 0) as u64;
        } else if r == 0 {
            // 1 / (1 - x^{p-1})^c
            answers[idx] = if d % (p - 1) == 0 { sb.lucas((c - 1 + d / (p - 1)) as u128, (d / (p - 1)) as u128) } else { 0 };
        } else {
            by_r[r as usize].push(idx);
        }
    }
    // Q_r = prod_{j=r+1}^{p-1} (1 - j x), from r = p - 1 (empty product) downwards
    let mut q = vec![0u64; p as usize];
    q[0] = 1;
    let mut deg = 0usize;
    for r in (1..p as usize).rev() {
        // here q = Q_r with degree p - 1 - r
        for &idx in &by_r[r] {
            let (n, k) = queries[idx];
            let (c, d) = (k / p, n - k);
            let t = (d % (p - 1)) as usize;
            let coef = if t <= deg { q[t] } else { 0 };
            let i = d / (p - 1);
            answers[idx] = coef * sb.lucas((c + i) as u128, i as u128) % p;
        }
        // Q_{r-1} = Q_r (1 - r x)
        for t in (1..=deg + 1).rev() {
            q[t] = (q[t] + p - q[t - 1] * r as u64 % p) % p;
        }
        deg += 1;
    }
    answers
}

/// `s(n, k) mod p` (signed first kind) for the prime `p` and `k <= n < 2^63`, for every query.
///
/// With `n = a p + b`: `x (x-1) ... (x-n+1) = (x^p - x)^a x (x-1) ... (x-b+1)` modulo `p`, so
/// `s(n, k) = sum_i C(a, i) (-1)^{a-i} s(b, k - a - (p-1) i)` over the at most two admissible `i`.
/// `O(p^2 + T log_p n)`.
#[allow(clippy::needless_range_loop)]
pub fn stirling1_small_prime(p: u64, queries: &[(u64, u64)]) -> Vec<u64> {
    assert!(p >= 2, "p must be prime");
    let sb = SmallBinom::new(p);
    let mut answers = vec![0u64; queries.len()];
    let mut by_b: Vec<Vec<usize>> = vec![Vec::new(); p as usize];
    for (idx, &(n, k)) in queries.iter().enumerate() {
        assert!(k <= n, "k must not exceed n");
        by_b[(n % p) as usize].push(idx);
    }
    // row[j] = s(b, j) mod p, advanced from b to b + 1
    let mut row = vec![0u64; p as usize + 1];
    row[0] = 1;
    for b in 0..p as usize {
        for &idx in &by_b[b] {
            let (n, k) = queries[idx];
            let a = (n / p) as i128;
            let k = k as i128;
            let pm1 = p as i128 - 1;
            let mut total = 0u64;
            if k >= a {
                // admissible i: 0 <= k - a - (p-1) i <= b and 0 <= i <= a
                let lo = ((k - a - b as i128 + pm1 - 1).div_euclid(pm1)).max(0);
                let hi = ((k - a).div_euclid(pm1)).min(a);
                let mut i = lo;
                while i <= hi {
                    let j = (k - a - pm1 * i) as usize;
                    let mut term = sb.lucas(a as u128, i as u128) * row[j] % p;
                    if (a - i) % 2 == 1 {
                        term = (p - term) % p;
                    }
                    total = (total + term) % p;
                    i += 1;
                }
            }
            answers[idx] = total;
        }
        // s(b + 1, j) = s(b, j - 1) - b s(b, j)
        for j in (0..=b + 1).rev() {
            let up = if j >= 1 { row[j - 1] } else { 0 };
            row[j] = (up + p - row[j] * b as u64 % p) % p;
        }
    }
    answers
}
