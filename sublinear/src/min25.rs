//! Min_25's sieve: `sum_{x <= n} f(x)` for a multiplicative `f` whose values on the primes are
//! known through prefix sums, in about `O(n^{3/4} / log n)`.

use super::quotient::Quotients;

/// `sum_{x=1}^n f(x) mod M`, where `prime_f` holds `sum_{p <= x} f(p) mod M` for all `x` in `Q_n`,
/// `primes` lists (at least) the primes up to `sqrt(n)` and `f_pow(p, e)` is `f(p^e) mod M`.
///
/// # Panics
/// Panics unless `M < 2^31` and `n < 2^50`.
pub fn multiplicative_prefix_sum<const M: u64>(
    n: u64,
    primes: &[u32],
    prime_f: &Quotients<u64>,
    f_pow: impl Fn(u64, u32) -> u64,
) -> u64 {
    assert_eq!(prime_f.n(), n, "table for a different n");
    assert!(M < 1 << 31, "the modulus must be below 2^31");
    assert!(n < 1 << 50, "n must be below 2^50");
    let np = primes.len();
    // at_prime[k] = sum of f over the first k primes
    let at_prime: Vec<u64> = std::iter::once(0).chain(primes.iter().map(|&p| prime_f.get(p as u64))).collect();
    let f1: Vec<u64> = primes.iter().map(|&p| f_pow(p as u64, 1)).collect();
    let f2: Vec<u64> = primes.iter().map(|&p| f_pow(p as u64, 2)).collect();
    // ps1[k] = sum_{t < k} f(p_t) at_prime[t + 1],  ps2[k] = sum_{t < k} f(p_t^2)
    let mut ps1 = vec![0u64; np + 1];
    let mut ps2 = vec![0u64; np + 1];
    for t in 0..np {
        ps1[t + 1] = (ps1[t] + f1[t] * at_prime[t + 1]) % M;
        ps2[t + 1] = (ps2[t] + f2[t]) % M;
    }
    let inv: Vec<f64> = primes.iter().map(|&p| 1.0 / p as f64).collect();
    struct Ctx<'a, const M: u64, F: Fn(u64, u32) -> u64> {
        r: u64,
        primes: &'a [u32],
        prime_f: &'a Quotients<u64>,
        at_prime: Vec<u64>,
        f1: Vec<u64>,
        ps1: Vec<u64>,
        ps2: Vec<u64>,
        inv: Vec<f64>,
        f_pow: F,
    }
    impl<const M: u64, F: Fn(u64, u32) -> u64> Ctx<'_, M, F> {
        /// `sum of f(p)` over the primes `p <= v` where `v = floor(n / m)`.
        #[inline(always)]
        fn prime_sum(&self, v: u64, m: u64) -> u64 {
            if v <= self.r {
                self.prime_f.small[v as usize]
            } else {
                self.prime_f.large[m as usize]
            }
        }

        /// Number of primes `p` (in the list) with `p^e <= v`.
        fn count_pow(&self, v: u64, e: u32) -> usize {
            self.primes.partition_point(|&p| (p as u64).pow(e) <= v)
        }

        /// Sum of `f(x)` over `2 <= x <= v`, `v = floor(n / m)`, all of whose prime factors are
        /// `>= primes[j]`.
        fn rec(&self, v: u64, m: u64, j: usize) -> u64 {
            let last = if j == 0 { 1 } else { self.primes[j - 1] as u64 };
            if last >= v {
                return 0;
            }
            // the primes in (primes[j - 1], v] give the terms f(x) with x prime
            let mut res = self.prime_sum(v, m) + M - self.at_prime[j];
            let kend = self.count_pow(v, 2); // primes with p^2 <= v
            let kcube = self.count_pow(v, 3).max(j).min(kend); // p^3 <= v: recursion needed
            for k in j..kcube {
                let p = self.primes[k] as u64;
                let mut pe = p;
                let mut e = 1u32;
                let mut acc = 0u64;
                let below = self.at_prime[k + 1];
                let next_sq = self.primes.get(k + 1).map_or(u64::MAX, |&q| q as u64 * q as u64);
                let mut v2 = ((v as f64 + 0.5) * self.inv[k]) as u64;
                while pe * p <= v {
                    let m2 = m * pe;
                    let fe = if e == 1 { self.f1[k] } else { (self.f_pow)(p, e) };
                    // the primes in (p, v2], plus the composites if p_{k+1}^2 <= v2
                    let sub = if next_sq <= v2 { self.rec(v2, m2, k + 1) } else { self.prime_sum(v2, m2) + M - below };
                    acc = (acc + fe * (sub % M) % M + (self.f_pow)(p, e + 1)) % M;
                    pe *= p;
                    e += 1;
                    v2 /= p;
                }
                res += acc;
            }
            // p^2 <= v < p^3: only e = 1, whose child holds primes only
            let kl = kcube;
            if kl < kend {
                let r = self.r;
                // k < ka: m * p <= r, the value v / p sits in large[m p]
                let ka = kl + self.primes[kl..kend].partition_point(|&p| m * p as u64 <= r);
                let big = 4 * M * M;
                let mut acc = 0u64;
                for k in kl..ka {
                    let idx = m as usize * self.primes[k] as usize;
                    acc += self.f1[k] * self.prime_f.large[idx];
                    acc = acc.min(acc.wrapping_sub(big));
                }
                for k in ka..kend {
                    let d = ((v as f64 + 0.5) * self.inv[k]) as usize;
                    acc += self.f1[k] * self.prime_f.small[d];
                    acc = acc.min(acc.wrapping_sub(big));
                }
                let sub = (self.ps1[kend] + M - self.ps1[kl]) % M;
                res += acc % M + M - sub + self.ps2[kend] + M - self.ps2[kl];
            }
            res % M
        }
    }
    let ctx = Ctx::<M, _> { r: prime_f.root(), primes, prime_f, at_prime, f1, ps1, ps2, inv, f_pow };
    (ctx.rec(n, 1, 0) + 1) % M
}
