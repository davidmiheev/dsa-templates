//! Binomial coefficients modulo composite moduli, primes and q-binomials.

use crate::barrett::Barrett;
use number_theory::{factorize, mod_inverse};

struct PrimePower {
    p: u64,
    e: u32,
    pe: u64,
    /// `prod[i]` = product of the integers in `1..=i` coprime to `p`, modulo `pe`
    prod: Vec<u32>,
}

impl PrimePower {
    /// `n!` with every factor `p` removed, modulo `p^e`.
    fn fact_without_p(&self, mut n: u64) -> u64 {
        let mut r = 1u64;
        let full = self.prod[self.pe as usize - 1] as u64; // +-1 modulo pe
        while n > 0 {
            if full != 1 && (n / self.pe) % 2 == 1 {
                r = r * full % self.pe;
            }
            r = r * self.prod[(n % self.pe) as usize] as u64 % self.pe;
            n /= self.p;
        }
        r
    }

    /// The exponent of `p` in `n!`.
    fn valuation(&self, mut n: u64) -> u64 {
        let mut v = 0;
        while n > 0 {
            n /= self.p;
            v += n;
        }
        v
    }

    fn binom(&self, n: u64, k: u64) -> u64 {
        let v = self.valuation(n) - self.valuation(k) - self.valuation(n - k);
        if v >= self.e as u64 {
            return 0;
        }
        let num = self.fact_without_p(n);
        let den = self.fact_without_p(k) * self.fact_without_p(n - k) % self.pe;
        let mut r = num * mod_inverse(den, self.pe) % self.pe;
        for _ in 0..v {
            r = r * self.p % self.pe;
        }
        r
    }
}

/// `C(n, k) mod m` for any `m <= 10^7` and `n` up to `2^63`: the prime power factors of `m` are
/// handled with the factorial-without-`p` tables and combined by the Chinese remainder theorem.
/// Preprocessing `O(m)`, a query `O(sum_p log_p n)`.
pub struct BinomialMod {
    m: u64,
    parts: Vec<PrimePower>,
}

impl BinomialMod {
    pub fn new(m: u64) -> Self {
        assert!(m >= 1, "the modulus must be positive");
        assert!(m <= 10_000_000, "the modulus is too large for the tables");
        let parts = factorize(m)
            .into_iter()
            .map(|(p, e)| {
                let pe = p.pow(e);
                let mut prod = vec![1u32; pe as usize];
                for i in 1..pe as usize {
                    prod[i] = if (i as u64).is_multiple_of(p) { prod[i - 1] } else { (prod[i - 1] as u64 * i as u64 % pe) as u32 };
                }
                PrimePower { p, e, pe, prod }
            })
            .collect();
        BinomialMod { m, parts }
    }

    /// `C(n, k) mod m` (0 when `k > n`).
    pub fn binom(&self, n: u64, k: u64) -> u64 {
        if k > n || self.m == 1 {
            return 0;
        }
        let (mut x, mut modulus) = (0u64, 1u64);
        for part in &self.parts {
            let r = part.binom(n, k);
            // x = x (mod modulus), x = r (mod pe)
            let t = (r + part.pe - x % part.pe) % part.pe * mod_inverse(modulus % part.pe, part.pe) % part.pe;
            x += modulus * t;
            modulus *= part.pe;
        }
        x
    }
}

/// `C(n, k) mod m` for a prime `m` and `n < m`, from factorial tables up to `limit`.
pub struct BinomialPrime {
    bar: Barrett,
    fact: Vec<u32>,
    inv_fact: Vec<u32>,
}

impl BinomialPrime {
    /// Tables for `n <= limit` (`limit < m`).
    pub fn new(m: u64, limit: usize) -> Self {
        assert!((limit as u64) < m, "the table must stay below the modulus");
        let bar = Barrett::new(m);
        let mut fact = vec![1u32; limit + 1];
        for i in 1..=limit {
            fact[i] = bar.mul(fact[i - 1] as u64, i as u64) as u32;
        }
        let mut inv_fact = vec![1u32; limit + 1];
        inv_fact[limit] = bar.pow(fact[limit] as u64, m - 2) as u32;
        for i in (1..=limit).rev() {
            inv_fact[i - 1] = bar.mul(inv_fact[i] as u64, i as u64) as u32;
        }
        BinomialPrime { bar, fact, inv_fact }
    }

    /// `C(n, k) mod m` (0 when `k > n`).
    #[inline]
    pub fn binom(&self, n: usize, k: usize) -> u64 {
        if k > n {
            return 0;
        }
        let x = self.bar.mul(self.fact[n] as u64, self.inv_fact[k] as u64);
        self.bar.mul(x, self.inv_fact[n - k] as u64)
    }

    /// `n!` from the table.
    pub fn fact(&self, n: usize) -> u64 {
        self.fact[n] as u64
    }
}

/// `[n choose k]_q mod m` for a prime `m`, `0 <= q < m` and `n < m`, by the q-Lucas theorem: with
/// `d` the multiplicative order of `q`, `[n choose k]_q = C(n / d, k / d) [n % d choose k % d]_q`.
pub struct QBinomialPrime {
    bar: Barrett,
    q_is_zero: bool,
    d: u64,
    /// ordinary binomials for `n / d`
    outer: BinomialPrime,
    q_fact: Vec<u32>,
    inv_q_fact: Vec<u32>,
}

impl QBinomialPrime {
    /// Tables for `n <= limit` (`limit < m`).
    pub fn new(m: u64, q: u64, limit: usize) -> Self {
        assert!(q < m && (limit as u64) < m, "q and the table limit must be below the modulus");
        let bar = Barrett::new(m);
        if q == 0 {
            return QBinomialPrime { bar, q_is_zero: true, d: 1, outer: BinomialPrime::new(m, 0), q_fact: Vec::new(), inv_q_fact: Vec::new() };
        }
        // order of q in the multiplicative group: a divisor of m - 1
        let mut d = m - 1;
        for (p, _) in factorize(m - 1) {
            while d.is_multiple_of(p) && bar.pow(q, d / p) == 1 {
                d /= p;
            }
        }
        let lim_q = (d - 1).min(limit as u64) as usize;
        let mut q_fact = vec![1u32; lim_q + 1];
        let mut bracket = 0u64; // [i]_q
        for i in 1..=lim_q {
            bracket = (bar.mul(bracket, q) + 1) % m;
            q_fact[i] = bar.mul(q_fact[i - 1] as u64, bracket) as u32;
        }
        let mut inv_q_fact = vec![1u32; lim_q + 1];
        inv_q_fact[lim_q] = bar.pow(q_fact[lim_q] as u64, m - 2) as u32;
        let mut bracket_at = vec![0u64; lim_q + 1];
        for i in 1..=lim_q {
            bracket_at[i] = (bar.mul(bracket_at[i - 1], q) + 1) % m;
        }
        for i in (1..=lim_q).rev() {
            inv_q_fact[i - 1] = bar.mul(inv_q_fact[i] as u64, bracket_at[i]) as u32;
        }
        let outer = BinomialPrime::new(m, (limit as u64 / d) as usize);
        QBinomialPrime { bar, q_is_zero: false, d, outer, q_fact, inv_q_fact }
    }

    /// `[n choose k]_q mod m` (0 when `k > n`).
    pub fn binom(&self, n: usize, k: usize) -> u64 {
        if k > n {
            return 0;
        }
        if self.q_is_zero {
            return 1;
        }
        let d = self.d as usize;
        let (r1, r2) = (n % d, k % d);
        if r2 > r1 {
            return 0;
        }
        let low = self.bar.mul(self.q_fact[r1] as u64, self.inv_q_fact[r2] as u64);
        let low = self.bar.mul(low, self.inv_q_fact[r1 - r2] as u64);
        self.bar.mul(self.outer.binom(n / d, k / d), low)
    }
}
