//! Precomputed modular binomial coefficients: `nCr mod p`.
//!
//! [`ModBinomial::new`] builds factorials and inverse factorials up to `max_n`
//! in `O(max_n)` time using Fermat's little theorem for the modular inverse
//! (`modulo` must be prime). Each subsequent [`ModBinomial::ncr`] call is `O(1)`.

/// Binomial-coefficient table modulo a prime `modulo`.
pub struct ModBinomial {
    fact: Vec<i64>,
    inv_fact: Vec<i64>,
    modulo: i64,
}

impl ModBinomial {
    /// Build the table for all `n` in `0..=max_n`. `modulo` must be prime
    /// (the common choice is `1_000_000_007` or `998_244_353`).
    pub fn new(max_n: usize, modulo: i64) -> Self {
        let mut fact = vec![0; max_n + 1];
        let mut inv_fact = vec![0; max_n + 1];

        fact[0] = 1;
        inv_fact[0] = 1;

        for i in 1..=max_n {
            fact[i] = (fact[i - 1] * i as i64) % modulo;
        }

        inv_fact[max_n] = Self::power(fact[max_n], modulo - 2, modulo);
        for i in (1..=max_n).rev() {
            inv_fact[i - 1] = (inv_fact[i] * i as i64) % modulo;
        }

        Self {
            fact,
            inv_fact,
            modulo,
        }
    }

    fn power(mut base: i64, mut exp: i64, modulo: i64) -> i64 {
        let mut res = 1;
        base %= modulo;
        while exp > 0 {
            if exp % 2 == 1 {
                res = (res * base) % modulo;
            }
            base = (base * base) % modulo;
            exp /= 2;
        }
        res
    }

    /// Compute `nCr mod modulo`. Returns `0` when `r > n`.
    pub fn ncr(&self, n: usize, r: usize) -> i64 {
        if r > n {
            return 0;
        }
        let denom = (self.inv_fact[r] * self.inv_fact[n - r]) % self.modulo;
        (self.fact[n] * denom) % self.modulo
    }
}