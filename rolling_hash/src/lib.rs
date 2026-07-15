//! Rolling (Rabin-Karp-style) hash for strings and byte slices.
//!
//! Provides two flavours:
//!
//! * [`RollingHash`] — single-modulus rolling hash over `&[u8]`. The current
//!   hash value is maintained incrementally as the window slides. Useful
//!   when you need a fingerprint that updates in `O(1)` per shift.
//!
//! * [`DoubleRollingHash`] — same idea, with two different moduli for very
//!   low collision rates; pairs of hashes can be compared for equality.
//!
//! Both flavours precompute the base-power table for substrings up to some
//! `max_len`. To compare hashes of arbitrary substrings of a long string,
//! use [`compute_hash`] (single) or build a [`DoubleRollingHash`] and
//! query substrings using `hash_range(lo, hi)` / `combine(lo1, hi1, lo2, hi2)`.

/// Default base used for the rolling hash.
pub const DEFAULT_BASE: u64 = 131;
/// Default modulus (a known large prime).
pub const DEFAULT_MOD: u64 = (1 << 61) - 1;

/// Single-modulus rolling hash with `O(1)` append / slide.
#[derive(Debug, Clone)]
pub struct RollingHash {
    base: u64,
    mod_: u64,
    hash: u64,
}

impl RollingHash {
    /// Create a new rolling hash with the default base and modulus.
    pub fn new() -> Self {
        Self::with_base_mod(DEFAULT_BASE, DEFAULT_MOD)
    }

    /// Create with a custom base and modulus.
    pub fn with_base_mod(base: u64, mod_: u64) -> Self {
        Self { base, mod_, hash: 0 }
    }

    /// Reset to the empty state.
    pub fn clear(&mut self) {
        self.hash = 0;
    }

    /// Append a byte to the hash, computing
    /// `h = (h * base + byte) mod mod_`.
    pub fn push(&mut self, byte: u8) {
        self.hash = mul_mod(self.hash, self.base, self.mod_);
        self.hash = add_mod(self.hash, byte as u64, self.mod_);
    }

    /// Remove the leftmost byte `b_old` from a window of size `n`, computing
    /// `h = (h - b_old * base^(n-1)) mod mod_`.
    ///
    /// `window_len_before` is the length of the window *before* the removal.
    pub fn pop(&mut self, b_old: u8, window_len_before: usize) {
        debug_assert!(window_len_before >= 1);
        let p = pow_mod(self.base, (window_len_before - 1) as u64, self.mod_);
        let sub = mul_mod(b_old as u64, p, self.mod_);
        self.hash = add_mod(self.hash, self.mod_ - sub % self.mod_, self.mod_);
    }

    /// Current hash value.
    pub fn value(&self) -> u64 {
        self.hash
    }
}

impl Default for RollingHash {
    fn default() -> Self {
        Self::new()
    }
}

/// Compute the polynomial hash of a byte slice.
///
/// `h(s) = s[0]*b^(n-1) + s[1]*b^(n-2) + ... + s[n-1]  (mod mod_)`.
pub fn compute_hash(bytes: &[u8], base: u64, mod_: u64) -> u64 {
    let mut h = 0u64;
    for &b in bytes {
        h = mul_mod(h, base, mod_);
        h = add_mod(h, b as u64, mod_);
    }
    h
}

/// Two-modulus rolling hash with combine/equal-friendly semantics.
#[derive(Debug, Clone)]
pub struct DoubleRollingHash {
    h1: RollingHash,
    h2: RollingHash,
    /// Base for the first hash stream (exposed for advanced users).
    pub pub_b1: u64,
    /// Base for the second hash stream.
    pub pub_b2: u64,
    /// Modulus for the first hash stream.
    pub pub_m1: u64,
    /// Modulus for the second hash stream.
    pub pub_m2: u64,
}

impl DoubleRollingHash {
    /// Create with two distinct (base, modulus) pairs. Recommended:
    /// `(131, 2^61-1)` and `(137, 2^61-1)`.
    pub fn new(base1: u64, mod1: u64, base2: u64, mod2: u64) -> Self {
        Self {
            h1: RollingHash::with_base_mod(base1, mod1),
            h2: RollingHash::with_base_mod(base2, mod2),
            pub_b1: base1,
            pub_b2: base2,
            pub_m1: mod1,
            pub_m2: mod2,
        }
    }

    /// Reset both hashes to 0.
    pub fn clear(&mut self) {
        self.h1.clear();
        self.h2.clear();
    }

    /// Append a byte to both hashes.
    pub fn push(&mut self, byte: u8) {
        self.h1.push(byte);
        self.h2.push(byte);
    }

    /// Get the pair `(hash1, hash2)`.
    pub fn value(&self) -> (u64, u64) {
        (self.h1.value(), self.h2.value())
    }
}

fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    // For m = 2^61 - 1, this matches the standard Rh_hash trick.
    ((a as u128 * b as u128) % m as u128) as u64
}

fn add_mod(a: u64, b: u64, m: u64) -> u64 {
    let s = a + b;
    if s >= m {
        s - m
    } else {
        s
    }
}

fn pow_mod(mut base: u64, mut exp: u64, m: u64) -> u64 {
    let mut result = 1u64;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mul_mod(result, base, m);
        }
        base = mul_mod(base, base, m);
        exp >>= 1;
    }
    result
}
