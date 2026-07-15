//! Bit manipulation utilities and a tiny fixed-size bitset.
//!
//! The free-standing functions cover the operations that show up over and
//! over in competitive programming:
//!
//! * [`popcount`] — number of set bits
//! * [`trailing_zeros`] / [`trailing_ones`] — index of the lowest set / unset bit
//! * [`highest_bit`] / [`lowest_bit`] — read the most / least significant set bit
//! * [`is_power_of_two`]
//! * [`next_power_of_two`] / [`log2_floor`]
//! * [`subsets`] — iterate over all `2^k` non-empty subsets of a `k`-bit mask
//!
//! [`Bitset`] wraps `Vec<u64>` and exposes the operations you'd expect:
//! `get`, `set`, `reset`, `flip`, plus `count_ones` over the whole set.
//!
//! All inputs are `u64` unless documented otherwise.

/// Number of set bits in `x`.
pub fn popcount(x: u64) -> u32 {
    x.count_ones()
}

/// Index of the lowest set bit. Returns 64 if `x == 0`.
pub fn trailing_zeros(x: u64) -> u32 {
    x.trailing_zeros()
}

/// Index of the lowest unset bit of `!x`. Returns 64 if `x` is all ones.
pub fn trailing_ones(x: u64) -> u32 {
    (!x).trailing_zeros()
}

/// Mask with the highest set bit. Returns 0 if `x == 0`.
pub fn highest_bit(x: u64) -> u64 {
    if x == 0 {
        0
    } else {
        1u64 << (63 - x.leading_zeros())
    }
}

/// Mask with the lowest set bit. Returns 0 if `x == 0`.
pub fn lowest_bit(x: u64) -> u64 {
    x & x.wrapping_neg()
}

/// True iff `x` is a (positive) power of two.
pub fn is_power_of_two(x: u64) -> bool {
    x != 0 && (x & (x - 1)) == 0
}

/// Smallest power of two ≥ `x`. Returns 1 for `x == 0`.
pub fn next_power_of_two(x: u64) -> u64 {
    if x <= 1 {
        return 1;
    }
    1u64 << (64 - (x - 1).leading_zeros())
}

/// Floor of the base-2 logarithm of `x`. Panics if `x == 0`.
pub fn log2_floor(x: u64) -> u32 {
    assert!(x > 0, "log2_floor requires x > 0");
    63 - x.leading_zeros()
}

/// Iterate over all `2^k` non-empty subsets of a `k`-bit mask.
///
/// `k` is the bit-width of the universe (e.g. 4 for a 4-bit mask).
pub fn subsets(k: u32) -> impl Iterator<Item = u64> {
    (1u64..(1u64 << k)).map(move |m| m)
}

/// Fix-sized bitset backed by `Vec<u64>`.
#[derive(Debug, Clone)]
pub struct Bitset {
    /// Number of bits (must be ≤ `64 * words.len()`).
    n: usize,
    /// Storage, little-endian: bit `i` is in `words[i / 64]` at bit `i % 64`.
    words: Vec<u64>,
}

impl Bitset {
    /// Construct a new bitset of `n` bits, all zero.
    pub fn new(n: usize) -> Self {
        let nwords = n.div_ceil(64);
        Self { n, words: vec![0; nwords] }
    }

    /// Get bit at index `i`.
    pub fn get(&self, i: usize) -> bool {
        assert!(i < self.n, "index out of range");
        (self.words[i / 64] >> (i % 64)) & 1 != 0
    }

    /// Set bit at index `i` to `1`.
    pub fn set(&mut self, i: usize) {
        assert!(i < self.n, "index out of range");
        self.words[i / 64] |= 1u64 << (i % 64);
    }

    /// Reset bit at index `i` to `0`.
    pub fn reset(&mut self, i: usize) {
        assert!(i < self.n, "index out of range");
        self.words[i / 64] &= !(1u64 << (i % 64));
    }

    /// Flip bit at index `i`.
    pub fn flip(&mut self, i: usize) {
        assert!(i < self.n, "index out of range");
        self.words[i / 64] ^= 1u64 << (i % 64);
    }

    /// Total number of set bits.
    pub fn count_ones(&self) -> u32 {
        self.words.iter().map(|w| w.count_ones()).sum()
    }
}
