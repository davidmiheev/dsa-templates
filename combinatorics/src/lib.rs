//! Enumerative combinatorics.
//!
//! * [`BlockFactorials`]: `n! mod M` for any `n < M` (NTT prime `M`) in `O(sqrt(M) log M)`
//!   preprocessing plus `O(B)` per query, by evaluating `prod_{i=1}^B (B x + i)` at consecutive
//!   integers with the sample-point shift.
//! * [`BinomialMod`] (any modulus, huge `n`), [`BinomialPrime`] (factorial tables) and
//!   [`QBinomialPrime`] (q-Lucas).
//! * [`stirling1_row`], [`stirling2_row`], [`stirling1_fixed_k`], [`stirling2_fixed_k`] and the
//!   small-prime / huge-`n` versions [`stirling1_small_prime`], [`stirling2_small_prime`].
//! * [`bell_numbers`], [`partition_numbers`], [`subset_sum_counts`], [`derangements_mod`],
//!   [`distinct_subsequences`].
//! * [`count_increasing_sequences`]: nondecreasing sequences between two bounding sequences.
//!
//! # Example
//!
//! ```
//! use combinatorics::{partition_numbers, BinomialMod};
//!
//! assert_eq!(partition_numbers::<998_244_353>(5).iter().map(|x| x.val()).collect::<Vec<_>>(), vec![1, 1, 2, 3, 5, 7]);
//! assert_eq!(BinomialMod::new(1000).binom(10, 5), 252);
//! ```

mod barrett;
mod binomial;
mod factorial;
mod lattice;
mod sequences;
mod stirling;

pub use binomial::{BinomialMod, BinomialPrime, QBinomialPrime};
pub use factorial::BlockFactorials;
pub use lattice::{count_increasing_sequences, count_with_base};
pub use sequences::{
    bell_numbers, derangements_mod, distinct_subsequences, partition_numbers, subset_sum_counts,
};
pub use stirling::{
    stirling1_fixed_k, stirling1_row, stirling1_small_prime, stirling2_fixed_k, stirling2_row,
    stirling2_small_prime,
};
