//! String algorithms.
//!
//! * [`suffix_array`] (SA-IS, `O(n)`), [`suffix_array_bytes`], [`lcp_array`] (Kasai),
//!   [`count_distinct_substrings`], [`longest_common_substring`].
//! * [`z_algorithm`], [`manacher`] (maximal palindromes around all `2n - 1` centres),
//!   [`lyndon_factorization`].
//! * [`AhoCorasick`] (trie with suffix links), [`Eertree`] (palindromic tree) and
//!   [`DequeEertree`] (palindromic tree with push / pop at both ends).
//! * [`runs`]: all maximal repetitions by the Lyndon-array method of Bannai et al.
//! * [`wildcard_match`]: matching with wildcards in both strings by three convolutions.
//!
//! # Example
//!
//! ```
//! use strings::*;
//!
//! assert_eq!(suffix_array_bytes(b"abracadabra"), vec![10, 7, 0, 3, 5, 8, 1, 4, 6, 9, 2]);
//! assert_eq!(z_algorithm(b"aabaa"), vec![5, 1, 0, 2, 1]);
//! assert_eq!(count_distinct_substrings(b"aab"), 5);
//! assert_eq!(lyndon_factorization(b"banana"), vec![0, 1, 3, 5, 6]);
//! ```

mod aho_corasick;
mod deque_eertree;
mod eertree;
mod misc;
mod runs;
mod suffix;

pub use aho_corasick::AhoCorasick;
pub use deque_eertree::DequeEertree;
pub use eertree::Eertree;
pub use misc::{lyndon_factorization, manacher, wildcard_match, z_algorithm};
pub use runs::runs;
pub use suffix::{
    count_distinct_substrings, lcp_array, longest_common_substring, suffix_array, suffix_array_bytes,
};
