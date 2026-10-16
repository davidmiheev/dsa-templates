//! Polynomials and formal power series over a prime field `Z/MZ` with an NTT-friendly modulus
//! (998244353 by default), as plain `Vec<StaticModInt<M>>` slices (index = degree).
//!
//! Series functions take the number of terms `n` wanted and return exactly `n` coefficients;
//! preconditions (`inv` needs a non-zero constant term, `exp` a zero one, ...) are asserted.
//! "No answer" is an empty vector (e.g. [`sqrt`] of a non-residue).
//!
//! * [`sparse`]: `O(N K)` versions of inv / log / exp / pow / sqrt for series with `K` terms.
//! * [`special`]: Taylor shift, shifting sampling points, geometric-sequence evaluation and
//!   interpolation, prefix sums of a polynomial, products of many polynomials.
//! * [`multipoint`]: subproduct tree, multipoint evaluation, interpolation, Newton basis.
//! * [`fps`]: [`fps::mul`], [`fps::inv`], [`fps::log`], [`fps::exp`], [`fps::pow`], [`fps::sqrt`],
//!   derivative / integral, division with remainder.
//!
//! # Example
//!
//! ```
//! use modint::ModInt998244353 as Mint;
//! use polynomial::fps;
//!
//! let f: Vec<Mint> = [1, 1, 0, 0, 0].iter().map(|&x| Mint::new(x)).collect(); // 1 + x
//! // 1 / (1 + x) = 1 - x + x^2 - ...
//! let g = fps::inv(&f, 5);
//! assert_eq!(g, [1, -1, 1, -1, 1].map(Mint::from).to_vec());
//! // exp(log(f)) = f
//! assert_eq!(fps::exp(&fps::log(&f, 5), 5), f);
//! ```

pub mod fps;
pub mod multipoint;
pub mod sparse;
pub mod special;

pub use modint::StaticModInt;
