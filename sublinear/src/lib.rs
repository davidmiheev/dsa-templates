//! Sums over the floor quotients `Q_n = { floor(n / i) }` in `O(n^{2/3})` – `O(n^{3/4})` time.
//!
//! * [`Quotients`]: a table indexed by the elements of `Q_n`.
//! * [`prime_pi`], [`prime_count_and_sum`]: Lucy_Hedgehog's algorithm.
//! * [`totient_sum`], [`mertens`], [`squarefree_count`].
//! * [`multiplicative_prefix_sum`]: Min_25's sieve.
//! * [`dirichlet_product`], [`dirichlet_inverse`]: prefix sums of Dirichlet products at all
//!   quotients, `O(n^{3/4})`; [`sum_by_division`]: `O(n^{2/3})` prefix sums given a convolution
//!   partner.

mod dirichlet;
mod division;
mod functions;
mod lucy;
mod min25;
mod quotient;

pub use dirichlet::{dirichlet_inverse, dirichlet_product};
pub use division::{sum_by_division, Ring};
pub use functions::{mertens, squarefree_count, totient_sum};
pub use lucy::{prime_count_and_sum, prime_pi};
pub use min25::multiplicative_prefix_sum;
pub use quotient::{isqrt, quotient_values, Quotients};
