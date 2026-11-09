//! Algebraic traits shared by the data-structure crates.
//!
//! * [`Monoid`]: an associative operation with an identity (sum, min, function composition...).
//! * [`MapMonoid`]: a monoid of *actions* `F` acting on the elements of a monoid `M`, as needed
//!   by lazy propagation (range add, range assign, range affine, ...).
//! * [`Group`]: a monoid with inverses (potentials, prefix differences).
//!
//! The traits are *types with associated functions* (no values), so each concrete monoid is a
//! zero-sized marker type and every call is resolved statically, with no closure plumbing.
//!
//! # Example
//!
//! ```
//! use algebra::Monoid;
//!
//! /// Range minimum over `i64`.
//! enum Min {}
//! impl Monoid for Min {
//!     type S = i64;
//!     fn identity() -> i64 { i64::MAX }
//!     fn op(a: &i64, b: &i64) -> i64 { *a.min(b) }
//! }
//! assert_eq!(Min::op(&3, &Min::identity()), 3);
//! ```

/// An associative binary operation with an identity element.
pub trait Monoid {
    /// The element type.
    type S: Clone;
    /// The identity: `op(identity(), x) == x == op(x, identity())`.
    fn identity() -> Self::S;
    /// The associative operation (not necessarily commutative).
    fn op(a: &Self::S, b: &Self::S) -> Self::S;
}

/// A monoid `F` of maps acting on the elements of the monoid `M`.
///
/// Laws: `apply(f, op(x, y)) == op(apply(f, x), apply(f, y))` (the action is a monoid
/// homomorphism), `apply(id_map(), x) == x`, and
/// `apply(compose(f, g), x) == apply(f, apply(g, x))`.
pub trait MapMonoid {
    /// The monoid being acted on.
    type M: Monoid;
    /// The action type.
    type F: Clone;
    /// The identity action.
    fn id_map() -> Self::F;
    /// Apply `f` to `x`.
    fn apply(f: &Self::F, x: &<Self::M as Monoid>::S) -> <Self::M as Monoid>::S;
    /// The action "do `g` first, then `f`" (that is, `f ∘ g`).
    fn compose(f: &Self::F, g: &Self::F) -> Self::F;
}

/// A monoid in which every element has an inverse.
pub trait Group: Monoid {
    /// `op(x, inv(x)) == identity() == op(inv(x), x)`.
    fn inv(x: &Self::S) -> Self::S;
}
