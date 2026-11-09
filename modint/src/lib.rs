//! `StaticModInt<M>`: an integer modulo the compile-time constant `M` (a prime for division to
//! make sense), stored as a `u32` in `0..M`. All the arithmetic operators are implemented, so
//! number-theoretic code reads like ordinary arithmetic.
//!
//! * `inv` / `/` assert that the divisor is invertible (nonzero modulo a prime); there is no
//!   `Option` in the API.
//! * [`ModInt998244353`] and [`ModInt1000000007`] are the usual aliases.
//!
//! # Example
//!
//! ```
//! use modint::ModInt998244353 as Mint;
//!
//! let a = Mint::new(5);
//! let b = Mint::new(3);
//! assert_eq!((a / b * b).val(), 5);
//! assert_eq!((a - b - b).val(), 998_244_352); // -1
//! assert_eq!(Mint::new(2).pow(10).val(), 1024);
//! assert_eq!(Mint::from(-1i64).val(), 998_244_352);
//! ```

use std::fmt;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Integer modulo `M`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
#[repr(transparent)]
pub struct StaticModInt<const M: u32>(u32);

/// Arithmetic modulo 998 244 353.
pub type ModInt998244353 = StaticModInt<998_244_353>;
/// Arithmetic modulo 1 000 000 007.
pub type ModInt1000000007 = StaticModInt<1_000_000_007>;

impl<const M: u32> StaticModInt<M> {
    /// The modulus.
    pub const MODULUS: u32 = M;

    /// `x mod M`.
    pub fn new(x: u64) -> Self {
        Self((x % M as u64) as u32)
    }

    /// Wrap an already reduced value (`x < M`).
    pub fn raw(x: u32) -> Self {
        assert!(x < M, "raw value {x} is not reduced modulo {M}");
        Self(x)
    }

    pub const fn zero() -> Self {
        Self(0)
    }

    pub const fn one() -> Self {
        Self(1)
    }

    /// The representative in `0..M`.
    pub fn val(self) -> u32 {
        self.0
    }

    /// `self^e` by square and multiply.
    pub fn pow(self, mut e: u64) -> Self {
        let (mut base, mut acc) = (self, Self(1 % M));
        while e > 0 {
            if e & 1 == 1 {
                acc *= base;
            }
            base *= base;
            e >>= 1;
        }
        acc
    }

    /// Multiplicative inverse (extended Euclid, so any modulus works as long as it is coprime).
    ///
    /// # Panics
    /// Panics if the value is not invertible modulo `M` (zero for a prime modulus).
    pub fn inv(self) -> Self {
        let (mut a, mut b) = (self.0 as i64, M as i64);
        let (mut x, mut y) = (1i64, 0i64);
        while b != 0 {
            let q = a / b;
            a -= q * b;
            std::mem::swap(&mut a, &mut b);
            x -= q * y;
            std::mem::swap(&mut x, &mut y);
        }
        assert!(a == 1, "{} is not invertible modulo {M}", self.0);
        Self(x.rem_euclid(M as i64) as u32)
    }
}

impl<const M: u32> fmt::Display for StaticModInt<M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<const M: u32> fmt::Debug for StaticModInt<M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

macro_rules! from_unsigned {
    ($($t:ty),*) => {$(
        impl<const M: u32> From<$t> for StaticModInt<M> {
            fn from(x: $t) -> Self { Self((x as u64 % M as u64) as u32) }
        }
    )*};
}
macro_rules! from_signed {
    ($($t:ty),*) => {$(
        impl<const M: u32> From<$t> for StaticModInt<M> {
            fn from(x: $t) -> Self { Self((x as i64).rem_euclid(M as i64) as u32) }
        }
    )*};
}
from_unsigned!(u8, u16, u32, u64, usize);
from_signed!(i8, i16, i32, i64, isize);

impl<const M: u32> Add for StaticModInt<M> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        let s = self.0 + o.0; // < 2M may overflow u32 only for M > 2^31; use u64 then
        Self(if M > (1 << 31) { ((self.0 as u64 + o.0 as u64) % M as u64) as u32 } else if s >= M { s - M } else { s })
    }
}
impl<const M: u32> Sub for StaticModInt<M> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self(if self.0 >= o.0 { self.0 - o.0 } else { (self.0 as u64 + M as u64 - o.0 as u64) as u32 })
    }
}
impl<const M: u32> Mul for StaticModInt<M> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self((self.0 as u64 * o.0 as u64 % M as u64) as u32)
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl<const M: u32> Div for StaticModInt<M> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        self * o.inv()
    }
}
impl<const M: u32> Neg for StaticModInt<M> {
    type Output = Self;
    fn neg(self) -> Self {
        Self(if self.0 == 0 { 0 } else { M - self.0 })
    }
}
macro_rules! assign_ops {
    ($($tr:ident $f:ident $op:tt),*) => {$(
        impl<const M: u32> $tr for StaticModInt<M> {
            fn $f(&mut self, o: Self) { *self = *self $op o; }
        }
    )*};
}
assign_ops!(AddAssign add_assign +, SubAssign sub_assign -, MulAssign mul_assign *, DivAssign div_assign /);

impl<const M: u32> Sum for StaticModInt<M> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self(0), |a, b| a + b)
    }
}
impl<const M: u32> Product for StaticModInt<M> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self(1 % M), |a, b| a * b)
    }
}
