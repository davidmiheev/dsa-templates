//! Prefix sums of a multiplicative function by the division trick, in `O(n^{2/3})`.

use modint::StaticModInt;
use std::ops::{Add, Mul, Sub};

/// A commutative ring for [`sum_by_division`].
pub trait Ring: Copy + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> {
    fn zero() -> Self;
    fn from_u64(x: u64) -> Self;
}

impl Ring for i64 {
    fn zero() -> Self {
        0
    }
    fn from_u64(x: u64) -> Self {
        x as i64
    }
}

impl<const M: u32> Ring for StaticModInt<M> {
    fn zero() -> Self {
        StaticModInt::new(0)
    }
    fn from_u64(x: u64) -> Self {
        StaticModInt::new(x)
    }
}

/// `F(n) = sum_{i <= n} f(i)`, given the prefix sums `small(i) = F(i)` for `i <= bound`
/// (`bound` at least about `n^{2/3}`), the prefix sums `g_prefix` of a function `g` with `g(1) = 1`
/// and the prefix sums `h_prefix` of the Dirichlet product `h = f * g`.
///
/// From `H(v) = sum_{d <= v} g(d) F(v / d)`: `F(v) = H(v) - sum_{d >= 2} g(d) F(v / d)` for
/// the `n / bound` large values `v = n / i`, `O(sqrt(v))` each.
pub fn sum_by_division<T: Ring>(
    n: u64,
    bound: u64,
    small: impl Fn(u64) -> T,
    g_prefix: impl Fn(u64) -> T,
    h_prefix: impl Fn(u64) -> T,
) -> T {
    if n <= bound {
        return small(n);
    }
    let top = n / (bound + 1); // v = n / i exceeds `bound` exactly for i <= top
    let mut big = vec![T::zero(); top as usize + 1];
    for i in (1..=top).rev() {
        let v = n / i;
        let mut res = h_prefix(v);
        let mut l = 2u64;
        let mut g_before = g_prefix(1);
        while l <= v {
            let q = v / l;
            let r = v / q;
            let g_now = g_prefix(r);
            let fq = if q <= bound { small(q) } else { big[(i * l) as usize] };
            res = res - (g_now - g_before) * fq;
            g_before = g_now;
            l = r + 1;
        }
        big[i as usize] = res;
    }
    big[1]
}
