//! Arbitrary-precision integers in a small base `BASE <= 2^20` (limbs are `u32`).
//!
//! * [`BigUint`]: addition, subtraction, multiplication (schoolbook below 32 limbs, otherwise an
//!   exact NTT with two primes, three for huge operands) and division with remainder (Knuth's algorithm D for short operands,
//!   otherwise Newton's reciprocal with Barrett-style quotient estimates, `O(M(n))`).
//! * [`BigInt`]: the signed wrapper.
//! * Conversion from and to digit strings of a radix with `BASE = radix^k`, e.g.
//!   [`Decimal`] (`BASE = 10^6`) and [`Hex`] (`BASE = 16^5`).
//!
//! # Example
//!
//! ```
//! use bigint::Decimal;
//!
//! let a = Decimal::parse(b"-123456789012345678901234567890");
//! let b = Decimal::parse(b"987654321");
//! assert_eq!(a.mul(&b).to_bytes(), b"-121932631124828532112482853211126352690");
//! let (q, r) = Decimal::parse(b"1000000000000000000000").divrem(&Decimal::parse(b"7777777"));
//! assert_eq!((q.to_bytes(), r.to_bytes()), (b"128571441428572".to_vec(), b"5555556".to_vec()));
//! ```

mod nat;

use std::cmp::Ordering;

/// A non-negative integer in base `BASE`.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct BigUint<const BASE: u32> {
    limbs: Vec<u32>,
}

fn radix_width<const BASE: u32>(radix: u32) -> usize {
    let (mut k, mut p) = (0usize, 1u64);
    while p < BASE as u64 {
        p *= radix as u64;
        k += 1;
    }
    assert!(p == BASE as u64 && k > 0, "BASE must be a power of the radix");
    k
}

fn digit_value(c: u8) -> u32 {
    match c {
        b'0'..=b'9' => (c - b'0') as u32,
        b'a'..=b'f' => (c - b'a') as u32 + 10,
        b'A'..=b'F' => (c - b'A') as u32 + 10,
        _ => panic!("invalid digit {:?}", c as char),
    }
}

impl<const BASE: u32> BigUint<BASE> {
    pub fn zero() -> Self {
        BigUint { limbs: Vec::new() }
    }

    pub fn from_u64(mut x: u64) -> Self {
        let mut limbs = Vec::new();
        while x > 0 {
            limbs.push((x % BASE as u64) as u32);
            x /= BASE as u64;
        }
        BigUint { limbs }
    }

    pub fn is_zero(&self) -> bool {
        self.limbs.is_empty()
    }

    /// Parses digits of the given radix (`BASE` must be a power of it), most significant first.
    pub fn from_digits(s: &[u8], radix: u32) -> Self {
        let k = radix_width::<BASE>(radix);
        let mut limbs = Vec::with_capacity(s.len() / k + 1);
        let mut end = s.len();
        while end > 0 {
            let start = end.saturating_sub(k);
            let v = s[start..end].iter().fold(0u32, |acc, &c| acc * radix + digit_value(c));
            limbs.push(v);
            end = start;
        }
        nat::trim(&mut limbs);
        BigUint { limbs }
    }

    /// The digits in the given radix (lowercase letters), most significant first; `"0"` for zero.
    pub fn to_digits(&self, radix: u32) -> Vec<u8> {
        let k = radix_width::<BASE>(radix);
        if self.is_zero() {
            return vec![b'0'];
        }
        let mut out = Vec::with_capacity(self.limbs.len() * k);
        let top = self.limbs.len() - 1;
        for (idx, &limb) in self.limbs.iter().enumerate().rev() {
            let mut buf = [b'0'; 32];
            let mut v = limb;
            for slot in buf[..k].iter_mut().rev() {
                *slot = char::from_digit(v % radix, radix).expect("digit") as u8;
                v /= radix;
            }
            let digits = &buf[..k];
            if idx == top {
                let skip = digits.iter().take_while(|&&c| c == b'0').count();
                out.extend_from_slice(&digits[skip..]);
            } else {
                out.extend_from_slice(digits);
            }
        }
        out
    }

    pub fn add(&self, other: &Self) -> Self {
        BigUint { limbs: nat::add::<BASE>(&self.limbs, &other.limbs) }
    }

    /// `self - other`, which must be non-negative.
    pub fn sub(&self, other: &Self) -> Self {
        assert!(self.cmp(other) != Ordering::Less, "subtraction would be negative");
        BigUint { limbs: nat::sub::<BASE>(&self.limbs, &other.limbs) }
    }

    pub fn mul(&self, other: &Self) -> Self {
        BigUint { limbs: nat::mul::<BASE>(&self.limbs, &other.limbs) }
    }

    /// `(floor(self / other), self mod other)`; `other` must be positive.
    pub fn divrem(&self, other: &Self) -> (Self, Self) {
        let (q, r) = nat::divrem::<BASE>(&self.limbs, &other.limbs);
        (BigUint { limbs: q }, BigUint { limbs: r })
    }
}

impl<const BASE: u32> Ord for BigUint<BASE> {
    fn cmp(&self, other: &Self) -> Ordering {
        nat::cmp(&self.limbs, &other.limbs)
    }
}

impl<const BASE: u32> PartialOrd for BigUint<BASE> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A signed integer.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct BigInt<const BASE: u32> {
    negative: bool,
    magnitude: BigUint<BASE>,
}

impl<const BASE: u32> BigInt<BASE> {
    fn new(negative: bool, magnitude: BigUint<BASE>) -> Self {
        BigInt { negative: negative && !magnitude.is_zero(), magnitude }
    }

    /// Parses an optional `-` followed by digits of the radix.
    pub fn from_digits(s: &[u8], radix: u32) -> Self {
        match s.first() {
            Some(b'-') => Self::new(true, BigUint::from_digits(&s[1..], radix)),
            _ => Self::new(false, BigUint::from_digits(s, radix)),
        }
    }

    pub fn to_digits(&self, radix: u32) -> Vec<u8> {
        let mut out = Vec::new();
        if self.negative {
            out.push(b'-');
        }
        out.extend(self.magnitude.to_digits(radix));
        out
    }

    pub fn add(&self, other: &Self) -> Self {
        if self.negative == other.negative {
            return Self::new(self.negative, self.magnitude.add(&other.magnitude));
        }
        match self.magnitude.cmp(&other.magnitude) {
            Ordering::Less => Self::new(other.negative, other.magnitude.sub(&self.magnitude)),
            _ => Self::new(self.negative, self.magnitude.sub(&other.magnitude)),
        }
    }

    pub fn mul(&self, other: &Self) -> Self {
        Self::new(self.negative != other.negative, self.magnitude.mul(&other.magnitude))
    }

    /// `(floor(self / other), self mod other)` for non-negative `self` and positive `other`.
    pub fn divrem(&self, other: &Self) -> (Self, Self) {
        assert!(!self.negative && !other.negative, "division needs non-negative operands");
        let (q, r) = self.magnitude.divrem(&other.magnitude);
        (Self::new(false, q), Self::new(false, r))
    }
}

/// Decimal numbers (`BASE = 10^6`).
pub type Decimal = BigInt<1_000_000>;
/// Hexadecimal numbers (`BASE = 16^5`).
pub type Hex = BigInt<1_048_576>;

impl Decimal {
    pub fn parse(s: &[u8]) -> Self {
        Self::from_digits(s, 10)
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.to_digits(10)
    }
}

impl Hex {
    pub fn parse(s: &[u8]) -> Self {
        Self::from_digits(s, 16)
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.to_digits(16)
    }
}
