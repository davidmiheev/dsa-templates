//! Multiplication modulo a runtime `2 <= m < 2^31`.

#[derive(Clone, Copy, Debug)]
pub(crate) struct Barrett {
    pub m: u64,
    inv: u64, // floor(2^64 / m)
}

impl Barrett {
    pub fn new(m: u64) -> Self {
        assert!((2..1 << 31).contains(&m), "modulus out of range");
        Barrett { m, inv: ((1u128 << 64) / m as u128) as u64 }
    }

    /// `x mod m` for any `x < 2^64`.
    #[inline(always)]
    pub fn reduce(&self, x: u64) -> u64 {
        let q = ((x as u128 * self.inv as u128) >> 64) as u64;
        let r = x - q * self.m;
        if r >= self.m {
            r - self.m
        } else {
            r
        }
    }

    #[inline(always)]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        self.reduce(a * b)
    }

    pub fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1 % self.m;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }
}
