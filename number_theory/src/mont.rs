//! Montgomery arithmetic modulo an odd `n < 2^64`.

#[derive(Clone, Copy)]
pub(crate) struct Mont {
    pub n: u64,
    inv: u64, // n * inv = 1 mod 2^64
    r2: u64,  // 2^128 mod n
    pub one: u64, // 2^64 mod n
}

impl Mont {
    pub fn new(n: u64) -> Self {
        debug_assert!(n % 2 == 1 && n > 1);
        let mut inv = n; // correct to 3 bits, each step doubles the precision
        for _ in 0..5 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
        }
        let one = ((1u128 << 64) % n as u128) as u64;
        let r2 = ((one as u128 * one as u128) % n as u128) as u64;
        Mont { n, inv, r2, one }
    }

    #[inline(always)]
    fn reduce(&self, t: u128) -> u64 {
        let m = (t as u64).wrapping_mul(self.inv);
        let mn_hi = ((m as u128 * self.n as u128) >> 64) as u64;
        let th = (t >> 64) as u64;
        if th >= mn_hi {
            th - mn_hi
        } else {
            th.wrapping_sub(mn_hi).wrapping_add(self.n)
        }
    }

    #[inline(always)]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        self.reduce(a as u128 * b as u128)
    }

    pub fn to_mont(self, x: u64) -> u64 {
        self.mul(x % self.n, self.r2)
    }

    /// `a^e` for `a` in Montgomery form.
    pub fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = self.one;
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
