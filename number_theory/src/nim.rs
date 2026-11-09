//! Nim multiplication of 64-bit numbers.

/// Nim products of `u64`s. Built once (a 64 KiB table of the 8-bit products), then each product
/// takes a few dozen table lookups.
pub struct NimMultiplier {
    table: Vec<u8>,
}

impl Default for NimMultiplier {
    fn default() -> Self {
        Self::new()
    }
}

/// The nim product of `a`, `b < 2^bits` (`bits` a power of two), straight from the recursion.
fn slow(a: u64, b: u64, bits: u32) -> u64 {
    if bits == 1 {
        return a & b;
    }
    let h = bits / 2;
    let (mask, x) = ((1u64 << h) - 1, 1u64 << (h - 1));
    let (a0, a1, b0, b1) = (a & mask, a >> h, b & mask, b >> h);
    let c = slow(a0, b0, h);
    let d = slow(a1, b1, h);
    let e = slow(a0 ^ a1, b0 ^ b1, h);
    ((e ^ c) << h) | (c ^ slow(d, x, h))
}

impl NimMultiplier {
    pub fn new() -> Self {
        let mut table = vec![0u8; 1 << 16];
        for a in 0..256u64 {
            for b in 0..256u64 {
                table[(a << 8 | b) as usize] = slow(a, b, 8) as u8;
            }
        }
        NimMultiplier { table }
    }

    /// The nim product `a ⊗ b`.
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        self.rec(a, b, 64)
    }

    fn rec(&self, a: u64, b: u64, bits: u32) -> u64 {
        if a <= 1 || b <= 1 {
            return a * b;
        }
        if bits == 8 {
            return self.table[(a << 8 | b) as usize] as u64;
        }
        let h = bits / 2;
        let (mask, x) = (u64::MAX >> (64 - h), 1u64 << (h - 1));
        let (a0, a1, b0, b1) = (a & mask, a >> h, b & mask, b >> h);
        let c = self.rec(a0, b0, h);
        let d = self.rec(a1, b1, h);
        let e = self.rec(a0 ^ a1, b0 ^ b1, h);
        ((e ^ c) << h) | (c ^ self.rec(d, x, h))
    }
}
