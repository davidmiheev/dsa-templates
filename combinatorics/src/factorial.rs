//! Factorials of huge arguments modulo an NTT prime.

use modint::StaticModInt;
use polynomial::special::shift_sampling_points;

type Z<const M: u32> = StaticModInt<M>;

/// `n! mod M` for every `n < M`.
///
/// Preprocessing computes `g(t) = prod_{i=1}^{B} (B t + i) = ((t+1)B)! / (tB)!` for
/// `t = 0..M/B` (the polynomial `g` of degree `B` has its first `B + 1` values found by doubling
/// the degree with the shift of sample points, the rest by chunked shifts), hence all `(jB)!`;
/// a query multiplies at most `B / 2` numbers.
pub struct BlockFactorials<const M: u32> {
    block: usize,
    /// `table[j] = (j B)!` for `j B < M`
    table: Vec<Z<M>>,
}

impl<const M: u32> BlockFactorials<M> {
    /// `block` is the block length `B`, a power of two (about `sqrt(M) / 32` is a good value).
    pub fn new(block: usize) -> Self {
        assert!(block.is_power_of_two(), "the block length must be a power of two");
        assert!((block as u64) < M as u64 / 2, "block too long for the modulus");
        let b = block;
        let v = Z::<M>::new(b as u64);
        let v_inv = v.inv();
        // g_d(x) = prod_{i=1}^{d} (v x + i) at x = 0..=d, doubled with g_{2d}(x) = g_d(x) g_d(x + d / v)
        let mut vals = vec![Z::<M>::one(), v + Z::one()];
        let mut d = 1usize;
        while d < b {
            let first = shift_sampling_points(&vals, Z::zero(), 2 * d + 1);
            let second = shift_sampling_points(&vals, Z::<M>::new(d as u64) * v_inv, 2 * d + 1);
            vals = first.iter().zip(&second).map(|(&x, &y)| x * y).collect();
            d *= 2;
        }
        // g_B(t) for t = 0..jmax in chunks of B + 1 points
        let jmax = (M as usize - 1) / b; // (jmax B)! is the last entry
        let mut table = Vec::with_capacity(jmax + 1);
        table.push(Z::<M>::one());
        let mut t = 0usize;
        while t < jmax {
            let chunk = shift_sampling_points(&vals, Z::<M>::new(t as u64), b + 1);
            for &g in chunk.iter().take(jmax - t) {
                let last = *table.last().expect("table is non-empty");
                table.push(last * g);
            }
            t += b + 1;
        }
        BlockFactorials { block: b, table }
    }

    /// `n!` modulo `M`.
    ///
    /// # Panics
    /// Panics if `n >= M`.
    pub fn fact(&self, n: u64) -> Z<M> {
        assert!(n < M as u64, "n must be below the modulus");
        let b = self.block as u64;
        let j = (n / b) as usize;
        let lo = j as u64 * b;
        if n - lo <= b / 2 || j + 1 >= self.table.len() {
            let mut r = self.table[j];
            for i in lo + 1..=n {
                r *= Z::new(i);
            }
            r
        } else {
            // go down from the next multiple of B
            let hi = (j as u64 + 1) * b;
            let mut prod = Z::<M>::one();
            for i in n + 1..=hi {
                prod *= Z::new(i);
            }
            self.table[j + 1] / prod
        }
    }
}
