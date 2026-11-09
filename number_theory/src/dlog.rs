//! Discrete logarithms: baby-step giant-step with an open-addressing table.

use super::{factorize, gcd, mod_inverse, NONE};

/// A map from nonzero `u32` keys to `u32` values (open addressing, power-of-two size); key and
/// value share a cache line.
struct Table {
    shift: u32,
    slots: Vec<[u32; 2]>,
}

const ABSENT: u32 = u32::MAX;

impl Table {
    fn with_capacity(n: usize) -> Self {
        let slots = (2 * n).next_power_of_two().max(4);
        Table { shift: 64 - slots.trailing_zeros(), slots: vec![[0, 0]; slots] }
    }

    #[inline(always)]
    fn slot(&self, key: u32) -> usize {
        ((key as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> self.shift) as usize
    }

    /// Hints the CPU to load the home slot of `key`.
    #[inline(always)]
    fn prefetch(&self, key: u32) {
        #[cfg(target_arch = "x86_64")]
        // SAFETY: prefetching has no architectural effect, and the address is in bounds.
        unsafe {
            use core::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
            _mm_prefetch::<_MM_HINT_T0>(self.slots.as_ptr().add(self.slot(key)) as *const i8);
        }
        #[cfg(not(target_arch = "x86_64"))]
        let _ = key;
    }

    /// Sets `key -> val`, overwriting.
    #[inline(always)]
    fn insert(&mut self, key: u32, val: u32) {
        let mask = self.slots.len() - 1;
        let mut i = self.slot(key);
        while self.slots[i][0] != 0 && self.slots[i][0] != key {
            i = (i + 1) & mask;
        }
        self.slots[i] = [key, val];
    }

    /// The value of `key`, or [`ABSENT`].
    #[inline(always)]
    fn get(&self, key: u32) -> u32 {
        let mask = self.slots.len() - 1;
        let mut i = self.slot(key);
        loop {
            let [k, v] = self.slots[i];
            if k == key {
                return v;
            }
            if k == 0 {
                return ABSENT;
            }
            i = (i + 1) & mask;
        }
    }
}

/// Multiplication modulo `2 <= m <= 2^32` by Barrett reduction.
#[derive(Clone, Copy)]
struct Barrett {
    m: u64,
    inv: u64, // floor(2^64 / m)
}

impl Barrett {
    fn new(m: u64) -> Self {
        assert!((2..=1 << 32).contains(&m));
        Barrett { m, inv: ((1u128 << 64) / m as u128) as u64 }
    }

    #[inline(always)]
    fn mul(&self, a: u64, b: u64) -> u64 {
        let x = a * b; // a, b < m <= 2^32
        let q = ((x as u128 * self.inv as u128) >> 64) as u64;
        let r = x - q * self.m;
        if r >= self.m {
            r - self.m
        } else {
            r
        }
    }

    fn pow(&self, mut a: u64, mut e: u64) -> u64 {
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

/// The smallest `k >= 0` with `x^k = y (mod m)`, or [`NONE`] (`0^0 = 1`).
///
/// # Panics
/// Panics unless `1 <= m <= 2^32`. Time and memory `O(sqrt(m))`.
pub fn discrete_log(x: u64, y: u64, m: u64) -> u64 {
    assert!((1..=1 << 32).contains(&m), "modulus out of range");
    let (x, mut y, mut m) = (x % m, y % m, m);
    if y == 1 % m {
        return 0;
    }
    // peel the common factors of x and m: x^k = y  <=>  x^(k - t) c = y' (mod m')
    let (mut c, mut t) = (1 % m, 0u64);
    loop {
        let g = gcd(x, m);
        if g == 1 {
            break;
        }
        if y % g != 0 {
            return NONE;
        }
        y /= g;
        m /= g;
        c = c * (x / g % m) % m;
        t += 1;
        if c == y % m {
            return t;
        }
    }
    let y = y % m;
    if m == 1 {
        return NONE; // unreachable: c == y == 0 returns above
    }
    let bar = Barrett::new(m);
    // find the smallest k >= 1 with c x^k = y (mod m), x invertible modulo m
    let s = super::roots::isqrt(m) + 1;
    let mut table = Table::with_capacity(s as usize);
    let mut cur = y;
    for j in 0..s {
        table.insert(cur as u32, j as u32);
        cur = bar.mul(cur, x);
    }
    let xs = bar.pow(x, s);
    let mut cur = c;
    for i in 1..=s + 1 {
        cur = bar.mul(cur, xs);
        let j = table.get(cur as u32);
        if j != ABSENT {
            return t + i * s - j as u64;
        }
    }
    NONE
}

/// Discrete logarithms modulo one fixed prime `p < 2^32` with a given primitive root `g`: a table
/// of `B` baby steps is built once and each query costs `p / B` lookups.
pub struct PrimeLog {
    p: u64,
    bar: Barrett,
    block: u64,
    giant: u64,
    table: Table,
}

impl PrimeLog {
    /// `queries` is the expected number of logarithms asked for (it only tunes the size of the
    /// table, about `sqrt(p * queries)` entries capped at `2^24`).
    pub fn new(p: u64, g: u64, queries: usize) -> Self {
        assert!((2..=1 << 32).contains(&p), "modulus out of range");
        let n = p - 1;
        let block = (((p as f64) * (queries.max(1) as f64)).sqrt() as u64).clamp(1, n).min(1 << 24);
        let bar = Barrett::new(p);
        let mut table = Table::with_capacity(block as usize);
        // insert in groups, prefetching the slots first to overlap the cache misses
        let mut cur = 1 % p;
        let mut j = 0u64;
        while j < block {
            let mut keys = [0u32; 16];
            let cnt = (block - j).min(16) as usize;
            for k in keys.iter_mut().take(cnt) {
                *k = cur as u32;
                table.prefetch(*k);
                cur = bar.mul(cur, g);
            }
            for (t, &k) in keys.iter().take(cnt).enumerate() {
                table.insert(k, (j + t as u64) as u32);
            }
            j += cnt as u64;
        }
        // cur = g^block
        let giant = mod_inverse(cur, p);
        PrimeLog { p, bar, block, giant, table }
    }

    /// The `b` in `0..p-1` with `g^b = a`, for `1 <= a < p`.
    pub fn log(&self, a: u64) -> u64 {
        assert!(a >= 1 && a < self.p, "a out of range");
        let mut cur = a;
        for i in 0..=(self.p - 1) / self.block {
            let j = self.table.get(cur as u32);
            if j != ABSENT {
                return i * self.block + j as u64;
            }
            cur = self.bar.mul(cur, self.giant);
        }
        panic!("g is not a primitive root of {}", self.p)
    }

    /// The logarithms of all the `queries` (each in `1..p`). Faster than calling [`PrimeLog::log`]
    /// one by one: a window of lookups is kept in flight so the cache misses overlap.
    pub fn log_many(&self, queries: &[u64]) -> Vec<u64> {
        const WINDOW: usize = 32;
        let mut ans = vec![0u64; queries.len()];
        let max_steps = (self.p - 1) / self.block;
        // per slot: query index, current element, giant steps taken; idx == usize::MAX when idle
        let mut idx = [usize::MAX; WINDOW];
        let mut cur = [0u64; WINDOW];
        let mut step = [0u64; WINDOW];
        let mut next = 0usize;
        let mut busy = 0usize;
        for s in 0..WINDOW {
            if next < queries.len() {
                assert!(queries[next] >= 1 && queries[next] < self.p, "a out of range");
                (idx[s], cur[s], step[s]) = (next, queries[next], 0);
                self.table.prefetch(cur[s] as u32);
                next += 1;
                busy += 1;
            }
        }
        while busy > 0 {
            for s in 0..WINDOW {
                if idx[s] == usize::MAX {
                    continue;
                }
                let j = self.table.get(cur[s] as u32);
                if j != ABSENT {
                    ans[idx[s]] = step[s] * self.block + j as u64;
                    if next < queries.len() {
                        assert!(queries[next] >= 1 && queries[next] < self.p, "a out of range");
                        (idx[s], cur[s], step[s]) = (next, queries[next], 0);
                        self.table.prefetch(cur[s] as u32);
                        next += 1;
                    } else {
                        idx[s] = usize::MAX;
                        busy -= 1;
                    }
                } else {
                    assert!(step[s] < max_steps, "g is not a primitive root of {}", self.p);
                    cur[s] = self.bar.mul(cur[s], self.giant);
                    step[s] += 1;
                    self.table.prefetch(cur[s] as u32);
                }
            }
        }
        ans
    }
}

/// The `b` in `0..p-1` with `g^b = a (mod p)` for a primitive root `g` of the prime `p < 2^32`
/// and `a != 0` (Pohlig-Hellman with baby-step giant-step on every prime power of `p - 1`).
pub fn discrete_log_prime(a: u64, g: u64, p: u64) -> u64 {
    assert!((2..=1 << 32).contains(&p), "modulus out of range");
    assert!(!a.is_multiple_of(p), "a must be nonzero");
    discrete_log_smooth(a, g, p, &factorize(p - 1))
}

/// The `b` in `0..order` with `g^b = a (mod p)`, where `g` has order exactly `order = prod q^e`
/// (`factors`, which may be a sub-product of the factorisation of `p - 1`) and `a` is a power of
/// `g`. Costs `O(sum e sqrt(q))` group operations.
pub(crate) fn discrete_log_smooth(a: u64, g: u64, p: u64, factors: &[(u64, u32)]) -> u64 {
    let order: u64 = factors.iter().map(|&(q, e)| q.pow(e)).product();
    let bar = Barrett::new(p);
    let (mut res, mut modulus) = (0u64, 1u64);
    for &(q, e) in factors {
        let qe = q.pow(e);
        let gamma = bar.pow(g, order / qe);
        let h = bar.pow(a % p, order / qe);
        let gq = bar.pow(gamma, qe / q); // of order q
        let s = super::roots::isqrt(q) + 1;
        let mut table = Table::with_capacity(s as usize);
        let mut cur = 1u64;
        for j in 0..s {
            table.insert(cur as u32, j as u32);
            cur = bar.mul(cur, gq);
        }
        let gqs_inv = mod_inverse(cur, p); // g_q^(-s)
        let gamma_inv = mod_inverse(gamma, p);
        // x = sum d_k q^k
        let (mut x, mut qk) = (0u64, 1u64);
        for _ in 0..e {
            let v = bar.mul(h, bar.pow(gamma_inv, x));
            let mut w = bar.pow(v, qe / (qk * q)); // = g_q^(d_k)
            let mut d = NONE;
            for i in 0..=s {
                let j = table.get(w as u32);
                if j != ABSENT {
                    d = i * s + j as u64;
                    break;
                }
                w = bar.mul(w, gqs_inv);
            }
            assert!(d != NONE, "g does not have the stated order modulo {p}");
            x += d * qk;
            qk *= q;
        }
        // combine res (mod modulus) with x (mod qe)
        let diff = (x + qe - res % qe) % qe;
        let t = diff * mod_inverse(modulus % qe, qe) % qe;
        res += modulus * t;
        modulus *= qe;
    }
    res
}
