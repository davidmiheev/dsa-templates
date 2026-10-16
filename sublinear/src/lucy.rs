//! Lucy_Hedgehog's prime counting: `pi(x)` and `sum of primes <= x` for every `x` in `Q_n`
//! in `O(n^{3/4})`.

use super::quotient::{isqrt, Quotients};

/// `pi(x)` for all `x` in `Q_n`.
#[allow(clippy::needless_range_loop)]
pub fn prime_pi(n: u64) -> Quotients<u64> {
    let r = isqrt(n);
    let ru = r as usize;
    let mut small: Vec<u64> = (0..=ru as u64).map(|v| v.saturating_sub(1)).collect();
    let mut large: Vec<u64> = (0..=ru as u64).map(|i| n.checked_div(i).map_or(0, |q| q - 1)).collect();
    for p in 2..=r {
        let sp = small[(p - 1) as usize];
        if small[p as usize] == sp {
            continue; // p is not prime
        }
        let p2 = p * p;
        let lim = r.min(n / p2) as usize;
        let split = (r / p).min(lim as u64) as usize; // i <= split: i * p <= r
        let pu = p as usize;
        for i in 1..=split {
            large[i] -= large[i * pu] - sp;
        }
        for i in split + 1..=lim {
            let d = i as u64 * p;
            large[i] -= small[(n / d) as usize] - sp;
        }
        if p2 <= r {
            for v in (p2..=r).rev() {
                small[v as usize] -= small[(v / p) as usize] - sp;
            }
        }
    }
    let mut q = Quotients::new(n, 0u64);
    q.small = small;
    q.large = large;
    q
}

/// `(pi(x), sum of primes <= x mod M)` for all `x` in `Q_n`, for a modulus `M < 2^32`.
pub fn prime_count_and_sum<const M: u64>(n: u64) -> (Quotients<u64>, Quotients<u64>) {
    assert!(M < 1 << 32, "the modulus must be below 2^32");
    let r = isqrt(n);
    let ru = r as usize;
    let tri = |v: u64| -> u32 { ((v as u128 * (v as u128 + 1) / 2 - 1) % M as u128) as u32 };
    let mut cs: Vec<u32> = (0..=r).map(|v| v.saturating_sub(1) as u32).collect();
    let mut cl: Vec<u64> = (0..=r).map(|i| n.checked_div(i).map_or(0, |q| q - 1)).collect();
    let mut ss: Vec<u32> = (0..=r).map(|v| if v == 0 { 0 } else { tri(v) }).collect();
    let mut sl: Vec<u32> = (0..=r).map(|i| n.checked_div(i).map_or(0, tri)).collect();
    // n / i as a float, for dividing by p without an integer division
    let vf: Vec<f64> = (0..=r).map(|i| n.checked_div(i).map_or(0.0, |q| q as f64 + 0.5)).collect();
    let m32 = M as u32;
    for p in 2..=r {
        let cp = cs[(p - 1) as usize];
        if cs[p as usize] == cp {
            continue;
        }
        let sp = ss[(p - 1) as usize];
        let p2 = p * p;
        let lim = r.min(n / p2) as usize;
        let split = (r / p).min(lim as u64) as usize;
        let pu = p as usize;
        let inv_p = 1.0 / p as f64;
        let cp64 = cp as u64;
        let spm = M as u32 - sp; // adding this is subtracting sp
        // SAFETY of the unchecked accesses below: i <= lim <= r, i * p <= r for i <= split, and
        // the divided values are floor(n / (i p)) <= r when i p > r; all arrays have r + 1 entries.
        for i in 1..=split {
            unsafe {
                *cl.get_unchecked_mut(i) -= *cl.get_unchecked(i * pu) - cp64;
                let t = ((*sl.get_unchecked(i * pu) + spm) as u64 * p % M) as u32;
                *sl.get_unchecked_mut(i) = sub_mod(*sl.get_unchecked(i), t, m32);
            }
        }
        for i in split + 1..=lim {
            unsafe {
                let d = (*vf.get_unchecked(i) * inv_p) as usize;
                debug_assert!(d <= ru);
                *cl.get_unchecked_mut(i) -= (*cs.get_unchecked(d) - cp) as u64;
                let t = ((*ss.get_unchecked(d) + spm) as u64 * p % M) as u32;
                *sl.get_unchecked_mut(i) = sub_mod(*sl.get_unchecked(i), t, m32);
            }
        }
        if p2 <= r {
            // v / p for v <= r < 2^32 by a 64-bit reciprocal (Lemire)
            let magic = u64::MAX / p + 1;
            for v in (p2..=r).rev() {
                let d = ((v as u128 * magic as u128) >> 64) as usize;
                let vu = v as usize;
                unsafe {
                    *cs.get_unchecked_mut(vu) -= *cs.get_unchecked(d) - cp;
                    let t = ((*ss.get_unchecked(d) + spm) as u64 * p % M) as u32;
                    *ss.get_unchecked_mut(vu) = sub_mod(*ss.get_unchecked(vu), t, m32);
                }
            }
        }
    }
    let mut count = Quotients::new(n, 0u64);
    count.small = cs.iter().map(|&x| x as u64).collect();
    count.large = cl;
    let mut sum = Quotients::new(n, 0u64);
    sum.small = ss.iter().map(|&x| x as u64).collect();
    sum.large = sl.iter().map(|&x| x as u64).collect();
    let _ = ru;
    (count, sum)
}

/// `(a - b) mod m` for `a, b < m`.
#[inline(always)]
fn sub_mod(a: u32, b: u32, m: u32) -> u32 {
    let (d, borrow) = a.overflowing_sub(b);
    if borrow {
        d.wrapping_add(m)
    } else {
        d
    }
}
