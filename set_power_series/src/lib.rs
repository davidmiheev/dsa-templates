//! Functions on the subsets of an `n`-element set, stored as arrays of length `2^n` (index = bit
//! mask) with entries modulo a prime `M < 2^31` (`u32` slices).
//!
//! * [`zeta_subset`] / [`mobius_subset`]: `b[S] = sum_{T subset S} a[T]` and its inverse;
//!   [`zeta_superset`] / [`mobius_superset`] for supersets. `O(n 2^n)`.
//! * [`walsh_hadamard`]: the XOR transform (self-inverse up to the factor `2^n`).
//! * [`and_convolution`], [`or_convolution`], [`xor_convolution`]: `c[k] = sum a[i] b[j]` over the
//!   pairs with `i & j = k` (resp. `|`, `^`).
//! * [`subset_convolution`]: `c[S] = sum_{T subset S} a[T] b[S \ T]`, by ranking the sets with
//!   their sizes: `O(n^2 2^n)`.
//!
//! # Example
//!
//! ```
//! use set_power_series::{subset_convolution, xor_convolution};
//!
//! const M: u32 = 998_244_353;
//! assert_eq!(xor_convolution::<M>(&[1, 2, 3, 4], &[5, 6, 7, 8]), vec![70, 68, 62, 60]);
//! // (1 + x_0)(1 + x_0) with x_0^2 = 0 is 1 + 2 x_0
//! assert_eq!(subset_convolution::<M>(&[1, 1], &[1, 1]), vec![1, 2]);
//! ```

fn pow_mod(mut b: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1;
    b %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % m;
        }
        b = b * b % m;
        e >>= 1;
    }
    r
}

#[inline(always)]
fn add<const M: u32>(a: u32, b: u32) -> u32 {
    let s = a + b;
    if s >= M {
        s - M
    } else {
        s
    }
}

#[inline(always)]
fn sub<const M: u32>(a: u32, b: u32) -> u32 {
    if a >= b {
        a - b
    } else {
        a + M - b
    }
}

fn check_len(len: usize) -> usize {
    assert!(len.is_power_of_two(), "length {len} is not a power of two");
    len
}

/// `a[S] <- sum_{T subset S} a[T]`.
pub fn zeta_subset<const M: u32>(a: &mut [u32]) {
    let n = check_len(a.len());
    let mut h = 1;
    while h < n {
        for block in a.chunks_exact_mut(2 * h) {
            let (lo, hi) = block.split_at_mut(h);
            for (x, &y) in hi.iter_mut().zip(lo.iter()) {
                *x = add::<M>(*x, y);
            }
        }
        h *= 2;
    }
}

/// Inverse of [`zeta_subset`].
pub fn mobius_subset<const M: u32>(a: &mut [u32]) {
    let n = check_len(a.len());
    let mut h = 1;
    while h < n {
        for block in a.chunks_exact_mut(2 * h) {
            let (lo, hi) = block.split_at_mut(h);
            for (x, &y) in hi.iter_mut().zip(lo.iter()) {
                *x = sub::<M>(*x, y);
            }
        }
        h *= 2;
    }
}

/// `a[S] <- sum_{T superset S} a[T]`.
pub fn zeta_superset<const M: u32>(a: &mut [u32]) {
    let n = check_len(a.len());
    let mut h = 1;
    while h < n {
        for block in a.chunks_exact_mut(2 * h) {
            let (lo, hi) = block.split_at_mut(h);
            for (x, &y) in lo.iter_mut().zip(hi.iter()) {
                *x = add::<M>(*x, y);
            }
        }
        h *= 2;
    }
}

/// Inverse of [`zeta_superset`].
pub fn mobius_superset<const M: u32>(a: &mut [u32]) {
    let n = check_len(a.len());
    let mut h = 1;
    while h < n {
        for block in a.chunks_exact_mut(2 * h) {
            let (lo, hi) = block.split_at_mut(h);
            for (x, &y) in lo.iter_mut().zip(hi.iter()) {
                *x = sub::<M>(*x, y);
            }
        }
        h *= 2;
    }
}

/// The Walsh–Hadamard transform `a[S] <- sum_T (-1)^{|S & T|} a[T]`; applying it twice multiplies
/// by `len`.
pub fn walsh_hadamard<const M: u32>(a: &mut [u32]) {
    let n = check_len(a.len());
    let mut h = 1;
    while h < n {
        for block in a.chunks_exact_mut(2 * h) {
            let (lo, hi) = block.split_at_mut(h);
            for (x, y) in lo.iter_mut().zip(hi.iter_mut()) {
                let (u, v) = (*x, *y);
                *x = add::<M>(u, v);
                *y = sub::<M>(u, v);
            }
        }
        h *= 2;
    }
}

fn pointwise<const M: u32>(a: &mut [u32], b: &[u32]) {
    for (x, &y) in a.iter_mut().zip(b) {
        *x = (*x as u64 * y as u64 % M as u64) as u32;
    }
}

/// `c[k] = sum_{i & j = k} a[i] b[j]`. Both slices need the same power-of-two length.
pub fn and_convolution<const M: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len(), "length mismatch");
    let (mut fa, mut fb) = (a.to_vec(), b.to_vec());
    zeta_superset::<M>(&mut fa);
    zeta_superset::<M>(&mut fb);
    pointwise::<M>(&mut fa, &fb);
    mobius_superset::<M>(&mut fa);
    fa
}

/// `c[k] = sum_{i | j = k} a[i] b[j]`.
pub fn or_convolution<const M: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len(), "length mismatch");
    let (mut fa, mut fb) = (a.to_vec(), b.to_vec());
    zeta_subset::<M>(&mut fa);
    zeta_subset::<M>(&mut fb);
    pointwise::<M>(&mut fa, &fb);
    mobius_subset::<M>(&mut fa);
    fa
}

/// `c[k] = sum_{i ^ j = k} a[i] b[j]` (requires `M` odd).
pub fn xor_convolution<const M: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len(), "length mismatch");
    assert!(M % 2 == 1, "the modulus must be odd");
    let (mut fa, mut fb) = (a.to_vec(), b.to_vec());
    walsh_hadamard::<M>(&mut fa);
    walsh_hadamard::<M>(&mut fb);
    pointwise::<M>(&mut fa, &fb);
    walsh_hadamard::<M>(&mut fa);
    let inv = pow_mod(a.len() as u64, M as u64 - 2, M as u64);
    for x in fa.iter_mut() {
        *x = (*x as u64 * inv % M as u64) as u32;
    }
    fa
}

/// `c[S] = sum_{T subset S} a[T] b[S \ T]`. Both slices need the same power-of-two length.
pub fn subset_convolution<const M: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len(), "length mismatch");
    let len = check_len(a.len());
    let n = len.trailing_zeros() as usize;
    if n == 0 {
        return vec![(a[0] as u64 * b[0] as u64 % M as u64) as u32];
    }
    // ranked[r][S] = a[S] if |S| = r else 0, then zeta transformed
    let rank = |v: &[u32]| -> Vec<Vec<u32>> {
        let mut out = vec![vec![0u32; len]; n + 1];
        for (s, &x) in v.iter().enumerate() {
            out[s.count_ones() as usize][s] = x;
        }
        for layer in out.iter_mut() {
            zeta_subset::<M>(layer);
        }
        out
    };
    let (ra, rb) = (rank(a), rank(b));
    let mut c = vec![0u32; len];
    let lim = 8 * (M as u64) * (M as u64);
    let mut acc = vec![0u64; len];
    let mut layer = vec![0u32; len];
    for k in 0..=n {
        acc.iter_mut().for_each(|x| *x = 0);
        for i in 0..=k {
            let (x, y) = (&ra[i], &rb[k - i]);
            for ((t, &u), &v) in acc.iter_mut().zip(x).zip(y) {
                *t += u as u64 * v as u64;
                if *t >= lim {
                    *t -= lim;
                }
            }
        }
        for (l, &t) in layer.iter_mut().zip(&acc) {
            *l = (t % M as u64) as u32;
        }
        mobius_subset::<M>(&mut layer);
        for (s, &v) in layer.iter().enumerate() {
            if s.count_ones() as usize == k {
                c[s] = v;
            }
        }
    }
    c
}
