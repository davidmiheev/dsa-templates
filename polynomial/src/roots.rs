//! Roots of a polynomial in `F_M`: `gcd(f, x^M - x)` collects the product of the distinct linear
//! factors, and Cantor–Zassenhaus splitting (`gcd(h, (x + a)^((M-1)/2) - 1)` for random `a`) peels
//! them apart. `O(n^2 + n log n log M)` expected for degree `n`.

use crate::fps::{div_rem, mul, trim};
use crate::gcd::gcd_naive;
use modint::StaticModInt;

type Mint<const M: u32> = StaticModInt<M>;

/// `h * k mod g`.
fn mul_mod<const M: u32>(h: &[Mint<M>], k: &[Mint<M>], g: &[Mint<M>]) -> Vec<Mint<M>> {
    let (_, mut r) = div_rem(&mul(h, k), g);
    trim(&mut r);
    r
}

/// `base^e mod g` (`g` monic of degree >= 1).
fn pow_mod_poly<const M: u32>(base: &[Mint<M>], mut e: u64, g: &[Mint<M>]) -> Vec<Mint<M>> {
    let (_, mut b) = div_rem(base, g);
    trim(&mut b);
    let mut result = vec![Mint::<M>::one()];
    while e > 0 {
        if e & 1 == 1 {
            result = mul_mod(&result, &b, g);
        }
        e >>= 1;
        if e > 0 {
            b = mul_mod(&b, &b, g);
        }
    }
    result
}

fn monic<const M: u32>(f: &[Mint<M>]) -> Vec<Mint<M>> {
    let inv = f.last().expect("non-zero polynomial").inv();
    f.iter().map(|&x| x * inv).collect()
}

/// All **distinct** roots of `f` in `F_M` (`M` prime), in increasing order.
///
/// # Panics
/// Panics if `f` is the zero polynomial.
pub fn find_roots<const M: u32>(f: &[Mint<M>]) -> Vec<Mint<M>> {
    let mut f = f.to_vec();
    trim(&mut f);
    assert!(!f.is_empty(), "the zero polynomial has every element as a root");
    let mut roots = Vec::new();
    // factor out x
    let zeros = f.iter().position(|x| x.val() != 0).unwrap();
    if zeros > 0 {
        roots.push(Mint::zero());
        f.drain(..zeros);
    }
    if f.len() > 1 {
        let g = monic(&f);
        // x^M mod g
        let x = vec![Mint::<M>::zero(), Mint::one()];
        let xp = pow_mod_poly(&x, M as u64, &g);
        let mut diff = xp;
        diff.resize(diff.len().max(2), Mint::zero());
        diff[1] -= Mint::one();
        trim(&mut diff);
        let h = if diff.is_empty() { g.clone() } else { gcd_naive(&g, &diff) };
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        split(&h, &mut roots, &mut state);
    }
    roots.sort_by_key(|r| r.val());
    roots
}

fn split<const M: u32>(h: &[Mint<M>], roots: &mut Vec<Mint<M>>, state: &mut u64) {
    // h is monic, squarefree and a product of distinct linear factors
    match h.len() {
        0 | 1 => {}
        2 => roots.push(-h[0]),
        _ => loop {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            let a = Mint::<M>::new(*state);
            let base = vec![a, Mint::one()];
            let mut pw = pow_mod_poly(&base, (M as u64 - 1) / 2, h);
            pw.resize(pw.len().max(1), Mint::zero());
            pw[0] -= Mint::one();
            trim(&mut pw);
            let d = if pw.is_empty() { h.to_vec() } else { gcd_naive(h, &pw) };
            if d.len() > 1 && d.len() < h.len() {
                let (q, _) = div_rem(h, &d);
                let q = monic(&q);
                split(&d, roots, state);
                split(&q, roots, state);
                return;
            }
        },
    }
}
