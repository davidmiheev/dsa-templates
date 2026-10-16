//! Polynomial greatest common divisors: the schoolbook Euclidean algorithm for small degrees and
//! the **half-gcd** (`O(n log^2 n)`) for large ones, with the extended version giving modular
//! inverses of polynomials.
//!
//! The Euclidean remainder sequence is `r_0 = a`, `r_1 = b`, `r_{i+1} = r_{i-1} - q_i r_i`; the
//! matrix `M_h = Q_h ... Q_1` with `Q_i = [[0, 1], [1, -q_i]]` maps `(a, b)` to `(r_h, r_{h+1})`.
//! [`half_gcd`] finds `M_h` for the first `h` with `deg r_h >= ceil(deg a / 2) > deg r_{h+1}`
//! recursively: the quotients of the top halves of `a` and `b` are already the true quotients up
//! to a suitable degree, so a half-sized call, one division step and a second half-sized call
//! on shifted operands advance the sequence by half its length.

use crate::fps::{div_rem, mul, trim};
use modint::StaticModInt;

type Mint<const M: u32> = StaticModInt<M>;
type Poly<const M: u32> = Vec<Mint<M>>;
type Matrix<const M: u32> = [[Poly<M>; 2]; 2];

/// Degree of a trimmed polynomial; `-1` stands for the zero polynomial.
fn degree<const M: u32>(a: &[Mint<M>]) -> isize {
    let mut n = a.len();
    while n > 0 && a[n - 1].val() == 0 {
        n -= 1;
    }
    n as isize - 1
}

fn add_poly<const M: u32>(a: &[Mint<M>], b: &[Mint<M>]) -> Poly<M> {
    let mut r = vec![Mint::<M>::zero(); a.len().max(b.len())];
    for (i, &x) in a.iter().enumerate() {
        r[i] += x;
    }
    for (i, &x) in b.iter().enumerate() {
        r[i] += x;
    }
    trim(&mut r);
    r
}

fn sub_poly<const M: u32>(a: &[Mint<M>], b: &[Mint<M>]) -> Poly<M> {
    let mut r = vec![Mint::<M>::zero(); a.len().max(b.len())];
    for (i, &x) in a.iter().enumerate() {
        r[i] += x;
    }
    for (i, &x) in b.iter().enumerate() {
        r[i] -= x;
    }
    trim(&mut r);
    r
}

fn shifted<const M: u32>(a: &[Mint<M>], k: usize) -> Poly<M> {
    if a.len() <= k { Vec::new() } else { a[k..].to_vec() }
}

fn identity<const M: u32>() -> Matrix<M> {
    [[vec![Mint::one()], Vec::new()], [Vec::new(), vec![Mint::one()]]]
}

fn matmul<const M: u32>(a: &Matrix<M>, b: &Matrix<M>) -> Matrix<M> {
    let cell = |i: usize, j: usize| -> Poly<M> {
        let mut t = add_poly(&mul(&a[i][0], &b[0][j]), &mul(&a[i][1], &b[1][j]));
        trim(&mut t);
        t
    };
    [[cell(0, 0), cell(0, 1)], [cell(1, 0), cell(1, 1)]]
}

fn apply<const M: u32>(m: &Matrix<M>, a: &[Mint<M>], b: &[Mint<M>]) -> (Poly<M>, Poly<M>) {
    (add_poly(&mul(&m[0][0], a), &mul(&m[0][1], b)), add_poly(&mul(&m[1][0], a), &mul(&m[1][1], b)))
}

/// One Euclid step matrix `[[0, 1], [1, -q]]`.
fn step_matrix<const M: u32>(q: &[Mint<M>]) -> Matrix<M> {
    let nq: Poly<M> = q.iter().map(|&x| -x).collect();
    [[Vec::new(), vec![Mint::one()]], [vec![Mint::one()], nq]]
}

const NAIVE: usize = 64;

fn is_identity<const M: u32>(m: &Matrix<M>) -> bool {
    m[0][0].len() == 1 && m[0][0][0] == Mint::one() && m[0][1].is_empty() && m[1][0].is_empty() && m[1][1].len() == 1 && m[1][1][0] == Mint::one()
}

/// The matrix `M_h` described in the module documentation (`deg a > deg b` unless `b = 0`).
pub fn half_gcd<const M: u32>(a: &[Mint<M>], b: &[Mint<M>]) -> Matrix<M> {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    trim(&mut a);
    trim(&mut b);
    let n = degree(&a);
    assert!(n > degree(&b), "half_gcd needs deg a > deg b");
    let m = ((n + 1) / 2) as isize;
    if degree(&b) < m {
        return identity();
    }
    if n as usize <= NAIVE {
        // schoolbook: run the Euclidean algorithm until the degree drops below m
        let mut mat = identity();
        let (mut c, mut d) = (a, b);
        while degree(&d) >= m {
            let (q, r) = div_rem(&c, &d);
            mat = matmul(&step_matrix(&q), &mat);
            c = std::mem::take(&mut d);
            d = r;
            trim(&mut d);
        }
        return mat;
    }
    let k = (n / 2) as usize;
    let r1 = half_gcd(&shifted(&a, k), &shifted(&b, k));
    let (c, d) = apply(&r1, &a, &b);
    if degree(&d) < m {
        return r1;
    }
    let (q, r) = div_rem(&c, &d);
    let (c2, d2) = (d, r);
    let q1 = step_matrix(&q);
    let r1q = matmul(&q1, &r1);
    if degree(&d2) < m {
        return r1q;
    }
    let n2 = degree(&c2);
    let s = (2 * m - n2) as usize;
    let r2 = half_gcd(&shifted(&c2, s), &shifted(&d2, s));
    matmul(&r2, &r1q)
}

/// Extended Euclid: `(g, s, t)` with `g = gcd(a, b)` (made monic), `g = s a + t b`.
pub fn ext_gcd<const M: u32>(a: &[Mint<M>], b: &[Mint<M>]) -> (Poly<M>, Poly<M>, Poly<M>) {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    trim(&mut a);
    trim(&mut b);
    let mut total = identity::<M>();
    if degree(&a) < degree(&b) {
        std::mem::swap(&mut a, &mut b);
        total = [[Vec::new(), vec![Mint::one()]], [vec![Mint::one()], Vec::new()]];
    }
    while !b.is_empty() {
        let mut mat = if degree(&a) as usize > NAIVE && degree(&a) > degree(&b) { half_gcd(&a, &b) } else { identity() };
        if is_identity(&mat) {
            // the half-gcd made no progress (or the operands are small): one division step
            let (q, _) = div_rem(&a, &b);
            mat = step_matrix(&q);
        }
        let (c, d) = apply(&mat, &a, &b);
        total = matmul(&mat, &total);
        a = c;
        b = d;
        trim(&mut a);
        trim(&mut b);
    }
    // total * (a0, b0) = (g, 0): the first row gives the cofactors
    let lead = *a.last().expect("gcd of two zero polynomials");
    let li = lead.inv();
    let scale = |p: &[Mint<M>]| -> Poly<M> { p.iter().map(|&x| x * li).collect() };
    (scale(&a), scale(&total[0][0]), scale(&total[0][1]))
}

/// The polynomial `h` of degree `< deg g` with `f h = 1 mod g`; an **empty vector** when no inverse
/// exists (`gcd(f, g) != 1`). For a constant `g` the inverse is the zero polynomial, returned as
/// `[0]`. `g` must be non-zero.
pub fn inverse_mod<const M: u32>(f: &[Mint<M>], g: &[Mint<M>]) -> Vec<Mint<M>> {
    let mut gt = g.to_vec();
    trim(&mut gt);
    assert!(!gt.is_empty(), "modulus must be non-zero");
    if gt.len() == 1 {
        return vec![Mint::zero()];
    }
    let (_, mut f_red) = div_rem(f, &gt);
    trim(&mut f_red);
    if f_red.is_empty() {
        return Vec::new();
    }
    let (d, _s, t) = ext_gcd(&gt, &f_red);
    // g = d is the gcd: d = s g + t f  =>  t / d is the inverse when d is a constant
    if degree(&d) != 0 {
        return Vec::new();
    }
    let (_, mut h) = div_rem(&t, &gt);
    h.resize(gt.len() - 1, Mint::zero());
    h.iter_mut().for_each(|x| *x /= d[0]);
    h
}

/// Euclid's algorithm on plain schoolbook arithmetic (monic result), for small degrees.
pub fn gcd_naive<const M: u32>(a: &[Mint<M>], b: &[Mint<M>]) -> Poly<M> {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    trim(&mut a);
    trim(&mut b);
    while !b.is_empty() {
        let lead_inv = b.last().unwrap().inv();
        let db = b.len();
        while a.len() >= db {
            let c = *a.last().unwrap() * lead_inv;
            let shift = a.len() - db;
            for i in 0..db {
                let t = c * b[i];
                a[shift + i] -= t;
            }
            trim(&mut a);
            if a.is_empty() {
                break;
            }
        }
        std::mem::swap(&mut a, &mut b);
    }
    if let Some(&l) = a.last() {
        let li = l.inv();
        a.iter_mut().for_each(|x| *x *= li);
    }
    a
}
