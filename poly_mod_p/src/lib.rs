//! Polynomials over `F_p` for a runtime prime `p < 2^31`, as `Vec<u64>` (index = degree, trimmed),
//! with schoolbook arithmetic — meant for moderate degrees — and **factorisation into monic
//! irreducibles**:
//!
//! 1. *square-free decomposition* (Yun's algorithm; a derivative that vanishes means the polynomial
//!    is a `p`-th power, whose root has the coefficients at multiples of `p`),
//! 2. *distinct-degree factorisation* with `gcd(f, x^(p^d) - x)`,
//! 3. *equal-degree splitting* (Cantor–Zassenhaus; for `p = 2` with the trace map).
//!
//! # Example
//!
//! ```
//! use poly_mod_p::factorize;
//!
//! // x^4 + 1 = (x + 1)^4 over F_2
//! assert_eq!(factorize(&[1, 0, 0, 0, 1], 2), vec![(vec![1, 1], 4)]);
//! // x^2 + 1 over F_5 = (x + 2)(x + 3)
//! let f = factorize(&[1, 0, 1], 5);
//! assert_eq!(f, vec![(vec![2, 1], 1), (vec![3, 1], 1)]);
//! ```

pub type Poly = Vec<u64>;

/// Remove trailing zeros.
pub fn trim(a: &mut Poly) {
    while a.last() == Some(&0) {
        a.pop();
    }
}

fn pow_mod(mut b: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1 % p;
    b %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % p;
        }
        b = b * b % p;
        e >>= 1;
    }
    r
}

fn inv(a: u64, p: u64) -> u64 {
    pow_mod(a, p - 2, p)
}

fn sub(a: &[u64], b: &[u64], p: u64) -> Poly {
    let mut r = vec![0u64; a.len().max(b.len())];
    for (i, &x) in a.iter().enumerate() {
        r[i] = x;
    }
    for (i, &x) in b.iter().enumerate() {
        r[i] = (r[i] + p - x) % p;
    }
    trim(&mut r);
    r
}

fn add(a: &[u64], b: &[u64], p: u64) -> Poly {
    let mut r = vec![0u64; a.len().max(b.len())];
    for (i, &x) in a.iter().enumerate() {
        r[i] = x;
    }
    for (i, &x) in b.iter().enumerate() {
        r[i] = (r[i] + x) % p;
    }
    trim(&mut r);
    r
}

fn mul(a: &[u64], b: &[u64], p: u64) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = (r[i + j] + x * y) % p;
        }
    }
    trim(&mut r);
    r
}

/// Quotient and remainder of `a / b` (`b != 0`).
pub fn div_rem(a: &[u64], b: &[u64], p: u64) -> (Poly, Poly) {
    assert!(!b.is_empty(), "division by the zero polynomial");
    let mut r = a.to_vec();
    trim(&mut r);
    if r.len() < b.len() {
        return (Vec::new(), r);
    }
    let lead_inv = inv(*b.last().unwrap(), p);
    let mut q = vec![0u64; r.len() - b.len() + 1];
    while r.len() >= b.len() {
        let c = *r.last().unwrap() * lead_inv % p;
        let shift = r.len() - b.len();
        q[shift] = c;
        for (i, &y) in b.iter().enumerate() {
            r[shift + i] = (r[shift + i] + p - c * y % p) % p;
        }
        trim(&mut r);
        if r.is_empty() {
            break;
        }
    }
    trim(&mut q);
    (q, r)
}

fn monic(a: &[u64], p: u64) -> Poly {
    let li = inv(*a.last().expect("non-zero polynomial"), p);
    a.iter().map(|&x| x * li % p).collect()
}

/// Monic greatest common divisor.
pub fn gcd(a: &[u64], b: &[u64], p: u64) -> Poly {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    trim(&mut a);
    trim(&mut b);
    while !b.is_empty() {
        let (_, r) = div_rem(&a, &b, p);
        a = std::mem::replace(&mut b, r);
    }
    if a.is_empty() { a } else { monic(&a, p) }
}

fn derivative(a: &[u64], p: u64) -> Poly {
    let mut r: Poly = (1..a.len()).map(|i| a[i] * (i as u64 % p) % p).collect();
    trim(&mut r);
    r
}

/// `base^e mod m` for polynomials.
fn pow_poly(base: &[u64], mut e: u64, m: &[u64], p: u64) -> Poly {
    let (_, mut b) = div_rem(base, m, p);
    let mut r: Poly = vec![1];
    while e > 0 {
        if e & 1 == 1 {
            r = div_rem(&mul(&r, &b, p), m, p).1;
        }
        e >>= 1;
        if e > 0 {
            b = div_rem(&mul(&b, &b, p), m, p).1;
        }
    }
    r
}

fn deg(a: &[u64]) -> usize {
    a.len().saturating_sub(1)
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

/// Square-free decomposition of a monic `f`: pairs `(g, e)` of pairwise coprime monic square-free
/// polynomials with `f = prod g^e`.
fn square_free(f: &[u64], p: u64) -> Vec<(Poly, usize)> {
    let mut res = Vec::new();
    if f.len() <= 1 {
        return res;
    }
    let d = derivative(f, p);
    if d.is_empty() {
        return p_th_root_parts(f, p);
    }
    let mut c = gcd(f, &d, p);
    let mut w = div_rem(f, &c, p).0;
    let mut i = 1;
    while w.len() > 1 {
        let y = gcd(&w, &c, p);
        let z = div_rem(&w, &y, p).0;
        if z.len() > 1 {
            res.push((z, i));
        }
        c = div_rem(&c, &y, p).0;
        w = y;
        i += 1;
    }
    if c.len() > 1 {
        res.extend(p_th_root_parts(&c, p));
    }
    res
}

/// `f = g(x)^p` with `g` read off the coefficients at multiples of `p`; factor `g`.
fn p_th_root_parts(f: &[u64], p: u64) -> Vec<(Poly, usize)> {
    let g: Poly = f.iter().step_by(p as usize).copied().collect();
    square_free(&monic(&g, p), p).into_iter().map(|(h, e)| (h, e * p as usize)).collect()
}

/// Distinct-degree factorisation of a monic square-free `f`: `(product of the irreducible factors of
/// degree d, d)`.
fn distinct_degree(f: &[u64], p: u64) -> Vec<(Poly, usize)> {
    let mut res = Vec::new();
    let mut rest = f.to_vec();
    let x: Poly = vec![0, 1];
    let mut h = x.clone();
    let mut d = 1;
    while deg(&rest) >= 2 * d {
        h = pow_poly(&h, p, &rest, p);
        let g = gcd(&rest, &sub(&h, &x, p), p);
        if g.len() > 1 {
            rest = div_rem(&rest, &g, p).0;
            h = div_rem(&h, &rest, p).1;
            res.push((g, d));
        }
        d += 1;
    }
    if rest.len() > 1 {
        let dg = deg(&rest);
        res.push((rest, dg));
    }
    res
}

/// Split `f`, a monic product of distinct irreducibles of degree `d`, into those irreducibles.
fn equal_degree(f: &[u64], d: usize, p: u64, rng: &mut Rng, out: &mut Vec<Poly>) {
    if deg(f) == d {
        out.push(f.to_vec());
        return;
    }
    loop {
        let r: Poly = (0..deg(f)).map(|_| rng.next() % p).collect();
        let mut r = r;
        trim(&mut r);
        if r.len() < 2 {
            continue;
        }
        let s = if p == 2 {
            // trace map: r + r^2 + ... + r^(2^(d-1))
            let (mut t, mut cur) = (r.clone(), r.clone());
            for _ in 1..d {
                cur = div_rem(&mul(&cur, &cur, p), f, p).1;
                t = add(&t, &cur, p);
            }
            gcd(f, &t, p)
        } else {
            // r^((p^d - 1) / 2) = (r^(1 + p + ... + p^(d-1)))^((p-1)/2)
            let mut a = div_rem(&r, f, p).1;
            for _ in 1..d {
                a = div_rem(&mul(&pow_poly(&a, p, f, p), &r, p), f, p).1;
            }
            let t = pow_poly(&a, (p - 1) / 2, f, p);
            gcd(f, &sub(&t, &[1], p), p)
        };
        if s.len() > 1 && s.len() < f.len() {
            let q = div_rem(f, &s, p).0;
            equal_degree(&s, d, p, rng, out);
            equal_degree(&monic(&q, p), d, p, rng, out);
            return;
        }
    }
}

/// Factorise the **monic** polynomial `f` over `F_p` (`p` prime) into monic irreducible factors with
/// multiplicities, sorted by (degree, coefficients). `f = 1` gives an empty list.
pub fn factorize(f: &[u64], p: u64) -> Vec<(Poly, usize)> {
    let mut f = f.to_vec();
    trim(&mut f);
    assert!(!f.is_empty(), "the zero polynomial cannot be factorised");
    assert_eq!(*f.last().unwrap(), 1, "the polynomial must be monic");
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    let mut result: Vec<(Poly, usize)> = Vec::new();
    for (g, e) in square_free(&f, p) {
        for (h, d) in distinct_degree(&g, p) {
            let mut parts = Vec::new();
            equal_degree(&h, d, p, &mut rng, &mut parts);
            for q in parts {
                result.push((q, e));
            }
        }
    }
    result.sort_by(|a, b| (a.0.len(), &a.0).cmp(&(b.0.len(), &b.0)));
    result
}
