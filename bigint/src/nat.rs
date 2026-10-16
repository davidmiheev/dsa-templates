//! Natural numbers as little-endian limb vectors in base `B <= 2^20` (no leading zero limbs, zero
//! is the empty vector).

use convolution::{convolution_mod, convolution_u64, P1, P3};
use std::cmp::Ordering;

pub fn trim(v: &mut Vec<u32>) {
    while v.last() == Some(&0) {
        v.pop();
    }
}

pub fn cmp(a: &[u32], b: &[u32]) -> Ordering {
    a.len().cmp(&b.len()).then_with(|| a.iter().rev().cmp(b.iter().rev()))
}

pub fn add<const B: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    let (a, b) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    let mut out = Vec::with_capacity(a.len() + 1);
    let mut carry = 0u32;
    for i in 0..a.len() {
        let mut s = a[i] + carry + if i < b.len() { b[i] } else { 0 };
        carry = (s >= B) as u32;
        if carry == 1 {
            s -= B;
        }
        out.push(s);
    }
    if carry == 1 {
        out.push(1);
    }
    out
}

/// `a -= b` for `a >= b`.
pub fn sub_in_place<const B: u32>(a: &mut Vec<u32>, b: &[u32]) {
    debug_assert!(cmp(a, b) != Ordering::Less, "subtraction underflow");
    let mut borrow = 0u32;
    for i in 0..a.len() {
        let t = if i < b.len() { b[i] } else { 0 } + borrow;
        if a[i] >= t {
            a[i] -= t;
            borrow = 0;
        } else {
            a[i] = a[i] + B - t;
            borrow = 1;
        }
        if i >= b.len() && borrow == 0 {
            break;
        }
    }
    trim(a);
}

pub fn sub<const B: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut r = a.to_vec();
    sub_in_place::<B>(&mut r, b);
    r
}

pub fn add_small<const B: u32>(a: &mut Vec<u32>, mut s: u32) {
    for x in a.iter_mut() {
        let t = *x + s;
        if t >= B {
            *x = t - B;
            s = 1;
        } else {
            *x = t;
            return;
        }
    }
    if s > 0 {
        a.push(s);
    }
}

pub fn mul_small<const B: u32>(a: &[u32], s: u32) -> Vec<u32> {
    let mut out = Vec::with_capacity(a.len() + 1);
    let mut carry = 0u64;
    for &x in a {
        let t = x as u64 * s as u64 + carry;
        out.push((t % B as u64) as u32);
        carry = t / B as u64;
    }
    while carry > 0 {
        out.push((carry % B as u64) as u32);
        carry /= B as u64;
    }
    trim(&mut out);
    out
}

pub fn divrem_small<const B: u32>(a: &[u32], s: u32) -> (Vec<u32>, u32) {
    let mut q = vec![0u32; a.len()];
    let mut rem = 0u64;
    for i in (0..a.len()).rev() {
        let cur = rem * B as u64 + a[i] as u64;
        q[i] = (cur / s as u64) as u32;
        rem = cur % s as u64;
    }
    trim(&mut q);
    (q, rem as u32)
}

fn carry_out<const B: u32>(mut acc: Vec<u64>) -> Vec<u32> {
    let mut out = Vec::with_capacity(acc.len() + 2);
    let mut carry = 0u64;
    for x in acc.drain(..) {
        let t = x + carry;
        out.push((t % B as u64) as u32);
        carry = t / B as u64;
    }
    while carry > 0 {
        out.push((carry % B as u64) as u32);
        carry /= B as u64;
    }
    trim(&mut out);
    out
}

/// The product; schoolbook below 32 limbs, otherwise exact via two NTT primes (three when the
/// coefficients may exceed their product).
pub fn mul<const B: u32>(a: &[u32], b: &[u32]) -> Vec<u32> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let short = a.len().min(b.len());
    if short <= 32 {
        let mut acc = vec![0u64; a.len() + b.len()];
        for (i, &x) in a.iter().enumerate() {
            if x == 0 {
                continue;
            }
            for (j, &y) in b.iter().enumerate() {
                acc[i + j] += x as u64 * y as u64;
            }
        }
        return carry_out::<B>(acc);
    }
    let bound = short as u128 * ((B - 1) as u128).pow(2);
    if bound >= P1 as u128 * P3 as u128 {
        assert!(bound < 1 << 64, "limbs too large for the exact convolution");
        let widen = |v: &[u32]| -> Vec<u64> { v.iter().map(|&x| x as u64).collect() };
        return carry_out::<B>(convolution_u64(&widen(a), &widen(b)));
    }
    let c1 = convolution_mod::<P1>(a, b);
    let c3 = convolution_mod::<P3>(a, b);
    let (m1, m3) = (P1 as u64, P3 as u64);
    // 1 / P1 modulo P3
    let inv = {
        let (mut r, mut e, mut base) = (1u64, m3 - 2, m1 % m3);
        while e > 0 {
            if e & 1 == 1 {
                r = r * base % m3;
            }
            base = base * base % m3;
            e >>= 1;
        }
        r
    };
    let acc: Vec<u64> = c1
        .iter()
        .zip(&c3)
        .map(|(&r1, &r3)| {
            let (r1, r3) = (r1 as u64, r3 as u64);
            let t = (r3 + m3 - r1 % m3) % m3 * inv % m3;
            r1 + m1 * t
        })
        .collect();
    carry_out::<B>(acc)
}

/// `a * B^k`.
pub fn shl(a: &[u32], k: usize) -> Vec<u32> {
    if a.is_empty() {
        return Vec::new();
    }
    let mut v = vec![0u32; k];
    v.extend_from_slice(a);
    v
}

/// `floor(a / B^k)`.
pub fn shr(a: &[u32], k: usize) -> &[u32] {
    if k >= a.len() {
        &[]
    } else {
        &a[k..]
    }
}

const SCHOOL: usize = 40;

/// Knuth's algorithm D, `O(|q| |b|)`; `b` must have at least two limbs.
fn knuth<const B: u32>(a: &[u32], b: &[u32]) -> (Vec<u32>, Vec<u32>) {
    let n = b.len();
    assert!(n >= 2 && b[n - 1] != 0, "divisor needs two limbs");
    if cmp(a, b) == Ordering::Less {
        return (Vec::new(), a.to_vec());
    }
    let base = B as u64;
    let d = B / (b[n - 1] + 1);
    let bn = mul_small::<B>(b, d);
    let mut an = mul_small::<B>(a, d);
    an.resize(a.len() + 1, 0);
    let m = a.len() - n;
    let mut q = vec![0u32; m + 1];
    let (b1, b2) = (bn[n - 1] as u64, bn[n - 2] as u64);
    for j in (0..=m).rev() {
        let num = an[j + n] as u64 * base + an[j + n - 1] as u64;
        let (mut qhat, mut rhat) = (num / b1, num % b1);
        while qhat >= base || qhat * b2 > rhat * base + an[j + n - 2] as u64 {
            qhat -= 1;
            rhat += b1;
            if rhat >= base {
                break;
            }
        }
        let (mut borrow, mut carry) = (0i64, 0u64);
        for i in 0..n {
            let p = qhat * bn[i] as u64 + carry;
            carry = p / base;
            let s = an[i + j] as i64 - (p % base) as i64 - borrow;
            if s < 0 {
                an[i + j] = (s + base as i64) as u32;
                borrow = 1;
            } else {
                an[i + j] = s as u32;
                borrow = 0;
            }
        }
        let s = an[j + n] as i64 - carry as i64 - borrow;
        if s < 0 {
            an[j + n] = (s + base as i64) as u32;
            qhat -= 1;
            let mut c = 0u32;
            for i in 0..n {
                let t = an[i + j] + bn[i] + c;
                c = (t >= B) as u32;
                an[i + j] = if c == 1 { t - B } else { t };
            }
            an[j + n] = (an[j + n] + c) % B;
        } else {
            an[j + n] = s as u32;
        }
        q[j] = qhat as u32;
    }
    trim(&mut q);
    an.truncate(n);
    trim(&mut an);
    let (r, rem) = divrem_small::<B>(&an, d);
    debug_assert_eq!(rem, 0);
    (q, r)
}

/// `~ floor(B^{2m} / b)` (within a few units) for a normalised `b` of `m` limbs (top limb `>= B / 2`),
/// by Newton iteration doubling the precision.
fn reciprocal<const B: u32>(b: &[u32]) -> Vec<u32> {
    let m = b.len();
    if m <= SCHOOL {
        let mut num = vec![0u32; 2 * m];
        num.push(1);
        return knuth::<B>(&num, b).0;
    }
    let h = m / 2 + 2;
    let y = shl(&reciprocal::<B>(&b[m - h..]), m - h);
    let t = mul::<B>(b, &y);
    let mut one = vec![0u32; 2 * m];
    one.push(1);
    let (positive, e) = if cmp(&t, &one) != Ordering::Greater { (true, sub::<B>(&one, &t)) } else { (false, sub::<B>(&t, &one)) };
    let corr = mul::<B>(&y, &e);
    let corr = shr(&corr, 2 * m);
    let mut y = y;
    if positive {
        y = add::<B>(&y, corr);
    } else {
        sub_in_place::<B>(&mut y, corr);
    }
    y
}

/// An estimate of `floor(a / b)` for `a < B^{2m}` from the reciprocal `r` of a normalised `b`.
fn estimate<const B: u32>(a: &[u32], b: &[u32], r: &[u32]) -> Vec<u32> {
    let m = b.len();
    let top = shr(a, m - 1);
    let mut prod = mul::<B>(top, r);
    if prod.len() <= m + 1 {
        return Vec::new();
    }
    prod.drain(..m + 1);
    prod
}

/// Turn an estimate `q` of `a / b` (off by a few units) into the exact `(q, a - q b)`.
fn correct<const B: u32>(a: &[u32], b: &[u32], mut q: Vec<u32>) -> (Vec<u32>, Vec<u32>) {
    let mut prod = mul::<B>(&q, b);
    while cmp(&prod, a) == Ordering::Greater {
        sub_in_place::<B>(&mut prod, b);
        sub_in_place::<B>(&mut q, &[1]);
    }
    let mut r = sub::<B>(a, &prod);
    while cmp(&r, b) != Ordering::Less {
        sub_in_place::<B>(&mut r, b);
        add_small::<B>(&mut q, 1);
    }
    (q, r)
}

/// `(floor(a / b), a mod b)` for `b > 0`.
pub fn divrem<const B: u32>(a: &[u32], b: &[u32]) -> (Vec<u32>, Vec<u32>) {
    assert!(!b.is_empty(), "division by zero");
    if cmp(a, b) == Ordering::Less {
        return (Vec::new(), a.to_vec());
    }
    if b.len() == 1 {
        let (q, r) = divrem_small::<B>(a, b[0]);
        return (q, if r == 0 { Vec::new() } else { vec![r] });
    }
    let (n, m) = (a.len(), b.len());
    let k = n - m + 1;
    if m <= SCHOOL || k <= SCHOOL {
        return knuth::<B>(a, b);
    }
    let d = B / (b[m - 1] + 1);
    let (a, b) = (mul_small::<B>(a, d), mul_small::<B>(b, d));
    let (n, m) = (a.len(), b.len());
    let (q, r) = if n <= 2 * m {
        let t = k + 2;
        if m > t {
            let s = m - t;
            let (a2, b2) = (&a[s..], &b[s..]);
            let est = estimate::<B>(a2, b2, &reciprocal::<B>(b2));
            correct::<B>(&a, &b, est)
        } else {
            let est = estimate::<B>(&a, &b, &reciprocal::<B>(&b));
            correct::<B>(&a, &b, est)
        }
    } else {
        // long quotient: chunks of m limbs from the top, one shared reciprocal
        let rec = reciprocal::<B>(&b);
        let mut q = vec![0u32; n - m + 1];
        let mut rem: Vec<u32> = Vec::new();
        let mut pos = n;
        while pos > 0 {
            let c = pos.min(m);
            let mut x = a[pos - c..pos].to_vec();
            if !rem.is_empty() {
                let mut hi = shl(&rem, c);
                hi[..c].copy_from_slice(&x);
                x = hi;
            }
            trim(&mut x);
            let (qc, r) = if cmp(&x, &b) == Ordering::Less { (Vec::new(), x) } else { correct::<B>(&x, &b, estimate::<B>(&x, &b, &rec)) };
            for (i, &v) in qc.iter().enumerate() {
                q[pos - c + i] = v;
            }
            rem = r;
            pos -= c;
        }
        trim(&mut q);
        (q, rem)
    };
    let (r, rem) = divrem_small::<B>(&r, d);
    debug_assert_eq!(rem, 0);
    (q, r)
}
