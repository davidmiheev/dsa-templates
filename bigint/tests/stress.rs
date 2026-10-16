//! Stress tests of the big-integer arithmetic against i128 and identities.

use bigint::{BigUint, Decimal, Hex};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

type D = BigUint<1_000_000>;
type H = BigUint<1_048_576>;

fn random_digits(rng: &mut Rng, len: usize, radix: u32, shape: usize) -> Vec<u8> {
    let digit = |v: u32| char::from_digit(v, radix).expect("digit") as u8;
    let mut s: Vec<u8> = (0..len)
        .map(|_| match shape {
            0 => digit(rng.below(radix as usize) as u32),
            1 => digit(radix - 1), // all nines: extreme carries
            2 => digit(if rng.chance(1, 8) { rng.below(radix as usize) as u32 } else { 0 }),
            _ => digit(if rng.chance(1, 2) { radix - 1 } else { 0 }),
        })
        .collect();
    if let Some(first) = s.first_mut() {
        if *first == b'0' {
            *first = digit(1 + rng.below(radix as usize - 1) as u32);
        }
    }
    s
}

fn naive_mul(a: &[u8], b: &[u8], radix: u32) -> Vec<u8> {
    let av: Vec<u32> = a.iter().rev().map(|&c| (c as char).to_digit(radix).expect("digit")).collect();
    let bv: Vec<u32> = b.iter().rev().map(|&c| (c as char).to_digit(radix).expect("digit")).collect();
    let mut acc = vec![0u64; av.len() + bv.len() + 1];
    for (i, &x) in av.iter().enumerate() {
        for (j, &y) in bv.iter().enumerate() {
            acc[i + j] += x as u64 * y as u64;
        }
    }
    let mut carry = 0u64;
    for x in acc.iter_mut() {
        let t = *x + carry;
        *x = t % radix as u64;
        carry = t / radix as u64;
    }
    while acc.len() > 1 && *acc.last().unwrap() == 0 {
        acc.pop();
    }
    acc.iter().rev().map(|&d| char::from_digit(d as u32, radix).unwrap() as u8).collect()
}

#[test]
fn small_values_against_i128() {
    let mut st = Stress::new("bigint/i128");
    for case in 0..20_000 * scale() {
        let mut rng = Rng::new(case as u64);
        let mask = |rng: &mut Rng| -> i128 {
            let v = (rng.next_u64() >> rng.below(64)) as i128;
            if rng.chance(1, 2) {
                -v
            } else {
                v
            }
        };
        let (a, b) = (mask(&mut rng), mask(&mut rng));
        let (da, db) = (Decimal::parse(a.to_string().as_bytes()), Decimal::parse(b.to_string().as_bytes()));
        assert_eq!(da.add(&db).to_bytes(), (a + b).to_string().into_bytes(), "{a} + {b}");
        assert_eq!(da.mul(&db).to_bytes(), (a * b).to_string().into_bytes(), "{a} * {b}");
        let hex = |x: i128| if x < 0 { format!("-{:x}", -x) } else { format!("{:x}", x) };
        let (ha, hb) = (Hex::parse(hex(a).as_bytes()), Hex::parse(hex(b).as_bytes()));
        assert_eq!(ha.add(&hb).to_bytes(), hex(a + b).into_bytes());
        assert_eq!(ha.mul(&hb).to_bytes(), hex(a * b).into_bytes());
        let (ua, ub) = (a.unsigned_abs(), b.unsigned_abs());
        if ub > 0 {
            let (q, r) = D::from_digits(ua.to_string().as_bytes(), 10).divrem(&D::from_digits(ub.to_string().as_bytes(), 10));
            assert_eq!((q.to_digits(10), r.to_digits(10)), ((ua / ub).to_string().into_bytes(), (ua % ub).to_string().into_bytes()));
        }
        st.case();
    }
    st.done();
}

#[test]
fn multiplication_matches_schoolbook() {
    let mut st = Stress::new("bigint/mul");
    for case in 0..300 * scale() {
        let mut rng = Rng::new(7_000 + case as u64);
        let (la, lb) = (1 + rng.below(700), 1 + rng.below(700));
        for radix in [10u32, 16] {
            let (sa, sb) = (rng.below(4), rng.below(4));
            let a = random_digits(&mut rng, la, radix, sa);
            let b = random_digits(&mut rng, lb, radix, sb);
            let want = naive_mul(&a, &b, radix);
            let got = if radix == 10 { D::from_digits(&a, 10).mul(&D::from_digits(&b, 10)).to_digits(10) } else { H::from_digits(&a, 16).mul(&H::from_digits(&b, 16)).to_digits(16) };
            assert_eq!(got, want, "radix {radix} lengths {la} {lb}");
        }
        st.case();
    }
    st.done();
}

fn check_division<const BASE: u32>(rng: &mut Rng, radix: u32, lb: usize, lq: usize) {
    type U<const B: u32> = BigUint<B>;
    let shape = rng.below(4);
    let b = U::<BASE>::from_digits(&random_digits(rng, lb, radix, shape), radix);
    let q_zero = rng.chance(1, 8);
    let shape = rng.below(4);
    let q = if q_zero { U::<BASE>::zero() } else { U::<BASE>::from_digits(&random_digits(rng, lq, radix, shape), radix) };
    // remainder: zero, one less than b, or random below b
    let r = match rng.below(3) {
        0 => U::<BASE>::zero(),
        1 => b.sub(&U::<BASE>::from_u64(1)),
        _ => {
            let len = 1 + rng.below(lb);
            let digits = random_digits(rng, len, radix, 0);
            U::<BASE>::from_digits(&digits, radix).divrem(&b).1
        }
    };
    let a = q.mul(&b).add(&r);
    let (gq, gr) = a.divrem(&b);
    assert!(gq == q && gr == r, "radix {radix}: |b| = {lb}, |q| = {lq}");
}

#[test]
fn division_identities() {
    let mut st = Stress::new("bigint/div");
    // lengths straddle the Knuth / Newton thresholds (40 limbs = 160 digits)
    let lens = [1usize, 3, 4, 5, 9, 40, 100, 159, 160, 161, 200, 330, 700, 1500, 4000];
    for case in 0..1500 * scale() {
        let mut rng = Rng::new(9_000 + case as u64);
        let lb = lens[rng.below(lens.len())] + rng.below(5);
        let lq = lens[rng.below(lens.len())] + rng.below(5);
        check_division::<1_000_000>(&mut rng, 10, lb, lq);
        check_division::<1_048_576>(&mut rng, 16, lb, lq);
        check_division::<100>(&mut rng, 10, lb / 2 + 1, lq / 2 + 1);
        check_division::<256>(&mut rng, 16, lb / 2 + 1, lq / 2 + 1);
        st.case();
    }
    st.done();
}

#[test]
fn parsing_roundtrip() {
    let mut rng = Rng::new(1);
    for len in [1usize, 2, 3, 4, 5, 8, 13, 100] {
        let s = random_digits(&mut rng, len, 10, 0);
        assert_eq!(Decimal::parse(&s).to_bytes(), s);
        let h = random_digits(&mut rng, len, 16, 0);
        assert_eq!(Hex::parse(&h).to_bytes(), h);
    }
    assert_eq!(Decimal::parse(b"0").to_bytes(), b"0");
    assert_eq!(Decimal::parse(b"-0").to_bytes(), b"0");
    assert_eq!(Decimal::parse(b"-5").add(&Decimal::parse(b"5")).to_bytes(), b"0");
}

#[test]
fn huge_operands_use_three_primes() {
    // (B^n - 1)^2 = B^{2n} - 2 B^n + 1 with B = 16^5: the coefficients exceed the two-prime range
    let n = 450_000usize;
    let a = vec![b'f'; 5 * n];
    let got = H::from_digits(&a, 16).mul(&H::from_digits(&a, 16)).to_digits(16);
    let mut want = vec![b'f'; 5 * (n - 1)];
    want.extend_from_slice(b"ffffe");
    want.extend(std::iter::repeat_n(b'0', 5 * (n - 1)));
    want.extend_from_slice(b"00001");
    assert_eq!(got, want);
}
