//! Stress test: `RollingHash` / `DoubleRollingHash` against the polynomial definition
//! `h(s) = sum s[i] * base^(n-1-i) mod m` evaluated with an independent power routine.

use rolling_hash::{compute_hash, DoubleRollingHash, RollingHash, DEFAULT_BASE, DEFAULT_MOD};
use std::collections::HashMap;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

/// (base, modulus) pairs; moduli stay above 256 so a single byte is already reduced.
const PARAMS: [(u64, u64); 5] = [
    (131, 2_305_843_009_213_693_951), // 2^61 - 1
    (911_382_323, 998_244_353),
    (31, 1_000_000_007),
    (257, 65_537),
    (5, 1_000_003),
];

fn pow(b: u64, mut e: u64, m: u64) -> u64 {
    let (mut acc, mut base) = (1u128, b as u128 % m as u128);
    while e > 0 {
        if e & 1 == 1 {
            acc = acc * base % m as u128;
        }
        base = base * base % m as u128;
        e >>= 1;
    }
    acc as u64
}

/// Definition, term by term (not Horner).
fn definition(s: &[u8], base: u64, m: u64) -> u64 {
    let n = s.len();
    s.iter()
        .enumerate()
        .fold(0u128, |acc, (i, &c)| (acc + c as u128 * pow(base, (n - 1 - i) as u64, m) as u128) % m as u128) as u64
}

fn bytes(rng: &mut Rng, len: usize) -> Vec<u8> {
    // Byte 0 is avoided: it contributes nothing to a polynomial hash (see `leading_zero_bytes_are_invisible`).
    let alphabet = rng.pick(&[2usize, 4, 26, 255]);
    (0..len).map(|_| (rng.below(alphabet) + if alphabet == 26 { 97 } else { 1 }) as u8).collect()
}

#[test]
fn compute_hash_and_push_match_definition() {
    let mut st = Stress::new("rolling_hash/push");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(case as u64);
        let (base, m) = PARAMS[case % PARAMS.len()];
        let len = rng.below(40);
        let s = bytes(&mut rng, len);
        let mut h = RollingHash::with_base_mod(base, m);
        st.case();
        assert_eq!(h.value(), 0, "empty hash is 0");
        for (i, &c) in s.iter().enumerate() {
            h.push(c);
            assert_eq!(h.value(), definition(&s[..=i], base, m), "prefix {} case {case}", i + 1);
            assert_eq!(h.value(), compute_hash(&s[..=i], base, m), "compute_hash prefix {}", i + 1);
            st.ops(1);
        }
        h.clear();
        assert_eq!(h.value(), 0, "clear resets");
    }
    st.done();
}

#[test]
fn sliding_window_matches_recomputation() {
    let mut st = Stress::new("rolling_hash/slide");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(3000 + case as u64);
        let (base, m) = PARAMS[case % PARAMS.len()];
        let len = rng.range(1, 60) as usize;
        let w = rng.range(1, len as i64) as usize;
        let s = bytes(&mut rng, len);
        let mut h = RollingHash::with_base_mod(base, m);
        for &c in &s[..w] {
            h.push(c);
        }
        st.case();
        assert_eq!(h.value(), definition(&s[..w], base, m));
        for i in 0..len - w {
            h.pop(s[i], w); // window shrinks to w-1 bytes...
            assert_eq!(h.value(), definition(&s[i + 1..i + w], base, m), "after pop, case {case}");
            h.push(s[i + w]); // ...then grows back to w
            assert_eq!(h.value(), definition(&s[i + 1..i + 1 + w], base, m), "window at {} case {case}", i + 1);
            st.ops(2);
        }
    }
    st.done();
}

#[test]
fn double_hash_pairs_two_single_hashes() {
    let mut st = Stress::new("rolling_hash/double");
    for case in 0..300 * scale() {
        let mut rng = Rng::new(9000 + case as u64);
        let ((b1, m1), (b2, m2)) = (PARAMS[case % 5], PARAMS[(case + 1) % 5]);
        let len = rng.below(50);
        let s = bytes(&mut rng, len);
        let mut d = DoubleRollingHash::new(b1, m1, b2, m2);
        assert_eq!((d.pub_b1, d.pub_m1, d.pub_b2, d.pub_m2), (b1, m1, b2, m2));
        for &c in &s {
            d.push(c);
        }
        assert_eq!(d.value(), (compute_hash(&s, b1, m1), compute_hash(&s, b2, m2)), "case {case}");
        d.clear();
        assert_eq!(d.value(), (0, 0));
        st.case();
        st.ops(s.len() + 2);
    }
    st.done();
}

#[test]
fn default_hash_has_no_collisions_among_random_strings() {
    // With a 61-bit modulus, 50,000 random strings colliding would mean a broken hash.
    let mut st = Stress::new("rolling_hash/collisions");
    let mut rng = Rng::new(17);
    let mut seen: HashMap<u64, Vec<u8>> = HashMap::new();
    st.case();
    for _ in 0..50_000 * scale() {
        let len = rng.range(1, 24) as usize;
        // ASCII only: DEFAULT_BASE (131) is below 256, so larger byte values would carry
        // into the next digit (see `bytes_above_the_base_collide`).
        let s: Vec<u8> = bytes(&mut rng, len).into_iter().map(|b| b % 100 + 1).collect();
        let mut h = RollingHash::new();
        for &c in &s {
            h.push(c);
        }
        assert_eq!(h.value(), compute_hash(&s, DEFAULT_BASE, DEFAULT_MOD));
        if let Some(prev) = seen.insert(h.value(), s.clone()) {
            assert_eq!(prev, s, "hash collision between distinct strings");
        }
        st.ops(1);
    }
    st.done();
}

/// Documents a property of any plain polynomial hash over bytes: a leading 0x00 byte adds
/// `0 * base^k`, so `[0, x]` and `[x]` hash equally. Text input never contains NUL, but
/// map bytes to `b + 1` first if binary data must be distinguished.
#[test]
fn leading_zero_bytes_are_invisible() {
    for (base, m) in PARAMS {
        assert_eq!(compute_hash(&[0, 1], base, m), compute_hash(&[1], base, m));
        assert_eq!(compute_hash(&[0, 0, 7, 9], base, m), compute_hash(&[7, 9], base, m));
        assert_ne!(compute_hash(&[1, 0], base, m), compute_hash(&[1], base, m), "trailing zeros do count");
    }
}

/// Documents a limitation of `DEFAULT_BASE = 131`: a byte above the base overflows into the
/// neighbouring digit, so short byte strings collide by construction, not by chance.
/// ASCII text (bytes < 128) is unaffected; choose a base above 256, ideally random, for binary data.
#[test]
fn bytes_above_the_base_collide() {
    let h = |s: &[u8]| compute_hash(s, DEFAULT_BASE, DEFAULT_MOD);
    assert_eq!(h(&[1, 1]), h(&[(DEFAULT_BASE + 1) as u8]), "1*131 + 1 == 132");
    assert_ne!(h(&[1, 1]), h(&[100]), "ASCII-range digits stay distinct");
}
