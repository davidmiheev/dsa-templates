//! Stress test: primality and factorisation against trial division, primitive roots by order.

use number_theory::{factorize, gcd, is_prime, mod_inverse, mod_sqrt, mul_mod, pow_mod, primitive_root, NONE};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn trial_factor(mut n: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut p = 2;
    while p * p <= n {
        if n.is_multiple_of(p) {
            let mut e = 0;
            while n.is_multiple_of(p) {
                n /= p;
                e += 1;
            }
            out.push((p, e));
        }
        p += 1;
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

#[test]
fn small_numbers_match_trial_division() {
    let mut st = Stress::new("number_theory/small");
    for n in 1..20_000u64 {
        let f = trial_factor(n);
        assert_eq!(factorize(n), f, "n={n}");
        assert_eq!(is_prime(n), f.len() == 1 && f[0].1 == 1, "n={n}");
        st.ops(1);
    }
    st.case();
    st.done();
}

#[test]
fn large_numbers() {
    let mut st = Stress::new("number_theory/large");
    let mut rng = Rng::new(1);
    for case in 0..300 * scale() {
        // products of random primes near 2^30..2^31 and arbitrary 64-bit numbers
        let n = match case % 3 {
            0 => rng.next_u64() | 1,
            1 => {
                let mut p = || loop {
                    let c = (rng.next_u64() >> 33) | (1 << 30);
                    if is_prime(c) {
                        return c;
                    }
                };
                p() * p()
            }
            _ => 1 + rng.next_u64() % 1_000_000_000_000,
        };
        let f = factorize(n);
        let prod = f.iter().fold(1u128, |acc, &(p, e)| acc * (p as u128).pow(e));
        assert_eq!(prod, n as u128);
        assert!(f.iter().all(|&(p, _)| is_prime(p)));
        assert!(f.windows(2).all(|w| w[0].0 < w[1].0));
        st.case();
    }
    // known values
    assert!(is_prime(18446744073709551557)); // largest prime below 2^64
    assert!(!is_prime(3215031751)); // strong pseudoprime to bases 2, 3, 5, 7
    assert_eq!(factorize(u64::MAX), vec![(3, 1), (5, 1), (17, 1), (257, 1), (641, 1), (65537, 1), (6700417, 1)]);
    st.done();
}

#[test]
fn modular_helpers_and_primitive_roots() {
    let mut rng = Rng::new(2);
    for _ in 0..1000 {
        let m = 2 + rng.next_u64() % 1_000_000_007;
        let a = rng.next_u64() % m;
        if gcd(a, m) == 1 {
            assert_eq!(mul_mod(a, mod_inverse(a, m), m), 1 % m);
        }
        assert_eq!(pow_mod(a, 0, m), 1 % m);
    }
    for p in (2..400u64).filter(|&p| is_prime(p)) {
        let g = primitive_root(p);
        let mut x = 1;
        let mut seen = std::collections::HashSet::new();
        for _ in 0..p - 1 {
            x = x * g % p;
            seen.insert(x);
        }
        assert_eq!(seen.len() as u64, p - 1, "p={p} g={g}");
    }
    assert_eq!(primitive_root(998_244_353), 3);
    assert_eq!(primitive_root(1_000_000_007), 5);
}

#[test]
fn square_roots() {
    for p in (2..300u64).filter(|&p| is_prime(p)).chain([998_244_353, 1_000_000_007, 469_762_049]) {
        let mut rng = Rng::new(p);
        let probes: Vec<u64> = if p < 300 { (0..p).collect() } else { (0..200).map(|_| rng.next_u64() % p).collect() };
        for a in probes {
            let r = mod_sqrt(a, p);
            let residue = a == 0 || p == 2 || pow_mod(a, (p - 1) / 2, p) == 1;
            if residue {
                assert_ne!(r, NONE, "a={a} p={p}");
                assert_eq!(mul_mod(r, r, p), a % p, "a={a} p={p}");
            } else {
                assert_eq!(r, NONE, "a={a} p={p}");
            }
        }
    }
}

use number_theory::*;

#[test]
fn roots_floor_sums_and_totient() {
    let mut st = Stress::new("number_theory/arith");
    for case in 0..20_000 * scale() {
        let mut rng = Rng::new(case as u64);
        // integer roots, including values right at the powers
        let k = rng.range(1, 65) as u32;
        let a = match rng.below(3) {
            0 => rng.next_u64(),
            1 => rng.next_u64() >> rng.below(64),
            _ => {
                let r = (rng.next_u64() >> 1) % 1_000_000;
                let x = (r as u128).pow(k.min(3)).min(u64::MAX as u128) as u64;
                x.saturating_add(rng.below(3) as u64).saturating_sub(1)
            }
        };
        let r = kth_root(a, k);
        let pw = |x: u64| (0..k).try_fold(1u128, |acc, _| acc.checked_mul(x as u128).filter(|&v| v <= u64::MAX as u128));
        assert!(pw(r).is_some_and(|v| v <= a as u128), "a={a} k={k} r={r}");
        assert!(r == u64::MAX || pw(r + 1).is_none_or(|v| v > a as u128), "a={a} k={k} r={r}");
        // floor sum and minimum of a linear function
        let n = rng.range(1, 40) as u64;
        let m = rng.range(1, 40) as u64;
        let (aa, bb) = (rng.below(100) as u64, rng.below(100) as u64);
        let sum: u128 = (0..n).map(|i| ((aa * i + bb) / m) as u128).sum();
        assert_eq!(floor_sum(n, m, aa, bb), sum);
        let mn = (0..n).map(|x| (aa * x + bb) % m).min().unwrap();
        assert_eq!(min_mod_linear(n, m, aa, bb), mn, "n={n} m={m} a={aa} b={bb}");
        // euler phi
        let v = rng.range(1, 2000) as u64;
        assert_eq!(euler_phi(v), (1..=v).filter(|&i| gcd(i, v) == 1).count() as u64);
        st.case();
        st.ops(1);
    }
    // huge arguments
    let big = (1u64 << 62) - 3;
    assert_eq!(isqrt(u64::MAX), u32::MAX as u64);
    assert_eq!(kth_root(u64::MAX, 2), u32::MAX as u64);
    assert_eq!(kth_root(u64::MAX, 64), 1);
    assert_eq!(kth_root(1 << 63, 63), 2);
    assert_eq!(min_mod_linear(big, big, big - 1, 5), 0);
    assert_eq!(floor_sum(1_000_000_000, 1_000_000_000, 999_999_999, 999_999_999), {
        let (n, m, a, b) = (1_000_000_000u128, 1_000_000_000u128, 999_999_999u128, 999_999_999u128);
        (0..n).map(|i| (a * i + b) / m).sum::<u128>()
    });
    st.done();
}

#[test]
fn discrete_logs_and_roots() {
    let mut st = Stress::new("number_theory/dlog");
    // general modulus against brute force
    for m in 1..=60u64 {
        for x in 0..m {
            for y in 0..m {
                let want = (0..=2 * m).find(|&k| (if k == 0 { 1 % m } else { pow_mod(x, k, m) }) == y);
                assert_eq!(discrete_log(x, y, m), want.unwrap_or(NONE), "x={x} y={y} m={m}");
            }
        }
        st.case();
    }
    let primes: Vec<u64> = (2..400u64).filter(|&n| is_prime(n)).chain([998_244_353, 1_000_000_007, 999_999_937, 65_537, 4_294_967_291]).collect();
    for &p in &primes {
        let g = primitive_root(p);
        let table = PrimeLog::new(p, g, 50);
        let mut rng = Rng::new(p);
        for _ in 0..40 {
            let a = 1 + rng.below(p as usize - 1) as u64;
            let b = discrete_log_prime(a, g, p);
            assert!(b < p - 1 && pow_mod(g, b, p) == a, "p={p} a={a}");
            assert_eq!(table.log(a), b);
            assert_eq!(discrete_log(g, a, p), b);
            // k-th roots
            let k = rng.below(2 * p as usize) as u64;
            let r = kth_root_mod(k, a, p);
            if r == NONE {
                if p < 400 {
                    assert!((0..p).all(|x| pow_mod(x, k, p) != a), "p={p} k={k} a={a}");
                } else {
                    let d = gcd(k, p - 1);
                    assert!(pow_mod(a, (p - 1) / d, p) != 1, "p={p} k={k} a={a}");
                }
            } else {
                assert_eq!(pow_mod(r, k, p), a, "p={p} k={k} a={a}");
            }
            st.ops(1);
        }
        st.case();
    }
    for p in [2u64, 3, 5] {
        for k in 0..6 {
            for y in 0..p {
                let r = kth_root_mod(k, y, p);
                let any = (0..p).any(|x| (if k == 0 { 1 } else { pow_mod(x, k, p) }) == y);
                assert_eq!(r != NONE, any, "p={p} k={k} y={y}");
                if r != NONE {
                    assert_eq!(if k == 0 { 1 } else { pow_mod(r, k, p) }, y);
                }
            }
        }
    }
    st.done();
}

fn tower_exact(a: u64, b: u64) -> Option<u128> {
    let mut v: u128 = 1;
    for _ in 0..b {
        v = (a as u128).checked_pow(u32::try_from(v).ok()?)?;
    }
    Some(v)
}

#[test]
fn tetration_matches_exact_towers() {
    for a in 0..8u64 {
        for b in 0..6u64 {
            if let Some(v) = tower_exact(a, b) {
                for m in 1..300u64 {
                    assert_eq!(tetration_mod(a, b, m), (v % m as u128) as u64, "{a}^^{b} mod {m}");
                }
            }
        }
    }
    // a known value: 3^^3 = 3^27 = 7625597484987
    assert_eq!(tetration_mod(3, 3, 1_000_000_007), 7_625_597_484_987 % 1_000_000_007);
    // towers of a fixed height agree once b is large (they stabilise modulo m)
    for m in [7u64, 100, 1_000_000_000, 999_999_937] {
        assert_eq!(tetration_mod(2, 100, m), tetration_mod(2, 101, m));
    }
}

#[test]
fn nim_product_field_axioms() {
    let nm = NimMultiplier::new();
    assert_eq!(nm.mul(2, 2), 3);
    assert_eq!(nm.mul(4, 4), 6);
    assert_eq!(nm.mul(8, 8), 13);
    assert_eq!(nm.mul(16, 16), 24);
    // the defining mex recursion for small values
    let n = 40;
    let mut t = vec![vec![0u64; n]; n];
    for a in 0..n {
        for b in 0..n {
            let mut seen = vec![false; 512];
            for a2 in 0..a {
                for b2 in 0..b {
                    let v = t[a2][b] ^ t[a][b2] ^ t[a2][b2];
                    if (v as usize) < seen.len() {
                        seen[v as usize] = true;
                    }
                }
            }
            t[a][b] = seen.iter().position(|&s| !s).unwrap() as u64;
            assert_eq!(nm.mul(a as u64, b as u64), t[a][b], "{a} x {b}");
        }
    }
    let mut rng = Rng::new(5);
    for _ in 0..20_000 {
        let (a, b, c) = (rng.next_u64(), rng.next_u64(), rng.next_u64());
        assert_eq!(nm.mul(a, b), nm.mul(b, a));
        assert_eq!(nm.mul(nm.mul(a, b), c), nm.mul(a, nm.mul(b, c)));
        assert_eq!(nm.mul(a, b ^ c), nm.mul(a, b) ^ nm.mul(a, c));
        assert_eq!(nm.mul(a, 1), a);
        assert_eq!(nm.mul(a, 0), 0);
    }
}

#[test]
fn gaussian_gcd_and_two_squares() {
    // g divides x: x / g is a Gaussian integer
    let divides = |g: (i64, i64), x: (i64, i64)| {
        let n = g.0 as i128 * g.0 as i128 + g.1 as i128 * g.1 as i128;
        let re = x.0 as i128 * g.0 as i128 + x.1 as i128 * g.1 as i128;
        let im = x.1 as i128 * g.0 as i128 - x.0 as i128 * g.1 as i128;
        re % n == 0 && im % n == 0
    };
    let mut rng = Rng::new(3);
    for _ in 0..20_000 {
        let lim = [3i64, 20, 1_000_000_000][rng.below(3)];
        let mut p = || rng.range(-lim, lim + 1);
        let (a, b, c, d) = ((p(), p()), (p(), p()), (p(), p()), (p(), p()));
        // common factor c*d... build two numbers with a known common factor
        let mul = |x: (i64, i64), y: (i64, i64)| -> (i64, i64) {
            (((x.0 as i128 * y.0 as i128 - x.1 as i128 * y.1 as i128) % (1 << 60)) as i64, ((x.0 as i128 * y.1 as i128 + x.1 as i128 * y.0 as i128) % (1 << 60)) as i64)
        };
        let (_, _) = (c, d);
        let g = gaussian_gcd(a, b);
        if a == (0, 0) && b == (0, 0) {
            assert_eq!(g, (0, 0));
            continue;
        }
        assert!(divides(g, a) && divides(g, b), "{a:?} {b:?} -> {g:?}");
        // norm of the gcd equals that of the gcd of the quotients times itself: g * gcd(a/g, b/g) is a unit
        let n = g.0 as i128 * g.0 as i128 + g.1 as i128 * g.1 as i128;
        let q = |x: (i64, i64)| -> (i64, i64) {
            (((x.0 as i128 * g.0 as i128 + x.1 as i128 * g.1 as i128) / n) as i64, ((x.1 as i128 * g.0 as i128 - x.0 as i128 * g.1 as i128) / n) as i64)
        };
        let h = gaussian_gcd(q(a), q(b));
        assert_eq!(h.0.abs() + h.1.abs(), 1, "gcd of the quotients must be a unit");
        // a multiple of g shares at least g
        let m1 = mul(a, (3, 4));
        if m1.0.abs() < 1 << 50 {
            let g2 = gaussian_gcd(m1, mul(b, (3, 4)));
            assert!(divides(g, g2) || g == (0, 0));
        }
    }
    // sums of two squares against brute force
    for n in 0..3000u64 {
        let want: Vec<(u64, u64)> = (0..=isqrt(n)).flat_map(|a| (0..=isqrt(n)).map(move |b| (a, b))).filter(|&(a, b)| a * a + b * b == n).collect();
        assert_eq!(two_squares(n), want, "n={n}");
    }
    let mut rng = Rng::new(9);
    for _ in 0..300 {
        let n = rng.next_u64() % 1_000_000_000_000_000_000;
        let sq = two_squares(n);
        let mut seen = std::collections::HashSet::new();
        for &(a, b) in &sq {
            assert_eq!(a as u128 * a as u128 + b as u128 * b as u128, n as u128);
            assert!(seen.insert((a, b)));
            assert!(seen.contains(&(a, b)));
        }
        // symmetric under swapping
        assert!(sq.iter().all(|&(a, b)| sq.binary_search(&(b, a)).is_ok()));
    }
}

#[test]
fn primality_at_the_top_of_u64() {
    assert!(is_prime(18_446_744_073_709_551_557)); // 2^64 - 59, the largest prime below 2^64
    assert!(!is_prime(u64::MAX));
    assert!(!is_prime(3_825_123_056_546_413_051)); // strong pseudoprime to the first nine primes
    assert!(!is_prime(4_294_967_291u64 * 4_294_967_291)); // (2^32 - 5)^2
    assert!(!is_prime(4_294_967_291u64 * 4_294_967_279));
    assert!(is_prime(4_294_967_291) && is_prime(4_294_967_279));
    for c in [561u64, 1105, 1729, 2465, 2821, 6601, 8911, 41041, 825_265, 321_197_185, 5_394_826_801] {
        assert!(!is_prime(c), "Carmichael number {c}");
    }
    let f = factorize(4_294_967_291u64 * 4_294_967_279);
    assert_eq!(f, vec![(4_294_967_279, 1), (4_294_967_291, 1)]);
    assert_eq!(factorize(4_294_967_291u64 * 4_294_967_291), vec![(4_294_967_291, 2)]);
}

#[test]
fn prime_log_batches_match_single_queries() {
    for p in [2u64, 3, 5, 7, 101, 65_537, 998_244_353, 1_000_000_007] {
        let g = primitive_root(p);
        let table = PrimeLog::new(p, g, 1000);
        let mut rng = Rng::new(p);
        let qs: Vec<u64> = (0..500).map(|_| 1 + rng.below(p as usize - 1) as u64).collect();
        let want: Vec<u64> = qs.iter().map(|&a| table.log(a)).collect();
        assert_eq!(table.log_many(&qs), want, "p={p}");
        assert!(want.iter().zip(&qs).all(|(&b, &a)| pow_mod(g, b, p) == a));
        assert!(table.log_many(&[]).is_empty());
    }
}
