//! Stress test: bit helpers against bit-by-bit loops (exhaustive for 16-bit values,
//! random and edge-case values for 64-bit), and `Bitset` against `Vec<bool>`.

use bit_manipulation::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn bit(x: u64, i: u32) -> bool {
    (x >> i) & 1 == 1
}
fn n_popcount(x: u64) -> u32 {
    (0..64).filter(|&i| bit(x, i)).count() as u32
}
fn n_trailing_zeros(x: u64) -> u32 {
    (0..64).find(|&i| bit(x, i)).unwrap_or(64)
}
fn n_trailing_ones(x: u64) -> u32 {
    (0..64).find(|&i| !bit(x, i)).unwrap_or(64)
}
fn n_highest_bit(x: u64) -> u64 {
    (0..64).rev().find(|&i| bit(x, i)).map_or(0, |i| 1u64 << i)
}
fn n_lowest_bit(x: u64) -> u64 {
    (0..64).find(|&i| bit(x, i)).map_or(0, |i| 1u64 << i)
}
fn n_next_pow2(x: u64) -> u64 {
    let mut p = 1u64;
    while p < x {
        p <<= 1;
    }
    p
}

fn check_value(x: u64) {
    assert_eq!(popcount(x), n_popcount(x), "popcount({x})");
    assert_eq!(trailing_zeros(x), n_trailing_zeros(x), "trailing_zeros({x})");
    assert_eq!(trailing_ones(x), n_trailing_ones(x), "trailing_ones({x})");
    assert_eq!(highest_bit(x), n_highest_bit(x), "highest_bit({x})");
    assert_eq!(lowest_bit(x), n_lowest_bit(x), "lowest_bit({x})");
    assert_eq!(is_power_of_two(x), n_popcount(x) == 1, "is_power_of_two({x})");
    if x <= 1 << 63 {
        assert_eq!(next_power_of_two(x), n_next_pow2(x), "next_power_of_two({x})");
    }
    if x > 0 {
        assert_eq!(log2_floor(x), 63 - n_highest_bit(x).leading_zeros(), "log2_floor({x})");
        assert_eq!(1u64 << log2_floor(x), n_highest_bit(x), "log2_floor is the position of the top bit ({x})");
    }
}

#[test]
fn free_functions_match_bit_loops() {
    let mut st = Stress::new("bit_manipulation/fns");
    st.case();
    for x in 0..=u16::MAX as u64 {
        check_value(x);
        st.ops(9);
    }
    let mut rng = Rng::new(1);
    let mut edges = vec![u64::MAX, u64::MAX - 1, 1 << 63, (1 << 63) + 1, (1 << 63) - 1, 1 << 32, (1 << 32) - 1];
    for k in 0..64 {
        edges.extend([1u64 << k, (1u64 << k).wrapping_sub(1), (1u64 << k).wrapping_add(1)]);
    }
    for x in edges {
        check_value(x);
        st.ops(9);
    }
    for _ in 0..100_000 * scale() {
        let x = rng.next_u64() >> rng.below(64); // varied magnitudes
        check_value(x);
        st.ops(9);
    }
    st.done();
}

#[test]
#[should_panic(expected = "log2_floor requires x > 0")]
fn log2_floor_of_zero_panics() {
    log2_floor(0);
}

#[test]
fn subsets_enumerates_all_nonempty_masks() {
    let mut st = Stress::new("bit_manipulation/subsets");
    for k in 0..=16u32 {
        let got: Vec<u64> = subsets(k).collect();
        st.case();
        // 2^k - 1 non-empty subsets (the crate docs say 2^k), each exactly once, in increasing order.
        assert_eq!(got.len() as u64, (1u64 << k) - 1, "count for k={k}");
        assert!(got.windows(2).all(|w| w[0] < w[1]), "strictly increasing for k={k}");
        assert!(got.iter().all(|&m| m != 0 && m < 1 << k), "range for k={k}");
        st.ops(got.len());
    }
    st.done();
}

#[test]
fn bitset_matches_vec_of_bool() {
    let mut st = Stress::new("bit_manipulation/bitset");
    for case in 0..300 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 200) as usize; // includes sizes that are not multiples of 64
        let mut bs = Bitset::new(n);
        let mut model = vec![false; n];
        st.case();
        for step in 0..200 {
            let i = rng.below(n);
            match rng.below(5) {
                0 => {
                    bs.set(i);
                    model[i] = true;
                }
                1 => {
                    bs.reset(i);
                    model[i] = false;
                }
                2 => {
                    bs.flip(i);
                    model[i] = !model[i];
                }
                3 => assert_eq!(bs.get(i), model[i], "get({i}) case {case} step {step}"),
                _ => assert_eq!(
                    bs.count_ones(),
                    model.iter().filter(|&&b| b).count() as u32,
                    "count_ones case {case} step {step}"
                ),
            }
            st.ops(1);
        }
        for i in 0..n {
            assert_eq!(bs.get(i), model[i], "final get({i}) case {case}");
        }
        st.ops(n);
    }
    st.done();
}

#[test]
fn bitset_out_of_range_panics() {
    let mut bs = Bitset::new(70);
    assert!(std::panic::catch_unwind(|| Bitset::new(70).get(70)).is_err());
    assert!(std::panic::catch_unwind(move || bs.set(70)).is_err());
    // bits 64..70 live in the second word; index 69 is valid
    let mut ok = Bitset::new(70);
    ok.set(69);
    assert!(ok.get(69) && ok.count_ones() == 1);
}
