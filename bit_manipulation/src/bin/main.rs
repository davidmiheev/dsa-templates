use bit_manipulation::{
    Bitset, highest_bit, is_power_of_two, log2_floor, lowest_bit, next_power_of_two, popcount,
    subsets, trailing_ones, trailing_zeros,
};

/// Run the bit-manipulation example.
pub fn run() {
    println!("bit_manipulation");

    let x: u64 = 0b1011_0110_u64;

    println!("{}", popcount(x));
    println!("{}", trailing_zeros(x));
    println!("{}", trailing_ones(x));
    println!("{}", highest_bit(x));
    println!("{}", lowest_bit(x));
    println!("{}", is_power_of_two(64));
    println!("{}", next_power_of_two(1000));
    println!("{}", log2_floor(1024));

    // Iterate over all 2^4 non-empty 4-bit subsets.
    let out: Vec<String> = subsets(4).map(|m| format!("{:04b}", m)).collect();
    println!("{}", out.join(" "));

    // Tiny bitset demo: set bits 1, 3, 7, 64, 130, then read them back.
    let mut bs = Bitset::new(200);
    for i in [1usize, 3, 7, 64, 130] {
        bs.set(i);
    }
    let total = bs.count_ones();
    println!("{}", total);
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn popcount_and_trailing_zeros_agree() {
        let x: u64 = 0b1011_0110_u64;
        assert_eq!(popcount(x), 5);
        assert_eq!(trailing_zeros(x), 1);
    }

    #[test]
    fn highest_and_lowest_bit_masks() {
        let x: u64 = 0b1011_0110_u64;
        assert_eq!(highest_bit(x), 0b1000_0000_u64);
        assert_eq!(lowest_bit(x), 0b10_u64);
    }

    #[test]
    fn power_of_two_predicate() {
        assert!(is_power_of_two(1));
        assert!(is_power_of_two(64));
        assert!(!is_power_of_two(0));
        assert!(!is_power_of_two(6));
    }

    #[test]
    fn next_power_and_log2_floor() {
        assert_eq!(next_power_of_two(1), 1);
        assert_eq!(next_power_of_two(2), 2);
        assert_eq!(next_power_of_two(1000), 1024);
        assert_eq!(log2_floor(1024), 10);
    }

    #[test]
    fn subsets_iterator_covers_full_powerset() {
        let k: u32 = 3;
        let collected: Vec<u64> = subsets(k).collect();
        let expected: usize = ((1u64 << k) - 1).try_into().unwrap();
        assert_eq!(collected.len(), expected);
        assert_eq!(collected.iter().sum::<u64>(), (1u64 << k) * ((1u64 << k) - 1) / 2);
    }

    #[test]
    fn bitset_set_get_reset_flip() {
        let mut bs = Bitset::new(200);
        for i in [1usize, 3, 7, 64, 130] {
            bs.set(i);
        }
        assert_eq!(bs.count_ones(), 5);
        assert!(bs.get(64));
        assert!(!bs.get(65));
        bs.flip(64);
        assert!(!bs.get(64));
        bs.flip(64);
        assert!(bs.get(64));
        bs.reset(1);
        assert!(!bs.get(1));
        assert_eq!(bs.count_ones(), 4);
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
