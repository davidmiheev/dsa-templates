use math::ModBinomial;

/// Driver for the modular binomial example.
pub fn run() {
    let modulo = 1_000_000_007; // Typical modulo 10^9 + 7
    let binom = ModBinomial::new(1000, modulo);

    println!("5 C 2 = {}", binom.ncr(5, 2));   // 10
    println!("10 C 3 = {}", binom.ncr(10, 3)); // 120
    println!("10 C 10 = {}", binom.ncr(10, 10)); // 1
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_smoke() {
        run();
    }

    #[test]
    fn ncr_small_values() {
        let modulo = 1_000_000_007;
        let binom = ModBinomial::new(20, modulo);
        assert_eq!(binom.ncr(5, 2), 10);
        assert_eq!(binom.ncr(10, 3), 120);
        assert_eq!(binom.ncr(10, 10), 1);
        assert_eq!(binom.ncr(10, 0), 1);
        assert_eq!(binom.ncr(10, 11), 0); // r > n
    }

    #[test]
    fn ncr_modulo_alternate_prime() {
        let binom = ModBinomial::new(50, 998_244_353);
        assert_eq!(binom.ncr(4, 2), 6);
        assert_eq!(binom.ncr(20, 10), 184756);
    }
}