use rolling_hash::{DoubleRollingHash, RollingHash, compute_hash};

/// Run the rolling-hash example.
pub fn run() {
    println!("rolling_hash");

    // Single-modulus rolling hash: slide a window of length 3 over a fixed
    // string and print both the window and the rolling hash value.
    let bytes: &[u8] = b"abcabcabc";
    let mut rh = RollingHash::new();
    for i in 0..3 {
        rh.push(bytes[i]);
    }
    println!("{} {}",
        std::str::from_utf8(&bytes[..3]).unwrap(),
        rh.value()
    );
    for i in 3..bytes.len() {
        rh.pop(bytes[i - 3], 3);
        rh.push(bytes[i]);
        println!("{} {}",
            std::str::from_utf8(&bytes[i - 2..=i]).unwrap(),
            rh.value()
        );
    }

    // Compute the full-string hash for an arbitrary substring, then
    // compare with the rolling updates. Both should agree.
    let full: u64 = compute_hash(&bytes[..3], 131, (1 << 61) - 1);
    println!("{}", full);

    // Double-modulus: print a 2-tuple for "abc".
    let mut drh = DoubleRollingHash::new(131, (1 << 61) - 1, 137, (1 << 61) - 1);
    for &b in &b"abc"[..] {
        drh.push(b);
    }
    let (a, b) = drh.value();
    println!("{} {}", a, b);
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rolling_hash_matches_naive_on_each_window() {
        let bytes: &[u8] = b"abcabcabc";
        let base = 131u64;
        let mod_ = (1u64 << 61) - 1;

        // Sliding window of length 3.
        let mut rh = RollingHash::with_base_mod(base, mod_);
        for i in 0..3 {
            rh.push(bytes[i]);
        }
        // hash for "abc"
        assert_eq!(rh.value(), compute_hash(b"abc", base, mod_));
        for i in 3..bytes.len() {
            rh.pop(bytes[i - 3], 3);
            rh.push(bytes[i]);
            let window = &bytes[i - 2..=i];
            assert_eq!(rh.value(), compute_hash(window, base, mod_),
                "mismatch at window starting at {}", i - 2);
        }
    }

    #[test]
    fn rolling_hash_advances_in_o1() {
        // The value after pushing a byte and immediately after a fresh
        // `compute_hash` of the same single-byte slice must agree.
        let byte = b'x';
        let mut rh = RollingHash::new();
        rh.push(byte);
        assert_eq!(rh.value(), compute_hash(&[byte], 131, (1 << 61) - 1));
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
