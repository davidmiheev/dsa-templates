use suffix_automaton::SuffixAutomaton;

/// Run the suffix automaton example.
pub fn run() {
    println!("suffix_automaton");

    let s = "ababa";
    let sam = SuffixAutomaton::from_str(s);

    println!("contains(\"aba\") = {}", sam.contains("aba"));
    println!("contains(\"bab\") = {}", sam.contains("bab"));
    println!("contains(\"bb\")  = {}", sam.contains("bb"));
    println!("contains(\"\")    = {}", sam.contains(""));
    println!("distinct substrings of \"{}\" = {}", s, sam.distinct_substring_count());
    println!("state count = {}", sam.states.len());
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn substring_membership() {
        let sam = SuffixAutomaton::from_str("ababa");
        assert!(sam.contains("a"));
        assert!(sam.contains("ab"));
        assert!(sam.contains("aba"));
        assert!(sam.contains("bab"));
        assert!(sam.contains("ababa"));
        assert!(!sam.contains("bb"));
        assert!(!sam.contains("abc"));
    }

    #[test]
    fn distinct_substrings_count() {
        // "ababa" has 9 distinct substrings: "", a, b, ab, ba, aba, bab, abab, ababa.
        let sam = SuffixAutomaton::from_str("ababa");
        assert_eq!(sam.distinct_substring_count(), 9);
    }

    #[test]
    fn state_count_bounded_by_2n_minus_1() {
        let s = "aaaaa";
        let sam = SuffixAutomaton::from_str(s);
        assert!(sam.states.len() <= 2 * s.len());
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
