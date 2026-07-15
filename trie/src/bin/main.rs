use trie::Trie;

/// Driver for the trie example.
///
/// Inserts `"apple"` then `"app"`, and exercises `search` and `starts_with`.
pub fn run() {
    let mut trie = Trie::new();
    trie.insert("apple");
    println!("{}", trie.search("apple"));    // true
    println!("{}", trie.search("app"));      // false
    println!("{}", trie.starts_with("app")); // true
    trie.insert("app");
    println!("{}", trie.search("app"));      // true
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
    fn insert_search_starts_with() {
        let mut trie = Trie::new();
        trie.insert("apple");
        assert!(trie.search("apple"));
        assert!(!trie.search("app"));
        assert!(trie.starts_with("app"));
        trie.insert("app");
        assert!(trie.search("app"));
        assert!(trie.starts_with("ap"));
        assert!(!trie.starts_with("banana"));
    }
}