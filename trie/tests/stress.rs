//! Stress test: `Trie` against a `HashSet<String>` (search) and a linear scan (prefixes).

use std::collections::HashSet;
use trie::Trie;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const ALPHABET: [char; 5] = ['a', 'b', 'c', 'é', '日'];

fn word(rng: &mut Rng, max_len: usize, letters: usize) -> String {
    (0..rng.below(max_len + 1)).map(|_| ALPHABET[rng.below(letters)]).collect()
}

#[test]
fn matches_hash_set() {
    let mut st = Stress::new("trie/naive");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let letters = rng.range(1, 5) as usize; // tiny alphabets force shared prefixes
        let mut trie = Trie::new();
        let mut words: HashSet<String> = HashSet::new();
        st.case();
        assert!(trie.starts_with(""), "empty prefix always matches");
        assert!(!trie.search(""), "empty trie holds no empty word");
        for step in 0..120 {
            let w = word(&mut rng, 6, letters);
            match rng.below(3) {
                0 => {
                    trie.insert(&w);
                    words.insert(w);
                }
                1 => assert_eq!(trie.search(&w), words.contains(&w), "search({w:?}) case {case} step {step}"),
                _ => assert_eq!(
                    trie.starts_with(&w),
                    w.is_empty() || words.iter().any(|x| x.starts_with(&w)),
                    "starts_with({w:?}) case {case} step {step}"
                ),
            }
            st.ops(1);
        }
        for w in &words {
            assert!(trie.search(w), "inserted word {w:?} case {case}");
            for cut in w.char_indices().map(|(i, _)| i) {
                assert!(trie.starts_with(&w[..cut]), "prefix of {w:?} case {case}");
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn many_long_words() {
    let mut st = Stress::new("trie/large");
    let mut rng = Rng::new(3);
    let mut trie = Trie::new();
    let mut words: HashSet<String> = HashSet::new();
    st.case();
    for _ in 0..20_000 * scale() {
        let w = word(&mut rng, 30, 3);
        trie.insert(&w);
        words.insert(w);
    }
    for _ in 0..20_000 * scale() {
        let w = word(&mut rng, 30, 3);
        assert_eq!(trie.search(&w), words.contains(&w));
        st.ops(1);
    }
    st.done();
}
