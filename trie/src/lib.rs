//! Trie (prefix tree) for storing and querying strings.
//!
//! Supports three operations, each `O(L)` where `L` is the length of the input:
//! - [`Trie::insert`] - add a word to the trie
//! - [`Trie::search`] - whether a word was previously inserted
//! - [`Trie::starts_with`] - whether any inserted word has the given prefix
//!
//! Children are stored in a `HashMap<char, TrieNode>` for `O(1)` descent.

use std::collections::HashMap;

#[derive(Default)]
struct TrieNode {
    children: HashMap<char, TrieNode>,
    is_end_of_word: bool,
}

/// Trie data structure.
#[derive(Default)]
pub struct Trie {
    root: TrieNode,
}

impl Trie {
    /// Construct an empty trie.
    pub fn new() -> Self {
        Default::default()
    }

    /// Insert `word` into the trie.
    pub fn insert(&mut self, word: &str) {
        let mut curr = &mut self.root;
        for ch in word.chars() {
            curr = curr.children.entry(ch).or_default();
        }
        curr.is_end_of_word = true;
    }

    /// Whether `word` was previously inserted (and not merely a prefix).
    pub fn search(&self, word: &str) -> bool {
        let mut curr = &self.root;
        for ch in word.chars() {
            match curr.children.get(&ch) {
                Some(node) => curr = node,
                None => return false,
            }
        }
        curr.is_end_of_word
    }

    /// Whether any inserted word has `prefix`.
    pub fn starts_with(&self, prefix: &str) -> bool {
        let mut curr = &self.root;
        for ch in prefix.chars() {
            match curr.children.get(&ch) {
                Some(node) => curr = node,
                None => return false,
            }
        }
        true
    }
}