use std::collections::HashMap;

#[derive(Default)]
struct TrieNode {
    children: HashMap<char, TrieNode>,
    is_end_of_word: bool,
}

#[derive(Default)]
pub struct Trie {
    root: TrieNode,
}

impl Trie {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn insert(&mut self, word: &str) {
        let mut curr = &mut self.root;
        for ch in word.chars() {
            curr = curr.children.entry(ch).or_default();
        }
        curr.is_end_of_word = true;
    }

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

fn main() {
    let mut trie = Trie::new();
    trie.insert("apple");
    println!("{}", trie.search("apple"));   // true
    println!("{}", trie.search("app"));     // false
    println!("{}", trie.starts_with("app")); // true
    trie.insert("app");
    println!("{}", trie.search("app"));     // true
}
