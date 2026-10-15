//! Binary trie over `BITS`-bit unsigned integers (a multiset), for xor queries.
//!
//! Keys are stored most significant bit first, so a walk that prefers the child whose bit makes
//! the xor with `x` smallest (or largest) finds `min(k ^ x)` (or `max(k ^ x)`) over all keys `k`
//! in `O(BITS)`. Every node counts the keys below it, so removal just decrements counts and
//! queries skip empty branches. Memory is at most `BITS` nodes per distinct inserted key.
//!
//! # Example
//!
//! ```
//! use binary_trie::BinaryTrie;
//!
//! let mut t = BinaryTrie::<30>::new();
//! t.insert(5);
//! t.insert(12);
//! assert_eq!(t.min_xor(4), 1); // 5 ^ 4
//! assert_eq!(t.max_xor(4), 12 ^ 4);
//! assert!(t.remove(5) && !t.remove(5));
//! assert_eq!(t.min_xor(4), 12 ^ 4);
//! ```

/// A multiset of integers below `2^BITS` with `O(BITS)` insert, remove and xor queries.
#[derive(Clone, Debug)]
pub struct BinaryTrie<const BITS: usize> {
    /// `children[node] = [child for bit 0, child for bit 1]`; 0 means "none" (node 0 is the root).
    children: Vec<[u32; 2]>,
    /// Number of keys stored below each node.
    count: Vec<u32>,
}

impl<const BITS: usize> Default for BinaryTrie<BITS> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const BITS: usize> BinaryTrie<BITS> {
    /// An empty trie.
    pub fn new() -> Self {
        assert!(BITS >= 1 && BITS <= 64, "BITS must be in 1..=64");
        Self { children: vec![[0, 0]], count: vec![0] }
    }

    #[inline]
    fn bit(x: u64, level: usize) -> usize {
        ((x >> (BITS - 1 - level)) & 1) as usize
    }

    fn check(x: u64) {
        assert!(BITS == 64 || x < (1u64 << BITS), "key {x} does not fit in {BITS} bits");
    }

    /// Number of keys (with multiplicity).
    pub fn len(&self) -> usize {
        self.count[0] as usize
    }

    pub fn is_empty(&self) -> bool {
        self.count[0] == 0
    }

    /// Insert one occurrence of `x`.
    ///
    /// # Panics
    /// Panics if `x >= 2^BITS`.
    pub fn insert(&mut self, x: u64) {
        Self::check(x);
        let mut node = 0;
        self.count[0] += 1;
        for level in 0..BITS {
            let b = Self::bit(x, level);
            if self.children[node][b] == 0 {
                self.children.push([0, 0]);
                self.count.push(0);
                self.children[node][b] = (self.children.len() - 1) as u32;
            }
            node = self.children[node][b] as usize;
            self.count[node] += 1;
        }
    }

    /// Number of occurrences of `x`.
    pub fn count(&self, x: u64) -> usize {
        if BITS < 64 && x >= (1u64 << BITS) {
            return 0;
        }
        let mut node = 0;
        for level in 0..BITS {
            let next = self.children[node][Self::bit(x, level)] as usize;
            if next == 0 {
                return 0;
            }
            node = next;
        }
        self.count[node] as usize
    }

    pub fn contains(&self, x: u64) -> bool {
        self.count(x) > 0
    }

    /// Remove one occurrence of `x`; returns whether it was present.
    pub fn remove(&mut self, x: u64) -> bool {
        if self.count(x) == 0 {
            return false;
        }
        let mut node = 0;
        self.count[0] -= 1;
        for level in 0..BITS {
            node = self.children[node][Self::bit(x, level)] as usize;
            self.count[node] -= 1;
        }
        true
    }

    /// `min(k ^ x)` over the stored keys `k`.
    ///
    /// # Panics
    /// Panics if the trie is empty.
    pub fn min_xor(&self, x: u64) -> u64 {
        self.walk(x, false)
    }

    /// `max(k ^ x)` over the stored keys `k`.
    ///
    /// # Panics
    /// Panics if the trie is empty.
    pub fn max_xor(&self, x: u64) -> u64 {
        self.walk(x, true)
    }

    fn walk(&self, x: u64, maximize: bool) -> u64 {
        assert!(!self.is_empty(), "xor query on an empty trie");
        let (mut node, mut result) = (0usize, 0u64);
        for level in 0..BITS {
            let want = Self::bit(x, level) ^ maximize as usize; // child that gives the preferred xor bit
            let preferred = self.children[node][want] as usize;
            let (b, next) = if preferred != 0 && self.count[preferred] > 0 {
                (want, preferred)
            } else {
                (want ^ 1, self.children[node][want ^ 1] as usize)
            };
            result |= ((Self::bit(x, level) ^ b) as u64) << (BITS - 1 - level);
            node = next;
        }
        result
    }
}
