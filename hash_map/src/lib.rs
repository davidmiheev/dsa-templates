//! `FastMap<V>`: a hash map from `u64` keys, tuned for speed on competitive-programming inputs.
//!
//! Open addressing with linear probing in a power-of-two table (load factor at most 1/2), a
//! multiplicative (Fibonacci) hash on the high bits, and backward-shift deletion, so there are no
//! tombstones and lookups stay short after any sequence of removals. The hash has a fixed
//! multiplier — fine for judges, not for adversarial inputs.
//!
//! Lookups that may fail are split into `contains` / `get` (asserting presence) / `get_or`
//! instead of returning `Option`.
//!
//! # Example
//!
//! ```
//! use hash_map::{pair_key, FastMap};
//!
//! let mut m = FastMap::new();
//! m.insert(10, "ten");
//! *m.entry(20, "twenty") = "twenty!";
//! assert_eq!(*m.get(10), "ten");
//! assert_eq!(*m.get_or(30, &"none"), "none");
//! assert!(m.remove(10));
//! assert!(!m.contains(10));
//! assert_eq!(m.len(), 1);
//! // Undirected edge keys:
//! assert_eq!(pair_key(5, 2), pair_key(2, 5));
//! ```

/// A hash map with `u64` keys.
#[derive(Clone, Debug)]
pub struct FastMap<V> {
    slots: Vec<Option<(u64, V)>>,
    len: usize,
    shift: u32,
}

/// Key for the unordered pair `{a, b}` of 32-bit values.
pub fn pair_key(a: usize, b: usize) -> u64 {
    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
    ((lo as u64) << 32) | hi as u64
}

impl<V> Default for FastMap<V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<V> FastMap<V> {
    pub fn new() -> Self {
        Self::with_capacity(0)
    }

    /// Room for `n` entries without growing.
    pub fn with_capacity(n: usize) -> Self {
        let size = (2 * n).next_power_of_two().max(8);
        let mut slots = Vec::new();
        slots.resize_with(size, || None);
        Self { slots, len: 0, shift: 64 - size.trailing_zeros() }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    fn home(&self, key: u64) -> usize {
        (key.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> self.shift) as usize
    }

    /// Index of `key`'s slot, or of the empty slot where it would go.
    fn probe(&self, key: u64) -> usize {
        let mask = self.slots.len() - 1;
        let mut i = self.home(key);
        loop {
            match &self.slots[i] {
                Some((k, _)) if *k != key => i = (i + 1) & mask,
                _ => return i,
            }
        }
    }

    fn grow(&mut self) {
        let size = self.slots.len() * 2;
        let mut bigger = Vec::new();
        bigger.resize_with(size, || None);
        let old = std::mem::replace(&mut self.slots, bigger);
        self.shift = 64 - size.trailing_zeros();
        for (k, v) in old.into_iter().flatten() {
            let i = self.probe(k);
            self.slots[i] = Some((k, v));
        }
    }

    pub fn contains(&self, key: u64) -> bool {
        self.slots[self.probe(key)].is_some()
    }

    /// Insert or overwrite.
    pub fn insert(&mut self, key: u64, value: V) {
        if 2 * (self.len + 1) > self.slots.len() {
            self.grow();
        }
        let i = self.probe(key);
        if self.slots[i].is_none() {
            self.len += 1;
        }
        self.slots[i] = Some((key, value));
    }

    /// The value of `key`.
    ///
    /// # Panics
    /// Panics if the key is absent.
    pub fn get(&self, key: u64) -> &V {
        match &self.slots[self.probe(key)] {
            Some((_, v)) => v,
            None => panic!("key {key} is not in the map"),
        }
    }

    /// The value of `key`, or `default` if absent.
    pub fn get_or<'a>(&'a self, key: u64, default: &'a V) -> &'a V {
        match &self.slots[self.probe(key)] {
            Some((_, v)) => v,
            None => default,
        }
    }

    /// Mutable access to the value of `key`, inserting `default` first if absent.
    pub fn entry(&mut self, key: u64, default: V) -> &mut V {
        if 2 * (self.len + 1) > self.slots.len() {
            self.grow();
        }
        let i = self.probe(key);
        if self.slots[i].is_none() {
            self.slots[i] = Some((key, default));
            self.len += 1;
        }
        &mut self.slots[i].as_mut().unwrap().1
    }

    /// Remove `key`; returns whether it was present.
    pub fn remove(&mut self, key: u64) -> bool {
        let mask = self.slots.len() - 1;
        let mut hole = self.probe(key);
        if self.slots[hole].is_none() {
            return false;
        }
        self.slots[hole] = None;
        self.len -= 1;
        // Backward shift: pull later entries of the cluster into the hole when allowed.
        let mut j = hole;
        loop {
            j = (j + 1) & mask;
            let Some((k, _)) = &self.slots[j] else { break };
            let home = self.home(*k);
            // the entry may move to `hole` iff its home is not cyclically in (hole, j]
            let stays = if hole <= j { hole < home && home <= j } else { hole < home || home <= j };
            if !stays {
                self.slots[hole] = self.slots[j].take();
                hole = j;
            }
        }
        true
    }

    /// Iterate over `(key, &value)` in unspecified order.
    pub fn iter(&self) -> impl Iterator<Item = (u64, &V)> {
        self.slots.iter().flatten().map(|(k, v)| (*k, v))
    }
}
