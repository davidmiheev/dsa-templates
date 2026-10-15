//! A bucketed sorted list (sorted multiset with rank queries), generic over `T: Ord + Copy`.
//!
//! Elements live in sorted `Vec` chunks ("buckets") of at most `2 * LOAD` elements, plus a
//! contiguous array with each bucket's maximum. This is the structure behind Python's
//! `sortedcontainers.SortedList`: a lookup is a binary search over the maxima followed by a
//! binary search inside one bucket, an insert or remove shifts at most `2 * LOAD` elements in
//! one cache-friendly `memmove`, and rank queries add up the sizes of the buckets before.
//! Everything is `O(log n + sqrt n)`-ish in theory but ~hundreds of nanoseconds in practice,
//! several times faster than a pointer-based balanced tree for `n` up to millions.
//!
//! Duplicates are allowed. Operations:
//! - [`SortedList::insert`] / [`SortedList::remove`] (one occurrence) / [`SortedList::contains`]
//! - [`SortedList::bisect_left`] / [`SortedList::bisect_right`] (ranks)
//! - [`SortedList::get`] (k-th smallest, asserts the index)
//!
//! # Example
//!
//! ```
//! use sorted_list::SortedList;
//!
//! let mut s: SortedList<i64> = [5, 3, 8, 3].into_iter().collect();
//! s.insert(1);
//! assert_eq!(s.get(0), 1);
//! assert_eq!(s.bisect_left(3), 1);  // keys < 3
//! assert_eq!(s.bisect_right(3), 3); // keys <= 3
//! assert!(s.remove(3) && s.contains(3));
//! assert_eq!(s.len(), 4);
//! ```

/// Target bucket size; a bucket is split when it grows past `2 * LOAD`.
const LOAD: usize = 512;

/// Sorted multiset of `T` with `O(log n + sqrt n)` insert, remove, rank and k-th queries.
#[derive(Clone, Debug)]
pub struct SortedList<T> {
    buckets: Vec<Vec<T>>,
    /// `maxes[i]` is the last (largest) element of `buckets[i]`; buckets are never empty.
    maxes: Vec<T>,
    len: usize,
}

impl<T: Ord + Copy> Default for SortedList<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Ord + Copy> FromIterator<T> for SortedList<T> {
    /// Build from any iterator in `O(n log n)` (a sort and one pass).
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut v: Vec<T> = iter.into_iter().collect();
        v.sort_unstable();
        let len = v.len();
        let buckets: Vec<Vec<T>> = v.chunks(LOAD).map(|c| c.to_vec()).collect();
        let maxes = buckets.iter().map(|b| *b.last().unwrap()).collect();
        SortedList { buckets, maxes, len }
    }
}

impl<T: Ord + Copy> SortedList<T> {
    /// Construct an empty sorted list.
    pub fn new() -> Self {
        SortedList { buckets: Vec::new(), maxes: Vec::new(), len: 0 }
    }

    /// Insert a key. Duplicates are allowed.
    pub fn insert(&mut self, key: T) {
        self.len += 1;
        if self.buckets.is_empty() {
            self.buckets.push(vec![key]);
            self.maxes.push(key);
            return;
        }
        // First bucket whose maximum is >= key (or the last bucket for a new global maximum).
        let b = self.maxes.partition_point(|&m| m < key).min(self.buckets.len() - 1);
        let bucket = &mut self.buckets[b];
        let pos = bucket.partition_point(|&x| x < key);
        bucket.insert(pos, key);
        if pos + 1 == bucket.len() {
            self.maxes[b] = key;
        }
        if bucket.len() > 2 * LOAD {
            let right = bucket.split_off(LOAD);
            self.maxes[b] = *self.buckets[b].last().unwrap();
            self.maxes.insert(b + 1, *right.last().unwrap());
            self.buckets.insert(b + 1, right);
        }
    }

    /// Remove one occurrence of `key`. Returns whether it was present.
    pub fn remove(&mut self, key: T) -> bool {
        let b = self.maxes.partition_point(|&m| m < key);
        if b == self.buckets.len() {
            return false;
        }
        let bucket = &mut self.buckets[b];
        let pos = bucket.partition_point(|&x| x < key);
        if bucket[pos] != key {
            return false;
        }
        bucket.remove(pos);
        self.len -= 1;
        if bucket.is_empty() {
            self.buckets.remove(b);
            self.maxes.remove(b);
        } else {
            self.maxes[b] = *bucket.last().unwrap();
        }
        true
    }

    /// Whether `key` occurs at least once.
    pub fn contains(&self, key: T) -> bool {
        let b = self.maxes.partition_point(|&m| m < key);
        b < self.buckets.len() && {
            let bucket = &self.buckets[b];
            bucket[bucket.partition_point(|&x| x < key)] == key
        }
    }

    /// Number of keys strictly less than `key`.
    pub fn bisect_left(&self, key: T) -> usize {
        let b = self.maxes.partition_point(|&m| m < key);
        let before: usize = self.buckets[..b].iter().map(Vec::len).sum();
        match self.buckets.get(b) {
            Some(bucket) => before + bucket.partition_point(|&x| x < key),
            None => before,
        }
    }

    /// Number of keys less than or equal to `key`.
    pub fn bisect_right(&self, key: T) -> usize {
        let b = self.maxes.partition_point(|&m| m <= key);
        let before: usize = self.buckets[..b].iter().map(Vec::len).sum();
        match self.buckets.get(b) {
            Some(bucket) => before + bucket.partition_point(|&x| x <= key),
            None => before,
        }
    }

    /// The element at `index` (0-based) in sorted order.
    ///
    /// # Panics
    /// Panics if `index >= len()`.
    pub fn get(&self, index: usize) -> T {
        assert!(index < self.len, "index {index} out of range for length {}", self.len);
        let mut rest = index;
        for bucket in &self.buckets {
            if rest < bucket.len() {
                return bucket[rest];
            }
            rest -= bucket.len();
        }
        unreachable!("index was checked against len")
    }

    /// All elements in ascending order.
    pub fn iter(&self) -> impl Iterator<Item = T> + '_ {
        self.buckets.iter().flatten().copied()
    }

    /// Number of elements in the list (with multiplicities).
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the list is empty.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}
