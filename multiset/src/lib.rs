//! A multiset that allows multiple occurrences of the same element.
//!
//! Backed by a [`std::collections::BTreeMap`], so elements are kept in sorted
//! order by key. `BTreeMap` is preferred over `HashMap` here to avoid
//! collision blowups and to allow `O(log n)` `first` / `last` queries.
//!
//! All operations are `O(log n)`; `len` and `is_empty` are `O(n)` (they sum
//! all counts) but can be cached if needed.

use std::collections::BTreeMap;

/// A multiset of values of type `T` (any `Ord + Clone`).
#[derive(Debug, Default)]
pub struct MultiSet<T> {
    elems: BTreeMap<T, usize>,
}

impl<T: Ord + Clone> MultiSet<T> {
    /// Construct an empty multiset.
    pub fn new() -> Self {
        MultiSet {
            elems: BTreeMap::new(),
        }
    }

    /// Insert one more occurrence of `el`.
    pub fn insert(&mut self, el: T) {
        *self.elems.entry(el).or_insert(0) += 1;
    }

    /// Remove one occurrence of `el`. Returns `false` if `el` was not present.
    pub fn remove(&mut self, el: &T) -> bool {
        if let Some(count) = self.elems.get_mut(el) {
            if *count > 1 {
                *count -= 1;
            } else {
                self.elems.remove(el);
            }
            true
        } else {
            false
        }
    }

    /// Remove all occurrences of `el`. Returns `false` if `el` was not present.
    pub fn remove_all(&mut self, el: &T) -> bool {
        self.elems.remove(el).is_some()
    }

    /// Whether `el` is present at least once.
    pub fn contains(&self, el: &T) -> bool {
        self.elems.contains_key(el)
    }

    /// Smallest element (by `Ord`), or `None` if the multiset is empty.
    pub fn first(&self) -> Option<T> {
        self.elems.keys().next().cloned()
    }

    /// Largest element (by `Ord`), or `None` if the multiset is empty.
    pub fn last(&self) -> Option<T> {
        self.elems.keys().next_back().cloned()
    }

    /// Iterate over `(element, count)` pairs in ascending order.
    pub fn iter(&self) -> impl Iterator<Item = (&T, &usize)> {
        self.elems.iter()
    }

    /// Consume the multiset and iterate over `(element, count)` pairs in
    /// ascending order.
    pub fn into_iter(self) -> impl Iterator<Item = (T, usize)> {
        self.elems.into_iter()
    }

    /// Iterate mutably over `(element, count)` pairs in ascending order.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (&T, &mut usize)> {
        self.elems.iter_mut()
    }

    /// Number of occurrences of `el` (0 if absent).
    pub fn count(&self, el: &T) -> usize {
        *self.elems.get(el).unwrap_or(&0)
    }

    /// Total number of elements (with multiplicities).
    pub fn len(&self) -> usize {
        self.elems.values().sum()
    }

    /// Whether the multiset has no elements.
    pub fn is_empty(&self) -> bool {
        self.elems.is_empty()
    }

    /// Remove all elements.
    pub fn clear(&mut self) {
        self.elems.clear();
    }
}

impl<T: Ord + Clone> FromIterator<T> for MultiSet<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut multiset = MultiSet::new();
        for el in iter {
            multiset.insert(el);
        }
        multiset
    }
}