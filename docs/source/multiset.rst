MultiSet
========

A multiset is a set that allows for multiple occurrences of the same element. It is similar to a list, but it does not allow for duplicate elements.
We will use multiset that has keys in sorted order

Rust Implementation
-------------------

.. code-block:: rust

    use std::clone::Clone;
    use std::iter::FromIterator;
    use std::collections::BTreeMap;


    #[derive(Debug, Default)]
    pub struct MultiSet<T> {
        elems: BTreeMap<T, usize>,
    }

    impl<T: Ord + Clone> MultiSet<T> {
        pub fn new() -> Self {
            MultiSet {
                elems: BTreeMap::new(),
            }
        }

        pub fn insert(&mut self, el: T) {
            *self.elems.entry(el).or_insert(0) += 1;
        }

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

        pub fn remove_all(&mut self, el: &T) -> bool {
            if self.elems.remove(el).is_some() {
                true
            } else {
                false
            }
        }

        pub fn contains(&self, el: &T) -> bool {
            self.elems.contains_key(el)
        }

        pub fn first(&self) -> Option<T> {
            self.elems.keys().next().cloned()
        }

        pub fn last(&self) -> Option<T> {
            self.elems.keys().next_back().cloned()
        }

        pub fn iter(&self) -> impl Iterator<Item = (&T, &usize)> {
            self.elems.iter()
        }

        pub fn into_iter(self) -> impl Iterator<Item = (T, usize)> {
            self.elems.into_iter()
        }

        pub fn iter_mut(&mut self) -> impl Iterator<Item = (&T, &mut usize)> {
            self.elems.iter_mut()
        }

        pub fn count(&self, el: &T) -> usize {
            *self.elems.get(el).unwrap_or(&0)
        }

        pub fn len(&self) -> usize {
            self.elems.values().sum()
        }

        pub fn is_empty(&self) -> bool {
            self.elems.is_empty()
        }

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

