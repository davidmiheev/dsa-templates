use std::clone::Clone;
use std::iter::FromIterator;
use std::collections::BTreeMap;

/// A multiset is a collection that allows multiple occurrences of the same element.
/// This implementation uses a `BTreeMap` to store the elements and their counts 
/// in sorted order by key. The elements must implement the `Ord` trait to maintain order
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

fn main() {
    let mut ms = MultiSet::from_iter(vec![1, 2, 2, 3, 4, 4, 4]);
    
    ms.insert(1);
    ms.insert(2);
    ms.insert(2);
    ms.insert(3);
    
    println!("Multiset: {:?}", ms);
    println!("Count of 2: {}", ms.count(&2));

    if let Some(first) = ms.first() {
        println!("First element: {}", first);
    }
    
    if let Some(last) = ms.last() {
        println!("Last element: {}", last);
    }
    
    ms.remove(&2);
    println!("After removing one occurrence of 2: {:?}", ms);
    
    ms.remove_all(&4);
    println!("After removing all occurrences of 3: {:?}", ms);
    
    if let Some(first) = ms.first() {
        println!("First element: {}", first);
    }
    
    if let Some(last) = ms.last() {
        println!("Last element: {}", last);
    }
}