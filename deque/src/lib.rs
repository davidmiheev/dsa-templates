//! `Deque<T>`: a growable ring buffer (capacity a power of two, so the wrap-around is a mask)
//! with `O(1)` amortised `push_front` / `push_back` / `pop_front` / `pop_back` and `O(1)`
//! indexing by position from the front. Safe code only: slots are plain `T` values, hence
//! `T: Copy + Default`.
//!
//! Operations with a precondition (`pop`, `front`, `back`, indexing) assert it instead of
//! returning an `Option`.
//!
//! # Example
//!
//! ```
//! use deque::Deque;
//!
//! let mut d = Deque::new();
//! d.push_back(2);
//! d.push_back(3);
//! d.push_front(1);
//! assert_eq!(d[0], 1);
//! assert_eq!(d[2], 3);
//! assert_eq!(d.pop_front(), 1);
//! assert_eq!(d.pop_back(), 3);
//! assert_eq!(d.len(), 1);
//! ```

use std::ops::{Index, IndexMut};

/// Ring-buffer deque. See the [module documentation](self).
#[derive(Clone, Debug)]
pub struct Deque<T> {
    buf: Vec<T>,
    head: usize,
    len: usize,
}

impl<T: Copy + Default> Default for Deque<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Copy + Default> Deque<T> {
    pub fn new() -> Self {
        Self::with_capacity(0)
    }

    /// Room for `n` elements without reallocating.
    pub fn with_capacity(n: usize) -> Self {
        Self { buf: vec![T::default(); n.next_power_of_two().max(4)], head: 0, len: 0 }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    fn mask(&self) -> usize {
        self.buf.len() - 1
    }

    fn grow(&mut self) {
        let mut bigger = vec![T::default(); self.buf.len() * 2];
        for i in 0..self.len {
            bigger[i] = self.buf[(self.head + i) & self.mask()];
        }
        self.buf = bigger;
        self.head = 0;
    }

    pub fn push_back(&mut self, x: T) {
        if self.len == self.buf.len() {
            self.grow();
        }
        let i = (self.head + self.len) & self.mask();
        self.buf[i] = x;
        self.len += 1;
    }

    pub fn push_front(&mut self, x: T) {
        if self.len == self.buf.len() {
            self.grow();
        }
        self.head = (self.head + self.mask()) & self.mask(); // head - 1 modulo the capacity
        self.buf[self.head] = x;
        self.len += 1;
    }

    /// Remove and return the first element.
    ///
    /// # Panics
    /// Panics if the deque is empty.
    pub fn pop_front(&mut self) -> T {
        assert!(self.len > 0, "pop_front on an empty deque");
        let x = self.buf[self.head];
        self.head = (self.head + 1) & self.mask();
        self.len -= 1;
        x
    }

    /// Remove and return the last element.
    ///
    /// # Panics
    /// Panics if the deque is empty.
    pub fn pop_back(&mut self) -> T {
        assert!(self.len > 0, "pop_back on an empty deque");
        self.len -= 1;
        self.buf[(self.head + self.len) & self.mask()]
    }

    /// The first element (panics if empty).
    pub fn front(&self) -> T {
        self[0]
    }

    /// The last element (panics if empty).
    pub fn back(&self) -> T {
        self[self.len - 1]
    }

    /// The elements from front to back.
    pub fn iter(&self) -> impl Iterator<Item = T> + '_ {
        (0..self.len).map(move |i| self[i])
    }
}

impl<T: Copy + Default> Index<usize> for Deque<T> {
    type Output = T;
    fn index(&self, i: usize) -> &T {
        assert!(i < self.len, "index {i} out of range for a deque of {} elements", self.len);
        &self.buf[(self.head + i) & self.mask()]
    }
}

impl<T: Copy + Default> IndexMut<usize> for Deque<T> {
    fn index_mut(&mut self, i: usize) -> &mut T {
        assert!(i < self.len, "index {i} out of range for a deque of {} elements", self.len);
        let mask = self.mask();
        &mut self.buf[(self.head + i) & mask]
    }
}
