//! A queue and a deque that keep the fold of a [`Monoid`] over all their elements, with
//! amortized `O(1)` operations ("sliding window aggregation").
//!
//! Both are built from stacks that store, next to every element, the aggregate of the stack up
//! to it; the fold of the whole structure is the product of the two stack tops. When one stack
//! runs empty the other is split (queue: moved over; deque: split in half), which amortizes to
//! constant time. The operation does not need to be commutative or invertible: the fold is
//! always the product of the elements in sequence order.
//!
//! * [`QueueAggregation`]: `push_back`, `pop_front`, `fold`.
//! * [`DequeAggregation`]: `push_front`, `push_back`, `pop_front`, `pop_back`, `fold`.
//!
//! # Example
//!
//! ```
//! use algebra::Monoid;
//! use queue_aggregation::QueueAggregation;
//!
//! /// Concatenation of digits as a decimal number (non-commutative): (value, number of digits).
//! enum Digits {}
//! impl Monoid for Digits {
//!     type S = (u64, u32);
//!     fn identity() -> (u64, u32) { (0, 0) }
//!     fn op(a: &(u64, u32), b: &(u64, u32)) -> (u64, u32) { (a.0 * 10u64.pow(b.1) + b.0, a.1 + b.1) }
//! }
//!
//! let mut q = QueueAggregation::<Digits>::new();
//! for d in [1, 2, 3] { q.push_back((d, 1)); }
//! assert_eq!(q.fold().0, 123);
//! q.pop_front();
//! assert_eq!(q.fold().0, 23);
//! ```

use algebra::Monoid;

/// Queue with an `O(1)` amortized fold over all elements.
pub struct QueueAggregation<M: Monoid> {
    /// Front stack, top = oldest element; each entry is `(element, fold from it to the stack bottom)`.
    front: Vec<(M::S, M::S)>,
    /// Back stack in push order, plus the fold over all of it.
    back: Vec<M::S>,
    back_fold: M::S,
}

impl<M: Monoid> Default for QueueAggregation<M> {
    fn default() -> Self {
        Self::new()
    }
}

impl<M: Monoid> QueueAggregation<M> {
    pub fn new() -> Self {
        Self { front: Vec::new(), back: Vec::new(), back_fold: M::identity() }
    }

    pub fn len(&self) -> usize {
        self.front.len() + self.back.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Append `x` at the back.
    pub fn push_back(&mut self, x: M::S) {
        self.back_fold = M::op(&self.back_fold, &x);
        self.back.push(x);
    }

    /// Remove and return the oldest element.
    ///
    /// # Panics
    /// Panics if the queue is empty.
    pub fn pop_front(&mut self) -> M::S {
        assert!(!self.is_empty(), "pop_front on an empty queue");
        if self.front.is_empty() {
            // Newest first, so that the oldest ends on top of the front stack.
            let mut fold = M::identity();
            while let Some(x) = self.back.pop() {
                fold = M::op(&x, &fold);
                self.front.push((x, fold.clone()));
            }
            self.back_fold = M::identity();
        }
        self.front.pop().unwrap().0
    }

    /// The product of all elements in order (the identity if empty).
    pub fn fold(&self) -> M::S {
        match self.front.last() {
            Some((_, f)) => M::op(f, &self.back_fold),
            None => self.back_fold.clone(),
        }
    }
}

/// Deque with an `O(1)` amortized fold over all elements.
pub struct DequeAggregation<M: Monoid> {
    /// Front stack, top = front of the deque; entry = `(element, fold from it towards the middle)`.
    front: Vec<(M::S, M::S)>,
    /// Back stack, top = back of the deque; entry = `(element, fold from the middle to it)`.
    back: Vec<(M::S, M::S)>,
}

impl<M: Monoid> Default for DequeAggregation<M> {
    fn default() -> Self {
        Self::new()
    }
}

impl<M: Monoid> DequeAggregation<M> {
    pub fn new() -> Self {
        Self { front: Vec::new(), back: Vec::new() }
    }

    pub fn len(&self) -> usize {
        self.front.len() + self.back.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn front_fold(&self) -> M::S {
        self.front.last().map_or_else(M::identity, |e| e.1.clone())
    }

    fn back_fold(&self) -> M::S {
        self.back.last().map_or_else(M::identity, |e| e.1.clone())
    }

    pub fn push_front(&mut self, x: M::S) {
        let fold = M::op(&x, &self.front_fold());
        self.front.push((x, fold));
    }

    pub fn push_back(&mut self, x: M::S) {
        let fold = M::op(&self.back_fold(), &x);
        self.back.push((x, fold));
    }

    /// Rebuild both stacks so that `items` (in deque order) is split in half, the first half on
    /// the front stack.
    fn rebalance(&mut self, items: Vec<M::S>) {
        let mid = (items.len() + 1) / 2;
        self.front.clear();
        self.back.clear();
        for x in items[..mid].iter().rev() {
            self.push_front(x.clone());
        }
        for x in &items[mid..] {
            self.push_back(x.clone());
        }
    }

    /// Remove and return the front element.
    ///
    /// # Panics
    /// Panics if the deque is empty.
    pub fn pop_front(&mut self) -> M::S {
        assert!(!self.is_empty(), "pop_front on an empty deque");
        if self.front.is_empty() {
            let items: Vec<M::S> = self.back.drain(..).map(|e| e.0).collect();
            self.rebalance(items);
        }
        self.front.pop().unwrap().0
    }

    /// Remove and return the back element.
    ///
    /// # Panics
    /// Panics if the deque is empty.
    pub fn pop_back(&mut self) -> M::S {
        assert!(!self.is_empty(), "pop_back on an empty deque");
        if self.back.is_empty() {
            let items: Vec<M::S> = self.front.drain(..).rev().map(|e| e.0).collect(); // deque order
            // Give the back stack the *second* half: split so that the back half is the larger one.
            let mid = items.len() / 2;
            self.front.clear();
            self.back.clear();
            for x in items[..mid].iter().rev() {
                self.push_front(x.clone());
            }
            for x in &items[mid..] {
                self.push_back(x.clone());
            }
        }
        self.back.pop().unwrap().0
    }

    /// The product of all elements in order (the identity if empty).
    pub fn fold(&self) -> M::S {
        M::op(&self.front_fold(), &self.back_fold())
    }
}
