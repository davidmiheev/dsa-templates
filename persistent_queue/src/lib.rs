//! A fully persistent queue.
//!
//! Every [`PersistentQueue::push_back`] and [`PersistentQueue::pop_front`] returns the id of a *new*
//! version and leaves the old one untouched, so any version can be extended again, branching
//! into a tree of histories.
//!
//! All pushed elements form a tree (each pushed node's parent is the tail of the version it was
//! pushed onto). A version is just `(tail, pops)`: its elements are the nodes on the path from
//! the root to `tail`, without the first `pops`. The front is therefore a *level ancestor*
//! query, answered in `O(log n)` with skew-binary jump pointers that cost one extra word per
//! node. Everything is online (no need to know future operations) and uses `O(1)` extra space
//! per version.
//!
//! # Example
//!
//! ```
//! use persistent_queue::PersistentQueue;
//!
//! let mut q = PersistentQueue::new();
//! let v1 = q.push_back(PersistentQueue::<i32>::EMPTY, 10);
//! let v2 = q.push_back(v1, 20);
//! let (front, v3) = q.pop_front(v2); // v2 is unchanged
//! assert_eq!((front, q.len(v3)), (10, 1));
//! let v4 = q.push_back(v1, 99); // branch from the older version
//! assert_eq!(q.front(v4), 10);
//! assert_eq!(q.len(v2), 2);
//! ```

/// A persistent queue of `T`.
pub struct PersistentQueue<T> {
    // Node ids start at 1; id 0 is the sentinel root (depth 0, no value).
    parent: Vec<u32>,
    jump: Vec<u32>,
    depth: Vec<u32>,
    value: Vec<T>, // value[id - 1]
    /// Version `v` = (tail node, number of pops).
    versions: Vec<(u32, u32)>,
}

impl<T: Clone> Default for PersistentQueue<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Clone> PersistentQueue<T> {
    /// Version id of the empty queue (it exists from the start).
    pub const EMPTY: usize = 0;

    pub fn new() -> Self {
        PersistentQueue { parent: vec![0], jump: vec![0], depth: vec![0], value: Vec::new(), versions: vec![(0, 0)] }
    }

    /// Number of versions created so far (including the empty one).
    pub fn num_versions(&self) -> usize {
        self.versions.len()
    }

    fn check(&self, version: usize) -> (u32, u32) {
        assert!(version < self.versions.len(), "unknown version {version}");
        self.versions[version]
    }

    /// Number of elements in `version`.
    pub fn len(&self, version: usize) -> usize {
        let (tail, pops) = self.check(version);
        (self.depth[tail as usize] - pops) as usize
    }

    pub fn is_empty(&self, version: usize) -> bool {
        self.len(version) == 0
    }

    /// The ancestor of `u` at depth `d` (`d <= depth[u]`), in `O(log n)`.
    fn ancestor(&self, mut u: u32, d: u32) -> u32 {
        while self.depth[u as usize] > d {
            let j = self.jump[u as usize];
            u = if self.depth[j as usize] >= d { j } else { self.parent[u as usize] };
        }
        u
    }

    /// The front element of `version`.
    ///
    /// # Panics
    /// Panics if the version is empty or unknown.
    pub fn front(&self, version: usize) -> T {
        let (tail, pops) = self.check(version);
        assert!(self.len(version) > 0, "front of an empty queue (version {version})");
        let node = self.ancestor(tail, pops + 1);
        self.value[node as usize - 1].clone()
    }

    /// A new version: `version` with `x` appended. Returns its id.
    pub fn push_back(&mut self, version: usize, x: T) -> usize {
        let (tail, pops) = self.check(version);
        let id = self.parent.len() as u32;
        let (p, d) = (tail as usize, self.depth[tail as usize] + 1);
        let jp = self.jump[p] as usize;
        // Skew-binary jump pointers: jump two levels' worth when the last two jumps are equal.
        let jump = if self.depth[p] - self.depth[jp] == self.depth[jp] - self.depth[self.jump[jp] as usize] {
            self.jump[jp]
        } else {
            tail
        };
        self.parent.push(tail);
        self.jump.push(jump);
        self.depth.push(d);
        self.value.push(x);
        self.versions.push((id, pops));
        self.versions.len() - 1
    }

    /// A new version without the front element of `version`; returns `(front, new version)`.
    ///
    /// # Panics
    /// Panics if the version is empty or unknown.
    pub fn pop_front(&mut self, version: usize) -> (T, usize) {
        let front = self.front(version);
        let (tail, pops) = self.versions[version];
        self.versions.push((tail, pops + 1));
        (front, self.versions.len() - 1)
    }
}
