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

struct Node<T> {
    parent: u32,
    jump: u32,
    depth: u32,
    value: Option<T>, // `None` only for the sentinel root
}

/// A persistent queue of `T`.
pub struct PersistentQueue<T> {
    // Node ids start at 1; id 0 is the sentinel root (depth 0, no value). The fields of a node
    // live together so that a walk up the tree touches one cache line per step.
    nodes: Vec<Node<T>>,
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
        PersistentQueue { nodes: vec![Node { parent: 0, jump: 0, depth: 0, value: None }], versions: vec![(0, 0)] }
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
        (self.nodes[tail as usize].depth - pops) as usize
    }

    pub fn is_empty(&self, version: usize) -> bool {
        self.len(version) == 0
    }

    /// The ancestor of `u` at depth `d` (`d <= depth[u]`), in `O(log n)`.
    fn ancestor(&self, mut u: u32, d: u32) -> u32 {
        loop {
            let node = &self.nodes[u as usize];
            if node.depth <= d {
                return u;
            }
            u = if self.nodes[node.jump as usize].depth >= d { node.jump } else { node.parent };
        }
    }

    /// The front element of `version`.
    ///
    /// # Panics
    /// Panics if the version is empty or unknown.
    pub fn front(&self, version: usize) -> T {
        let (tail, pops) = self.check(version);
        assert!(self.len(version) > 0, "front of an empty queue (version {version})");
        let node = self.ancestor(tail, pops + 1);
        self.nodes[node as usize].value.clone().expect("the root holds no value")
    }

    /// A new version: `version` with `x` appended. Returns its id.
    pub fn push_back(&mut self, version: usize, x: T) -> usize {
        let (tail, pops) = self.check(version);
        let id = self.nodes.len() as u32;
        let p = &self.nodes[tail as usize];
        let jp = &self.nodes[p.jump as usize];
        // Skew-binary jump pointers: jump two levels' worth when the last two jumps are equal.
        let jump = if p.depth - jp.depth == jp.depth - self.nodes[jp.jump as usize].depth { jp.jump } else { tail };
        let depth = p.depth + 1;
        self.nodes.push(Node { parent: tail, jump, depth, value: Some(x) });
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
