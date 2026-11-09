//! Offline dynamic connectivity: a graph whose edges appear and disappear over the times
//! `0..times`, with something to evaluate at every time.
//!
//! Every edge is alive during a half-open interval of times. The intervals are hung on the nodes
//! of a segment tree over time, and a depth-first walk applies the edges of each node to a
//! rollback union-find on the way down and undoes them on the way up, so at the leaf of time `t`
//! the structure contains exactly the edges alive at `t`. Total `O((m + times) log times *
//! log n)`.
//!
//! * [`EdgeIntervals`] turns "toggle edge `{u, v}` at time `t`" events into lifetimes.
//! * [`offline_dynamic_connectivity`] runs the walk over any [`Rollback`] structure —
//!   [`union_find::RollbackUnionFind`] for plain connectivity, or
//!   [`union_find::RollbackSumUnionFind`] when components carry (changing) sums.
//!
//! # Example
//!
//! ```
//! use offline_connectivity::{offline_dynamic_connectivity, EdgeIntervals};
//! use union_find::RollbackUnionFind;
//!
//! // time 0: add 0-1, time 1: add 1-2, time 2: remove 0-1; ask "is 0 ~ 2?" at every time.
//! let mut ev = EdgeIntervals::new();
//! ev.toggle(0, 0, 1);
//! ev.toggle(1, 1, 2);
//! ev.toggle(2, 0, 1);
//! let mut answers = vec![false; 3];
//! let mut uf = RollbackUnionFind::new(3);
//! offline_dynamic_connectivity(3, &ev.finish(3), &mut uf, |t, uf| answers[t] = uf.same(0, 2));
//! assert_eq!(answers, vec![false, true, false]);
//! ```

use hash_map::{pair_key, FastMap};
use union_find::{RollbackSumUnionFind, RollbackUnionFind};

/// A union-find that can be rolled back.
pub trait Rollback {
    /// Merge the components of `u` and `v`.
    fn union(&mut self, u: usize, v: usize);
    /// A marker for the current state.
    fn snapshot(&self) -> usize;
    /// Return to a state previously taken with [`Rollback::snapshot`].
    fn rollback(&mut self, snapshot: usize);
}

impl Rollback for RollbackUnionFind {
    fn union(&mut self, u: usize, v: usize) {
        RollbackUnionFind::union(self, u, v);
    }
    fn snapshot(&self) -> usize {
        RollbackUnionFind::snapshot(self)
    }
    fn rollback(&mut self, snapshot: usize) {
        RollbackUnionFind::rollback(self, snapshot)
    }
}

impl<T: Copy + std::ops::Add<Output = T> + std::ops::Sub<Output = T>> Rollback for RollbackSumUnionFind<T> {
    fn union(&mut self, u: usize, v: usize) {
        RollbackSumUnionFind::union(self, u, v);
    }
    fn snapshot(&self) -> usize {
        RollbackSumUnionFind::snapshot(self)
    }
    fn rollback(&mut self, snapshot: usize) {
        RollbackSumUnionFind::rollback(self, snapshot)
    }
}

/// One edge `{u, v}` alive for the times `from..to`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Lifetime {
    pub u: usize,
    pub v: usize,
    pub from: usize,
    pub to: usize,
}

/// Collects edge lifetimes from toggle events given in non-decreasing time order.
#[derive(Default)]
pub struct EdgeIntervals {
    open: FastMap<(usize, usize, usize)>, // pair key -> (u, v, start)
    done: Vec<Lifetime>,
}

impl EdgeIntervals {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add the edge `{u, v}` at time `t` if it is absent, remove it if it is present. The edge is
    /// alive from the time it is added (inclusive) to the time it is removed (exclusive).
    pub fn toggle(&mut self, t: usize, u: usize, v: usize) {
        let key = pair_key(u, v);
        if self.open.contains(key) {
            let (a, b, from) = *self.open.get(key);
            self.open.remove(key);
            self.done.push(Lifetime { u: a, v: b, from, to: t });
        } else {
            self.open.insert(key, (u, v, t));
        }
    }

    /// Close all remaining edges at `times` and return every lifetime (empty ones dropped).
    pub fn finish(mut self, times: usize) -> Vec<Lifetime> {
        let rest: Vec<(usize, usize, usize)> = self.open.iter().map(|(_, &x)| x).collect();
        for (u, v, from) in rest {
            self.done.push(Lifetime { u, v, from, to: times });
        }
        self.done.retain(|l| l.from < l.to);
        self.done
    }
}

/// Walk the times `0..times` in order; at each time `t`, `dsu` holds exactly the edges alive at
/// `t` and `visit(t, &mut dsu)` is called (it may read or update `dsu`'s own data, but must leave
/// the partition as it found it).
pub fn offline_dynamic_connectivity<D: Rollback>(
    times: usize,
    lifetimes: &[Lifetime],
    dsu: &mut D,
    mut visit: impl FnMut(usize, &mut D),
) {
    if times == 0 {
        return;
    }
    let size = times.next_power_of_two();
    let mut node_edges: Vec<Vec<(usize, usize)>> = vec![Vec::new(); 2 * size];
    for l in lifetimes {
        let (mut lo, mut hi) = (l.from + size, l.to.min(times) + size);
        while lo < hi {
            if lo & 1 == 1 {
                node_edges[lo].push((l.u, l.v));
                lo += 1;
            }
            if hi & 1 == 1 {
                hi -= 1;
                node_edges[hi].push((l.u, l.v));
            }
            lo >>= 1;
            hi >>= 1;
        }
    }
    // Depth-first walk with explicit enter/leave events.
    let mut stack: Vec<(usize, Option<usize>)> = vec![(1, None)];
    while let Some((node, undo_to)) = stack.pop() {
        if let Some(snap) = undo_to {
            dsu.rollback(snap);
            continue;
        }
        let snap = dsu.snapshot();
        for &(u, v) in &node_edges[node] {
            dsu.union(u, v);
        }
        if node >= size {
            let t = node - size;
            if t < times {
                visit(t, dsu);
            }
            dsu.rollback(snap);
        } else {
            stack.push((node, Some(snap)));
            // left child must be visited first
            let (l, r) = (2 * node, 2 * node + 1);
            stack.push((r, None));
            stack.push((l, None));
        }
    }
}
