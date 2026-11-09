//! The permutation tree (common interval decomposition tree) of a permutation.
//!
//! A *common interval* of a permutation `P` is a range of positions `[l, r]` whose values form
//! a range of consecutive integers. The *strong* common intervals (those that do not overlap
//! any other common interval) form a tree under inclusion: the permutation tree. Its leaves are
//! the single positions. An inner node is **linear** when any run of consecutive children is
//! again a common interval (then its children are monotone: increasing or decreasing blocks),
//! and **prime** otherwise (then no union of two or more but not all children is a common
//! interval).
//!
//! [`PermutationTree::new`] builds the tree in `O(n log n)` with a monotone stack of nodes and a
//! lazy segment tree that finds the nearest `l` with `max(P[l..=r]) - min(P[l..=r]) = r - l`.
//!
//! Node ids: the leaf of position `i` is node `i`; inner nodes are numbered after the leaves.
//!
//! # Example
//!
//! ```
//! use permutation_tree::{Kind, PermutationTree};
//!
//! let t = PermutationTree::new(&[2, 0, 3, 1]); // the smallest prime permutation (2413)
//! assert_eq!(t.kind[t.root], Kind::Prime);
//! assert_eq!(t.children(t.root).len(), 4);
//!
//! let t = PermutationTree::new(&[0, 1, 2]);
//! assert_eq!(t.kind[t.root], Kind::LinearIncreasing);
//! assert_eq!(t.children(t.root).len(), 3);
//! ```

use algebra::{MapMonoid, Monoid};
use lazy_segment_tree::LazySegmentTree;

/// The type of a node.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// A single position.
    Leaf,
    /// A linear node whose children have increasing values from left to right.
    LinearIncreasing,
    /// A linear node whose children have decreasing values from left to right.
    LinearDecreasing,
    /// A prime node.
    Prime,
}

/// The permutation tree. See the [module documentation](self).
#[derive(Debug, Clone)]
pub struct PermutationTree {
    /// Parent of each node ([`PermutationTree::NONE`] for the root).
    pub parent: Vec<usize>,
    /// Left end (position) of each node's interval.
    pub left: Vec<usize>,
    /// Right end (position, inclusive) of each node's interval.
    pub right: Vec<usize>,
    /// Type of each node.
    pub kind: Vec<Kind>,
    /// The root node.
    pub root: usize,
    child_start: Vec<usize>,
    child_list: Vec<usize>,
}

const NIL: u32 = u32::MAX;

enum MinAdd {}
impl Monoid for MinAdd {
    type S = i64;
    fn identity() -> i64 {
        i64::MAX / 2
    }
    fn op(a: &i64, b: &i64) -> i64 {
        *a.min(b)
    }
}
enum RangeAdd {}
impl MapMonoid for RangeAdd {
    type M = MinAdd;
    type F = i64;
    fn id_map() -> i64 {
        0
    }
    fn apply(f: &i64, x: &i64) -> i64 {
        x + f
    }
    fn compose(f: &i64, g: &i64) -> i64 {
        f + g
    }
}

impl PermutationTree {
    /// The parent of the root.
    pub const NONE: usize = usize::MAX;

    /// Build the tree of `p`, a permutation of `0..p.len()`.
    ///
    /// # Panics
    /// Panics if `p` is empty or not a permutation of `0..p.len()`.
    pub fn new(p: &[usize]) -> Self {
        let n = p.len();
        assert!(n > 0, "the permutation must not be empty");
        assert!(n < u32::MAX as usize / 2, "permutation too long");
        let mut seen = vec![false; n];
        for &x in p {
            assert!(x < n && !seen[x], "the input is not a permutation of 0..{n}");
            seen[x] = true;
        }
        // Node storage during the build: children as singly linked lists.
        let cap = 2 * n;
        let mut l = vec![0u32; cap];
        let mut r = vec![0u32; cap];
        let mut mn = vec![0u32; cap];
        let mut mx = vec![0u32; cap];
        let mut kind = vec![Kind::Leaf; cap];
        let mut first = vec![NIL; cap];
        let mut last = vec![NIL; cap];
        let mut next = vec![NIL; cap];
        for (i, &x) in p.iter().enumerate() {
            l[i] = i as u32;
            r[i] = i as u32;
            mn[i] = x as u32;
            mx[i] = x as u32;
        }
        let mut count = n;
        // g(l) = max(P[l..=i]) - min(P[l..=i]) + l >= i, with equality exactly when [l, i] is a
        // common interval.
        let init: Vec<i64> = (0..n as i64).collect();
        let mut seg = LazySegmentTree::<RangeAdd>::from_slice(&init);
        let mut max_stack: Vec<(usize, i64)> = Vec::new(); // (first l, value)
        let mut min_stack: Vec<(usize, i64)> = Vec::new();
        let mut stack: Vec<u32> = Vec::new();
        for i in 0..n {
            let v = p[i] as i64;
            // Maintain max(P[l..=i]) and min(P[l..=i]) per l with range adds on g.
            let mut new_start = i;
            while let Some(&(s, val)) = max_stack.last() {
                if val < v {
                    // segment [s, new_start) currently has max `val`
                    seg.apply_range(s, new_start, v - val);
                    new_start = s;
                    max_stack.pop();
                } else {
                    break;
                }
            }
            max_stack.push((new_start, v));
            let mut new_start = i;
            while let Some(&(s, val)) = min_stack.last() {
                if val > v {
                    seg.apply_range(s, new_start, val - v);
                    new_start = s;
                    min_stack.pop();
                } else {
                    break;
                }
            }
            min_stack.push((new_start, v));

            // Merge the new leaf into the stack of strong intervals.
            let mut cur = i as u32;
            loop {
                if let Some(&top) = stack.last() {
                    let t = top as usize;
                    let c = cur as usize;
                    let inc = mx[t] + 1 == mn[c];
                    let dec = mx[c] + 1 == mn[t];
                    if inc || dec {
                        stack.pop();
                        let dir = if inc { Kind::LinearIncreasing } else { Kind::LinearDecreasing };
                        let res = if kind[t] == dir && kind[c] == dir {
                            // Concatenate the children lists.
                            next[last[t] as usize] = first[c];
                            last[t] = last[c];
                            top
                        } else if kind[t] == dir {
                            next[last[t] as usize] = cur;
                            last[t] = cur;
                            top
                        } else if kind[c] == dir {
                            next[t] = first[c];
                            first[c] = top;
                            cur
                        } else {
                            let node = count as u32;
                            count += 1;
                            kind[node as usize] = dir;
                            first[node as usize] = top;
                            next[t] = cur;
                            last[node as usize] = cur;
                            node
                        };
                        let rs = res as usize;
                        l[rs] = l[t];
                        r[rs] = r[c];
                        mn[rs] = mn[t].min(mn[c]);
                        mx[rs] = mx[t].max(mx[c]);
                        cur = res;
                        continue;
                    }
                }
                // Is there a common interval [l*, i] reaching left of `cur`? Take the nearest one
                // (the rightmost such l*), which is the lowest prime node above `cur`; larger
                // ones are found by the following iterations.
                let target = i as i64;
                let cur_left = l[cur as usize] as usize;
                let lp = seg.min_left(cur_left, |m| *m > target);
                if lp > 0 {
                    let l_star = lp - 1;
                    let node = count as u32;
                    count += 1;
                    kind[node as usize] = Kind::Prime;
                    // Children: the stack nodes starting at or after l*, then `cur`.
                    let mut popped: Vec<u32> = Vec::new();
                    while let Some(&top) = stack.last() {
                        if l[top as usize] as usize >= l_star {
                            popped.push(top);
                            stack.pop();
                        } else {
                            break;
                        }
                    }
                    debug_assert_eq!(l[*popped.last().expect("a node starts at l*") as usize] as usize, l_star);
                    popped.reverse();
                    popped.push(cur);
                    for w in popped.windows(2) {
                        next[w[0] as usize] = w[1];
                    }
                    let (a, z) = (popped[0], *popped.last().unwrap());
                    first[node as usize] = a;
                    last[node as usize] = z;
                    let (mut lo, mut hi) = (u32::MAX, 0);
                    for &c in &popped {
                        lo = lo.min(mn[c as usize]);
                        hi = hi.max(mx[c as usize]);
                    }
                    l[node as usize] = l[a as usize];
                    r[node as usize] = r[z as usize];
                    mn[node as usize] = lo;
                    mx[node as usize] = hi;
                    cur = node;
                    continue;
                }
                break;
            }
            stack.push(cur);
        }
        assert_eq!(stack.len(), 1, "the whole permutation is always one common interval");
        let root = stack[0] as usize;
        // Keep only the nodes reachable from the root (a linear node that was absorbed into
        // another is unreachable). Leaves keep their ids `0..n`; inner nodes are renumbered from
        // `n` in BFS order.
        let mut order = vec![root];
        let mut parent = vec![Self::NONE; count];
        let mut head = 0;
        while head < order.len() {
            let v = order[head];
            head += 1;
            let mut c = first[v];
            while c != NIL {
                parent[c as usize] = v;
                order.push(c as usize);
                c = next[c as usize];
            }
        }
        let mut new_id = vec![usize::MAX; count];
        for (i, id) in new_id.iter_mut().enumerate().take(n) {
            *id = i;
        }
        let mut m = n;
        for &v in &order {
            if v >= n {
                new_id[v] = m;
                m += 1;
            }
        }
        let mut parent2 = vec![Self::NONE; m];
        let (mut left, mut right, mut kind2) = (vec![0usize; m], vec![0usize; m], vec![Kind::Leaf; m]);
        for &v in &order {
            let id = new_id[v];
            left[id] = l[v] as usize;
            right[id] = r[v] as usize;
            kind2[id] = kind[v];
            parent2[id] = if parent[v] == Self::NONE { Self::NONE } else { new_id[parent[v]] };
        }
        // CSR children in left-to-right order.
        let mut start = vec![0usize; m + 1];
        for id in 0..m {
            if parent2[id] != Self::NONE {
                start[parent2[id] + 1] += 1;
            }
        }
        for id in 0..m {
            start[id + 1] += start[id];
        }
        let mut fill = start.clone();
        let mut list = vec![0usize; start[m]];
        // Iterating by increasing left end puts siblings in left-to-right order.
        let mut by_left: Vec<usize> = (0..m).filter(|&id| parent2[id] != Self::NONE).collect();
        by_left.sort_unstable_by_key(|&id| (left[id], std::cmp::Reverse(right[id])));
        for id in by_left {
            let p = parent2[id];
            list[fill[p]] = id;
            fill[p] += 1;
        }
        PermutationTree { parent: parent2, left, right, kind: kind2, root: new_id[root], child_start: start, child_list: list }
    }

    /// Number of nodes (leaves included).
    pub fn len(&self) -> usize {
        self.parent.len()
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    /// The children of `v`, from left to right.
    pub fn children(&self, v: usize) -> &[usize] {
        &self.child_list[self.child_start[v]..self.child_start[v + 1]]
    }

    /// Whether `v` is a linear node in the sense of the problem statement (leaves count as
    /// linear).
    pub fn is_linear(&self, v: usize) -> bool {
        self.kind[v] != Kind::Prime
    }
}
