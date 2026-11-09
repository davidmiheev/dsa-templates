//! Offline 3D dominance sums by CDQ divide and conquer.
//!
//! Given a sequence of events in *time order* — weighted point additions and point queries —
//! [`dominance_sums`] answers every query: the total weight of the points added **earlier**
//! whose `x` and `y` are both **strictly smaller** than the query's. Rectangle sums with point
//! additions (four queries per rectangle), point reads after rectangle additions (four
//! additions per rectangle), and similar "time, x, y" problems all reduce to this.
//!
//! It runs in `O(n log^2 n)` with sequential memory access, which makes it several times faster
//! than [`Fenwick2D`](crate::Fenwick2D) when all events are known in advance.
//!
//! # Example
//!
//! ```
//! use fenwick::{dominance_sums, Event};
//!
//! let events = [
//!     Event::Add { x: 1, y: 1, w: 10 },
//!     Event::Query { x: 5, y: 5 },            // sees the point at (1, 1)
//!     Event::Add { x: 2, y: 3, w: 5 },
//!     Event::Query { x: 3, y: 4 },            // sees both
//!     Event::Query { x: 2, y: 4 },            // only (1, 1): x must be strictly smaller
//!     Event::Add { x: 0, y: 0, w: 1000 },     // later than every query
//! ];
//! assert_eq!(dominance_sums::<i32, i32, i64>(&events), vec![10, 15, 10]);
//! ```

use crate::FenwickTree;
use std::ops::{AddAssign, Sub};

/// One step of the timeline given to [`dominance_sums`].
#[derive(Clone, Copy, Debug)]
pub enum Event<X, Y, W> {
    /// A point of weight `w` appears at `(x, y)`.
    Add { x: X, y: Y, w: W },
    /// Ask for the weight of the earlier points with `px < x` and `py < y`.
    Query { x: X, y: Y },
}

#[derive(Clone, Copy)]
struct Item<X, W> {
    x: X,
    y: u32,
    /// The weight of an addition, or the output slot of a query.
    w: W,
    slot: u32,
}

const ADD: u32 = u32::MAX;

/// The answers of all [`Event::Query`] events, in the order they appear in `events`.
///
/// # Panics
/// Panics if there are `2^32` or more events.
pub fn dominance_sums<X, Y, W>(events: &[Event<X, Y, W>]) -> Vec<W>
where
    X: Ord + Copy,
    Y: Ord + Copy,
    W: Copy + Default + AddAssign + Sub<Output = W>,
{
    assert!(events.len() < u32::MAX as usize, "too many events");
    let mut all_ys: Vec<Y> = events
        .iter()
        .map(|e| match *e {
            Event::Add { y, .. } | Event::Query { y, .. } => y,
        })
        .collect();
    all_ys.sort_unstable();
    all_ys.dedup();
    let rank = |y: Y| all_ys.partition_point(|v| *v < y) as u32;
    let mut queries = 0u32;
    let mut items: Vec<Item<X, W>> = events
        .iter()
        .map(|e| match *e {
            Event::Add { x, y, w } => Item { x, y: rank(y), w, slot: ADD },
            Event::Query { x, y } => {
                queries += 1;
                Item { x, y: rank(y), w: W::default(), slot: queries - 1 }
            }
        })
        .collect();
    let mut answers = vec![W::default(); queries as usize];
    let mut fenwick = FenwickTree::<W>::with_len(all_ys.len());
    let mut buf = Vec::with_capacity(items.len());
    cdq(&mut items, &mut buf, &mut fenwick, &mut answers);
    answers
}

/// Sorts `items` by `x` while adding, for each query in the right half, the contributions of
/// the additions in the left half.
fn cdq<X, W>(items: &mut [Item<X, W>], buf: &mut Vec<Item<X, W>>, fenwick: &mut FenwickTree<W>, answers: &mut [W])
where
    X: Ord + Copy,
    W: Copy + Default + AddAssign + Sub<Output = W>,
{
    let n = items.len();
    if n <= 1 {
        return;
    }
    let mid = n / 2;
    {
        let (left, right) = items.split_at_mut(mid);
        cdq(left, buf, fenwick, answers);
        cdq(right, buf, fenwick, answers);
    }
    // Both halves are sorted by x; feed the left additions in x order as the queries advance.
    let mut i = 0;
    for j in mid..n {
        if items[j].slot == ADD {
            continue;
        }
        while i < mid && items[i].x < items[j].x {
            if items[i].slot == ADD {
                fenwick.update(items[i].y as usize, items[i].w);
            }
            i += 1;
        }
        let y = items[j].y as usize;
        if y > 0 {
            answers[items[j].slot as usize] += fenwick.sum_first(y);
        }
    }
    for item in items.iter().take(i) {
        if item.slot == ADD {
            // Undo with `delta = 0 - w` so that only `Sub` is needed.
            fenwick.update(item.y as usize, W::default() - item.w);
        }
    }
    buf.clear();
    let (mut a, mut b) = (0, mid);
    while a < mid && b < n {
        if items[a].x <= items[b].x {
            buf.push(items[a]);
            a += 1;
        } else {
            buf.push(items[b]);
            b += 1;
        }
    }
    buf.extend_from_slice(&items[a..mid]);
    buf.extend_from_slice(&items[b..n]);
    items.copy_from_slice(buf);
}
