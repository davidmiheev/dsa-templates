//! Mo's algorithm for offline range queries.
//!
//! Mo's algorithm solves queries of the form "given a value array `a`, and
//! a list of range queries `[lo, hi)`, answer each query using only
//! `O(|hi - lo|)` moves between consecutive queries (with `O(1)` amortized
//! per move) once the queries have been reordered.
//!
//! The classic instance is counting distinct elements in a range. This
//! crate ships [`count_distinct_in_range`], which handles that case
//! directly:
//!
//! ```no_run
//! use mo_algorithm::{count_distinct_in_range, Query};
//!
//! let a = vec![1, 1, 2, 1, 3];
//! let queries = vec![Query { lo: 0, hi: 3 }, Query { lo: 2, hi: 5 }];
//! let answers = count_distinct_in_range(&a, &queries);
//! assert_eq!(answers, vec![2, 3]);
//! ```
//!
//! For more elaborate query types build your own solver using
//! [`solve_blocked`] (add and remove), or [`solve_rollback`] when the state can only be
//! extended cheaply (modes, maxima, ...). [`range_inversions`] counts inversions in a range with
//! "second-order" Mo in `O((n + q) sqrt n)`, with `O(1)` work per pointer move.

mod inversions;
mod rollback;

pub use inversions::range_inversions;
pub use rollback::solve_rollback;

/// Block size that balances the two pointers: about `n / sqrt(q)` (at least 1).
pub(crate) fn block_size(n: usize, q: usize) -> usize {
    ((n as f64) / (q.max(1) as f64).sqrt()).ceil().max(1.0) as usize
}

/// A single query on a half-open interval `[lo, hi)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Query {
    /// Inclusive lower bound of the query range.
    pub lo: usize,
    /// Exclusive upper bound of the query range.
    pub hi: usize,
}

impl Query {
    fn block(&self, block_size: usize) -> usize {
        self.lo / block_size
    }
}

/// Run a generic Mo's algorithm over `n` elements answering `queries`.
///
/// All moving window bookkeeping is done in this function; the caller
/// supplies a single `state` value passed to each callback. Interior
/// mutability (e.g. `RefCell`) lets one state object serve all three
/// callbacks if you need shared mutable bookkeeping.
pub fn solve_blocked<State, FAdd, FRemove, FVisit>(
    n: usize,
    queries: &[Query],
    state: &mut State,
    mut add: FAdd,
    mut remove: FRemove,
    mut visit: FVisit,
) where
    State: ?Sized,
    FAdd: FnMut(&mut State, usize),
    FRemove: FnMut(&mut State, usize),
    FVisit: FnMut(&mut State, usize),
{
    if queries.is_empty() {
        return;
    }
    let block_size = block_size(n, queries.len());

    // Sort queries by Mo's ordering.
    let mut order: Vec<usize> = (0..queries.len()).collect();
    order.sort_by_key(|&i| {
        let q = queries[i];
        let b = q.block(block_size);
        (b, if b & 1 == 0 { q.hi as i64 } else { -(q.hi as i64) })
    });

    let mut cur_lo = queries[order[0]].lo;
    let mut cur_hi = queries[order[0]].lo;

    for &qi in &order {
        let q = queries[qi];
        // Expand / contract window to cover [q.lo, q.hi).
        while cur_lo > q.lo {
            cur_lo -= 1;
            add(state, cur_lo);
        }
        while cur_hi < q.hi {
            add(state, cur_hi);
            cur_hi += 1;
        }
        while cur_lo < q.lo {
            remove(state, cur_lo);
            cur_lo += 1;
        }
        while cur_hi > q.hi {
            cur_hi -= 1;
            remove(state, cur_hi);
        }
        visit(state, qi);
    }
}

/// Solve "number of distinct values in range" queries using Mo's
/// algorithm.
///
/// `a` is the value array; `queries` are the half-open intervals.
/// Returns one answer per query.
pub fn count_distinct_in_range<T: Eq + std::hash::Hash + Copy>(
    a: &[T],
    queries: &[Query],
) -> Vec<usize> {
    // Map every value to a small id once, so the hot loop is plain array indexing.
    let mut ids: std::collections::HashMap<T, u32> = std::collections::HashMap::new();
    let id: Vec<u32> = a.iter().map(|&x| {
        let next = ids.len() as u32;
        *ids.entry(x).or_insert(next)
    }).collect();
    struct State {
        freq: Vec<u32>,
        distinct: usize,
    }

    let mut state = State { freq: vec![0; ids.len()], distinct: 0 };
    let mut answers: Vec<usize> = vec![0; queries.len()];

    solve_blocked(
        a.len(),
        queries,
        &mut state,
        |st, idx| {
            let f = &mut st.freq[id[idx] as usize];
            st.distinct += (*f == 0) as usize;
            *f += 1;
        },
        |st, idx| {
            let f = &mut st.freq[id[idx] as usize];
            *f -= 1;
            st.distinct -= (*f == 0) as usize;
        },
        |st, qi| {
            answers[qi] = st.distinct;
        },
    );

    answers
}
