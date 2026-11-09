use crate::{block_size, Query};

/// Mo's algorithm without removals ("rollback Mo"), for states that can only be extended.
///
/// `add(state, i)` adds element `i` to the window; `remove(state, i)` undoes the most recent
/// still-active `add(state, i)` (it only has to restore the data, not any aggregate that
/// `snapshot` captured). `snapshot` takes a copy of the aggregate (e.g. the current best
/// value) and `restore` puts it back after a batch of undone adds; the aggregate of the empty
/// window is the snapshot taken before the first `add`. `visit(state, q)` is called
/// when the window equals query `q`. Total work is `O((n + q) sqrt n)` calls.
///
/// ```
/// use mo_algorithm::{solve_rollback, Query};
///
/// // Maximum of the window: adding is trivial, removing is not.
/// let a = [3, 1, 4, 1, 5, 9, 2, 6];
/// let queries = [Query { lo: 0, hi: 3 }, Query { lo: 2, hi: 8 }, Query { lo: 5, hi: 6 }];
/// let mut best = 0;
/// let mut answers = vec![0; queries.len()];
/// solve_rollback(
///     a.len(),
///     &queries,
///     &mut best,
///     |best, i| *best = (*best).max(a[i]),
///     |_, _| {}, // nothing to undo besides the snapshot
///     |best| *best,
///     |best, old| *best = old,
///     |best, q| answers[q] = *best,
/// );
/// assert_eq!(answers, vec![4, 9, 9]);
/// ```
#[allow(clippy::too_many_arguments)]
pub fn solve_rollback<State, Snap: Clone, FAdd, FRemove, FSnap, FRestore, FVisit>(
    n: usize,
    queries: &[Query],
    state: &mut State,
    mut add: FAdd,
    mut remove: FRemove,
    mut snapshot: FSnap,
    mut restore: FRestore,
    mut visit: FVisit,
) where
    FAdd: FnMut(&mut State, usize),
    FRemove: FnMut(&mut State, usize),
    FSnap: FnMut(&State) -> Snap,
    FRestore: FnMut(&mut State, Snap),
    FVisit: FnMut(&mut State, usize),
{
    let empty = snapshot(state);
    let b = block_size(n, queries.len());
    let blocks = n / b + 1;
    // Bucket queries by the block of `lo`; inside a block sort by `hi`.
    let mut order: Vec<usize> = (0..queries.len()).collect();
    order.sort_unstable_by_key(|&i| (queries[i].lo / b, queries[i].hi));
    let mut idx = 0;
    for block in 0..blocks {
        let end = (block + 1) * b; // elements `end..` belong to later blocks
        // The permanent part of the window is `end..r`.
        let mut r = end;
        while idx < order.len() && queries[order[idx]].lo / b == block {
            let qi = order[idx];
            idx += 1;
            let Query { lo, hi } = queries[qi];
            if hi <= end {
                // The whole query lies inside the block: brute force it.
                let snap = snapshot(state);
                for i in lo..hi {
                    add(state, i);
                }
                visit(state, qi);
                for i in lo..hi {
                    remove(state, i);
                }
                restore(state, snap);
                continue;
            }
            while r < hi {
                add(state, r);
                r += 1;
            }
            let snap = snapshot(state);
            for i in (lo..end).rev() {
                add(state, i);
            }
            visit(state, qi);
            for i in lo..end {
                remove(state, i);
            }
            restore(state, snap);
        }
        // Drop the permanent part before the next block starts.
        for i in end..r {
            remove(state, i);
        }
        restore(state, empty.clone());
        if idx >= order.len() {
            break;
        }
    }
}
