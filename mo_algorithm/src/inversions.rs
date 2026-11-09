use crate::{block_size, Query};

/// Sqrt decomposition over `0..v` with `O(sqrt v)` insertion and `O(1)` "how many inserted
/// values are `<= x`".
struct PrefixCounter {
    block: usize,
    /// `before_block[b]` = inserted values in blocks `< b`.
    before_block: Vec<u32>,
    /// `within[x]` = inserted values `y <= x` in the same block as `x`.
    within: Vec<u32>,
}

impl PrefixCounter {
    fn new(v: usize) -> Self {
        let block = ((v as f64).sqrt().ceil() as usize).max(1);
        PrefixCounter { block, before_block: vec![0; v / block + 2], within: vec![0; v] }
    }

    fn insert(&mut self, x: usize) {
        let b = x / self.block;
        for c in &mut self.before_block[b + 1..] {
            *c += 1;
        }
        let end = ((b + 1) * self.block).min(self.within.len());
        for c in &mut self.within[x..end] {
            *c += 1;
        }
    }

    /// Inserted values `<= x`.
    fn count_le(&self, x: usize) -> u32 {
        self.before_block[x / self.block] + self.within[x]
    }

    /// Inserted values `< x`.
    fn count_lt(&self, x: usize) -> u32 {
        if x == 0 {
            0
        } else {
            self.count_le(x - 1)
        }
    }
}

/// One contiguous run of pointer moves: `delta = constant + sign * sum` where `sum` is the
/// total over `j in lo..hi` of the number of elements `a_i` (`i < k`) that are greater than
/// (`greater`) or smaller than `a_j`.
struct Run {
    k: usize,
    lo: usize,
    hi: usize,
    greater: bool,
    sign: i64,
    constant: i64,
}

/// For every query `[lo, hi)`, the number of pairs `lo <= i < j < hi` with `a[i] > a[j]`.
///
/// Uses second-order Mo: the pointer moves are collected as `O(q)` runs, and the sums that
/// each run needs are answered offline with an `O(1)`-query / `O(sqrt n)`-update structure, so
/// the total is `O((n + q) sqrt n)` with a tiny constant (no logarithm).
///
/// ```
/// use mo_algorithm::{range_inversions, Query};
///
/// let a = [3, 1, 2, 2];
/// let q = [Query { lo: 0, hi: 4 }, Query { lo: 1, hi: 4 }, Query { lo: 2, hi: 2 }];
/// assert_eq!(range_inversions(&a, &q), vec![3, 0, 0]);
/// ```
///
/// # Panics
/// Panics if a query is outside `0..=a.len()` or has `lo > hi`.
pub fn range_inversions<T: Ord + Copy>(a: &[T], queries: &[Query]) -> Vec<u64> {
    let n = a.len();
    for q in queries {
        assert!(q.lo <= q.hi && q.hi <= n, "invalid query [{}, {}) for length {n}", q.lo, q.hi);
    }
    // Compress values (equal values share a rank).
    let mut sorted = a.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let v: Vec<usize> = a.iter().map(|x| sorted.partition_point(|y| y < x)).collect();

    // c1[j] = sum_{t < j} #{i < t : v_i > v_t}; c2[j] = sum_{t < j} #{i < t : v_i < v_t}.
    let (mut c1, mut c2) = (vec![0i64; n + 1], vec![0i64; n + 1]);
    {
        let mut pc = PrefixCounter::new(sorted.len().max(1));
        for t in 0..n {
            c1[t + 1] = c1[t] + (t as u32 - pc.count_le(v[t])) as i64;
            c2[t + 1] = c2[t] + pc.count_lt(v[t]) as i64;
            pc.insert(v[t]);
        }
    }

    // Mo order and the runs of pointer moves, in processing order.
    let b = block_size(n, queries.len());
    let mut order: Vec<usize> = (0..queries.len()).collect();
    order.sort_unstable_by_key(|&i| {
        let q = queries[i];
        let blk = q.lo / b;
        (blk, if blk & 1 == 0 { q.hi as i64 } else { -(q.hi as i64) })
    });
    let mut runs: Vec<Run> = Vec::with_capacity(4 * queries.len());
    let mut runs_before_query: Vec<usize> = Vec::with_capacity(queries.len());
    let (mut l, mut r) = (0usize, 0usize);
    for &qi in &order {
        let Query { lo, hi } = queries[qi];
        if lo < l {
            // Add elements lo..l on the left, window [j + 1, r).
            runs.push(Run { k: r, lo, hi: l, greater: false, sign: 1, constant: -(c2[l] - c2[lo]) });
            l = lo;
        }
        if hi > r {
            // Add elements r..hi on the right, window [l, j).
            runs.push(Run { k: l, lo: r, hi, greater: true, sign: -1, constant: c1[hi] - c1[r] });
            r = hi;
        }
        if lo > l {
            // Remove elements l..lo on the left, window [j, r).
            runs.push(Run { k: r, lo: l, hi: lo, greater: false, sign: -1, constant: c2[lo] - c2[l] });
            l = lo;
        }
        if hi < r {
            // Remove elements hi..r on the right, window [l, j + 1).
            runs.push(Run { k: l, lo: hi, hi: r, greater: true, sign: 1, constant: -(c1[r] - c1[hi]) });
            r = hi;
        }
        runs_before_query.push(runs.len());
    }

    // Evaluate every run offline: bucket by `k`, sweep `k` while inserting a[0..k].
    let mut head = vec![0usize; n + 2];
    for run in &runs {
        head[run.k + 1] += 1;
    }
    for i in 0..=n {
        head[i + 1] += head[i];
    }
    let mut by_k = vec![0usize; runs.len()];
    let mut fill = head.clone();
    for (id, run) in runs.iter().enumerate() {
        by_k[fill[run.k]] = id;
        fill[run.k] += 1;
    }
    let mut sums = vec![0i64; runs.len()];
    let mut pc = PrefixCounter::new(sorted.len().max(1));
    for k in 0..=n {
        for &id in &by_k[head[k]..head[k + 1]] {
            let run = &runs[id];
            let mut s = 0u64;
            if run.greater {
                for &x in &v[run.lo..run.hi] {
                    s += (k as u32 - pc.count_le(x)) as u64;
                }
            } else {
                for &x in &v[run.lo..run.hi] {
                    s += pc.count_lt(x) as u64;
                }
            }
            sums[id] = s as i64;
        }
        if k < n {
            pc.insert(v[k]);
        }
    }

    // Replay: the running total at each query boundary is its answer.
    let mut answers = vec![0u64; queries.len()];
    let (mut total, mut next) = (0i64, 0usize);
    for (pos, &qi) in order.iter().enumerate() {
        while next < runs_before_query[pos] {
            total += runs[next].constant + runs[next].sign * sums[next];
            next += 1;
        }
        answers[qi] = total as u64;
    }
    answers
}
