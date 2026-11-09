//! Longest increasing subsequence in `O(n log n)` (patience sorting with back pointers).
//!
//! # Example
//!
//! ```
//! use subsequence::longest_increasing_subsequence;
//!
//! assert_eq!(longest_increasing_subsequence(&[2, 1, 4, 3, 5], true), vec![1, 3, 4]); // values 1, 3, 5 (one of several)
//! assert_eq!(longest_increasing_subsequence(&[2, 2, 2], false), vec![0, 1, 2]);
//! ```

/// The indices of a longest increasing subsequence of `a` (strictly increasing if `strict`,
/// non-decreasing otherwise).
pub fn longest_increasing_subsequence<T: Ord>(a: &[T], strict: bool) -> Vec<usize> {
    let n = a.len();
    // tails[l]: index of the smallest possible last element of an increasing subsequence of
    // length l + 1
    let mut tails: Vec<usize> = Vec::new();
    let mut prev = vec![usize::MAX; n];
    for i in 0..n {
        let pos = tails.partition_point(|&t| if strict { a[t] < a[i] } else { a[t] <= a[i] });
        if pos > 0 {
            prev[i] = tails[pos - 1];
        }
        if pos == tails.len() {
            tails.push(i);
        } else {
            tails[pos] = i;
        }
    }
    let mut out = Vec::with_capacity(tails.len());
    let mut cur = tails.last().copied().unwrap_or(usize::MAX);
    while cur != usize::MAX {
        out.push(cur);
        cur = prev[cur];
    }
    out.reverse();
    out
}
