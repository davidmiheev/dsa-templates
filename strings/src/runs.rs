//! Maximal repetitions (runs).

use crate::suffix::{lcp_array, suffix_array};
use sparse_table::SparseTable;

fn min2(a: usize, b: usize) -> usize {
    a.min(b)
}

type MinTable = SparseTable<usize, fn(usize, usize) -> usize>;

/// Longest common extension queries from a suffix array and a sparse table over the LCP array.
struct Lce {
    n: usize,
    rank: Vec<usize>,
    table: Option<MinTable>,
}

impl Lce {
    fn new(s: &[usize], upper: usize) -> (Self, Vec<usize>) {
        let n = s.len();
        let sa = suffix_array(s, upper);
        let mut rank = vec![0usize; n];
        for (i, &p) in sa.iter().enumerate() {
            rank[p] = i;
        }
        let lcp = lcp_array(s, &sa);
        let table = if lcp.is_empty() { None } else { Some(SparseTable::new(&lcp, min2 as fn(usize, usize) -> usize)) };
        (Lce { n, rank, table }, sa)
    }

    /// The longest common prefix of the suffixes at `i` and `j`.
    fn lce(&self, i: usize, j: usize) -> usize {
        if i == j {
            return self.n - i;
        }
        if i >= self.n || j >= self.n {
            return 0;
        }
        let (a, b) = (self.rank[i].min(self.rank[j]), self.rank[i].max(self.rank[j]));
        self.table.as_ref().expect("n >= 2 here").query(a, b)
    }
}

/// All runs `(t, l, r)` of the string: `s[l..r]` has minimal period `t`, `r - l >= 2 t`, and it is
/// maximal, sorted lexicographically. By the runs theorem each run is found from the longest
/// Lyndon words under one of the two alphabet orders, `O(n log n)` here (suffix-array LCE).
#[allow(clippy::needless_range_loop)]
pub fn runs(s: &[u8]) -> Vec<(usize, usize, usize)> {
    let n = s.len();
    let mut out = Vec::new();
    if n < 2 {
        return out;
    }
    let forward: Vec<usize> = s.iter().map(|&c| c as usize).collect();
    let backward: Vec<usize> = s.iter().rev().map(|&c| c as usize).collect();
    let (lce_f, _) = Lce::new(&forward, 255);
    let (lce_b, _) = Lce::new(&backward, 255);
    for order in 0..2 {
        // under the reversed order the suffix comparison uses the complemented alphabet
        let t: Vec<usize> = if order == 0 { forward.clone() } else { forward.iter().map(|&c| 255 - c).collect() };
        let sa = suffix_array(&t, 255);
        let mut rank = vec![0usize; n];
        for (i, &p) in sa.iter().enumerate() {
            rank[p] = i;
        }
        // next smaller suffix gives the longest Lyndon word starting at each position
        let mut nss = vec![n; n];
        let mut stack: Vec<usize> = Vec::new();
        for i in (0..n).rev() {
            while let Some(&top) = stack.last() {
                if rank[top] > rank[i] {
                    stack.pop();
                } else {
                    break;
                }
            }
            nss[i] = stack.last().copied().unwrap_or(n);
            stack.push(i);
        }
        for i in 0..n {
            let j = nss[i];
            let p = j - i;
            let right = if j < n { lce_f.lce(i, j) } else { 0 };
            // longest common suffix of s[..i] and s[..j]
            let left = if i == 0 {
                0
            } else {
                // reversed index of the last character of s[..i] is n - i
                lce_b.lce(n - i, n - j)
            };
            if left + right >= p {
                out.push((p, i - left, j + right));
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}
