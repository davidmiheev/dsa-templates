//! Seaweed ("sticky braid") machinery and longest increasing subsequences of subarrays.
//!
//! * [`steady_ant`] multiplies two permutations in the *unit-Monge* sense (the min-plus
//!   product of their distribution matrices) in `O(n log n)`.
//! * [`range_lis`] answers many queries "length of the longest increasing subsequence of
//!   `a[l..r]`" offline in `O((n + q) log^2 n)`: the whole sequence is turned into a single
//!   seaweed permutation by divide and conquer over positions (every merge is one
//!   [`steady_ant`] product), after which each query is a 2D dominance count.
//!
//! # Example
//!
//! ```
//! use seaweed::range_lis;
//!
//! let a = [3, 1, 2, 5, 4, 6];
//! // LIS of [3,1,2,5,4,6] is 1 2 4 6; of [1,2,5] it is 3; of [5,4] it is 1.
//! assert_eq!(range_lis(&a, &[(0, 6), (1, 4), (3, 5), (2, 2)], true), vec![4, 3, 1, 0]);
//! ```

use fenwick::FenwickTree;

const NONE: u32 = u32::MAX;

/// The sticky product `p ⊡ q` of two permutations of `0..n` (`p[i]` is the image of `i`).
///
/// With the distribution matrix `Σ_p(i, j) = #{ r >= i : p[r] < j }`, the result is the
/// permutation `r` with `Σ_r(i, k) = min_j (Σ_p(i, j) + Σ_q(j, k))`. This is the product of
/// seaweed (sticky braid) permutations; it is associative and has the identity permutation as
/// unit. Runs in `O(n log n)`.
///
/// # Panics
/// Panics if the slices have different lengths. They must be permutations of `0..n`
/// (not checked).
pub fn steady_ant(p: &[u32], q: &[u32]) -> Vec<u32> {
    assert_eq!(p.len(), q.len(), "both permutations must have the same length");
    let n = p.len();
    let mut out = vec![0u32; n];
    let mut scratch = vec![0u32; scratch_len(n)];
    ant(p, q, &mut out, &mut scratch);
    out
}

/// Scratch words `ant` needs for permutations of size `n`.
fn scratch_len(n: usize) -> usize {
    20 * n + 64
}

/// Writes `p ⊡ q` into `out`. `scratch` must hold at least `scratch_len(n)` words.
fn ant(p: &[u32], q: &[u32], out: &mut [u32], scratch: &mut [u32]) {
    let n = p.len();
    match n {
        0 => return,
        1 => {
            out[0] = 0;
            return;
        }
        2 => {
            // The identity is the unit, and swap ⊡ swap = swap.
            out.copy_from_slice(if p[0] == 0 { q } else { p });
            return;
        }
        _ => {}
    }
    let h = n / 2;
    let (rows_lo, s) = scratch.split_at_mut(h);
    let (rows_hi, s) = s.split_at_mut(n - h);
    let (p_lo, s) = s.split_at_mut(h);
    let (p_hi, s) = s.split_at_mut(n - h);
    let (cols_lo, s) = s.split_at_mut(h);
    let (cols_hi, s) = s.split_at_mut(n - h);
    let (rank, s) = s.split_at_mut(n);
    let (q_lo, s) = s.split_at_mut(h);
    let (q_hi, s) = s.split_at_mut(n - h);
    let (r_lo, s) = s.split_at_mut(h);
    let (r_hi, s) = s.split_at_mut(n - h);
    let (lo_col, s) = s.split_at_mut(n);
    let (hi_col, s) = s.split_at_mut(n);
    let (lo_row, s) = s.split_at_mut(n);
    let (hi_row, s) = s.split_at_mut(n);

    // Split the middle index at `h`: p by column value, q by row index.
    let (mut nl, mut nh) = (0, 0);
    for (i, &c) in p.iter().enumerate() {
        if (c as usize) < h {
            rows_lo[nl] = i as u32;
            p_lo[nl] = c;
            nl += 1;
        } else {
            rows_hi[nh] = i as u32;
            p_hi[nh] = c - h as u32;
            nh += 1;
        }
    }
    rank.fill(0);
    for &c in &q[..h] {
        rank[c as usize] = 1; // 1 = belongs to the low half
    }
    let (mut a, mut b) = (0u32, 0u32);
    for (c, r) in rank.iter_mut().enumerate() {
        if *r == 1 {
            *r = a;
            cols_lo[a as usize] = c as u32;
            a += 1;
        } else {
            *r = b;
            cols_hi[b as usize] = c as u32;
            b += 1;
        }
    }
    for (x, &c) in q[..h].iter().enumerate() {
        q_lo[x] = rank[c as usize];
    }
    for (x, &c) in q[h..].iter().enumerate() {
        q_hi[x] = rank[c as usize];
    }
    ant(p_lo, q_lo, r_lo, s);
    ant(p_hi, q_hi, r_hi, s);

    // The two sub-products placed back into the n x n grid (partial permutations), by row and
    // by column.
    lo_col.fill(NONE);
    hi_col.fill(NONE);
    lo_row.fill(NONE);
    hi_row.fill(NONE);
    for (x, &i) in rows_lo.iter().enumerate() {
        let c = cols_lo[r_lo[x] as usize];
        lo_col[i as usize] = c;
        lo_row[c as usize] = i;
    }
    for (x, &i) in rows_hi.iter().enumerate() {
        let c = cols_hi[r_hi[x] as usize];
        hi_col[i as usize] = c;
        hi_row[c as usize] = i;
    }

    // D(i, k) = Σ_hi(i, k) + a(i) - Σ_lo(i, k) - b(k) is non-increasing in i and k; the true
    // product follows the sub-products where D >= 0 (lo) resp. D < 0 (hi), plus fixes on the
    // staircase between the regions. `t` is the smallest k with D(i, k) < 0 (n + 1 if none).
    // dcol(i, k) = D(i, k + 1) - D(i, k); drow(i, k) = D(i + 1, k) - D(i, k).
    let dcol = |i: usize, k: usize| -> i64 {
        if hi_row[k] != NONE {
            (hi_row[k] as usize >= i) as i64 - 1
        } else {
            -((lo_row[k] as usize >= i) as i64)
        }
    };
    let drow = |i: usize, k: usize| -> i64 {
        if hi_col[i] != NONE {
            -(((hi_col[i] as usize) < k) as i64)
        } else {
            ((lo_col[i] as usize) < k) as i64 - 1
        }
    };
    out.fill(NONE);
    let mut t = n + 1;
    let mut hi_seen = 0i64;
    for r in 0..n {
        if hi_col[r] != NONE {
            hi_seen += 1;
        }
        // The staircase column of the next row, and D(r + 1, tn).
        let (mut k, mut e);
        if t == n + 1 {
            e = -hi_seen; // D(r + 1, n)
            k = if e < 0 { n } else { n + 1 };
        } else {
            e = -1 + drow(r, t);
            k = t;
        }
        if k <= n {
            while k > 0 {
                let prev = e - dcol(r + 1, k - 1);
                if prev < 0 {
                    k -= 1;
                    e = prev;
                } else {
                    break;
                }
            }
        }
        let tn = k;
        if lo_col[r] != NONE && (lo_col[r] as usize) + 1 < tn {
            out[r] = lo_col[r];
        }
        if hi_col[r] != NONE && hi_col[r] as usize >= t {
            out[r] = hi_col[r];
        }
        if tn <= n {
            // Cells of the band: columns tn - 1 ..= min(t - 1, n - 1).
            let mut d_prev = 0i64; // D(r + 1, c) for c >= tn
            for c in tn - 1..=(t - 1).min(n - 1) {
                let gc = if c >= tn { d_prev } else { 0 };
                let gc1 = if c + 1 == tn { e } else { d_prev + dcol(r + 1, c) };
                let nz = (lo_col[r] == c as u32) as i64 - (c + 1 == t) as i64 - gc1 + gc;
                if nz == 1 {
                    out[r] = c as u32;
                }
                d_prev = gc1;
            }
        }
        t = tn;
    }
}

/// Kernel of the block of columns `lo..hi` of the grid (values are the rows): the seaweed
/// permutation on `2 * (hi - lo)` strands, and the sorted values of the block.
fn build(perm: &[u32], lo: usize, hi: usize, scratch: &mut [u32]) -> (Vec<u32>, Vec<u32>) {
    if hi - lo == 1 {
        return (vec![perm[lo]], vec![0, 1]);
    }
    let mid = (lo + hi) / 2;
    let (vl, pl) = build(perm, lo, mid, scratch);
    let (vr, pr) = build(perm, mid, hi, scratch);
    let (nl, nr) = (mid - lo, hi - mid);
    let k = vl.len() + vr.len();
    // Merge the value sets; `map_*[row]` = rank of the block's `row`-th value in the union.
    let (mut v, mut map_l, mut map_r) = (Vec::with_capacity(k), Vec::with_capacity(vl.len()), Vec::with_capacity(vr.len()));
    let (mut a, mut b) = (0, 0);
    while a < vl.len() || b < vr.len() {
        if b == vr.len() || (a < vl.len() && vl[a] < vr[b]) {
            map_l.push(v.len() as u32);
            v.push(vl[a]);
            a += 1;
        } else {
            map_r.push(v.len() as u32);
            v.push(vr[b]);
            b += 1;
        }
    }
    let left = lift(&pl, vl.len(), nl, &map_l, k);
    let right = lift(&pr, vr.len(), nr, &map_r, k);
    // (left ⊕ identity) ⊡ (identity ⊕ right)
    let total = k + nl + nr;
    let mut x = Vec::with_capacity(total);
    x.extend_from_slice(&left);
    x.extend((k + nl) as u32..total as u32);
    let mut y = Vec::with_capacity(total);
    y.extend(0..nl as u32);
    y.extend(right.iter().map(|&e| e + nl as u32));
    let mut out = vec![0u32; total];
    ant(&x, &y, &mut out, scratch);
    (v, out)
}

/// Re-index a kernel with `k` compressed rows and `n` columns to `big_k` rows, the compressed
/// row `rho` becoming row `map[rho]`; the new rows are crossed straight through.
fn lift(pi: &[u32], k: usize, n: usize, map: &[u32], big_k: usize) -> Vec<u32> {
    let mut out = vec![NONE; big_k + n];
    let mut used = vec![false; big_k];
    for &m in map {
        used[m as usize] = true;
    }
    let start = |s: usize| if s < k { big_k - 1 - map[k - 1 - s] as usize } else { big_k + (s - k) };
    let end = |e: usize| if e < n { e } else { n + big_k - 1 - map[k - 1 - (e - n)] as usize };
    for s in 0..k + n {
        out[start(s)] = end(pi[s] as usize) as u32;
    }
    for (row, &u) in used.iter().enumerate() {
        if !u {
            out[big_k - 1 - row] = (n + big_k - 1 - row) as u32;
        }
    }
    out
}

/// Longest increasing subsequence of `a[l..r]` for every query `(l, r)`, offline.
///
/// `strict = true` asks for strictly increasing subsequences, `false` for non-decreasing
/// ones. The empty range has answer 0.
///
/// # Panics
/// Panics if a query is not `l <= r <= a.len()` or if `a.len() >= 2^31`.
pub fn range_lis<T: Ord + Copy>(a: &[T], queries: &[(usize, usize)], strict: bool) -> Vec<usize> {
    let n = a.len();
    assert!(n < (1 << 31), "sequence too long");
    for &(l, r) in queries {
        assert!(l <= r && r <= n, "invalid query [{l}, {r}) for length {n}");
    }
    if n == 0 {
        return vec![0; queries.len()];
    }
    // Turn `a` into a permutation: ties are ranked so that they can never be chained.
    let mut order: Vec<u32> = (0..n as u32).collect();
    if strict {
        order.sort_unstable_by(|&i, &j| a[i as usize].cmp(&a[j as usize]).then(j.cmp(&i)));
    } else {
        order.sort_unstable_by(|&i, &j| a[i as usize].cmp(&a[j as usize]).then(i.cmp(&j)));
    }
    let mut perm = vec![0u32; n];
    for (rank, &i) in order.iter().enumerate() {
        perm[i as usize] = rank as u32;
    }
    let mut scratch = vec![0u32; scratch_len(2 * n)];
    let (_, pi) = build(&perm, 0, n, &mut scratch);
    // Top strand `c` (start n + c) that ends at bottom column d < n: the point (c, d).
    let mut col_at = vec![NONE; n]; // by d
    for c in 0..n {
        let d = pi[n + c] as usize;
        if d < n {
            col_at[d] = c as u32;
        }
    }
    // LIS(l, r) = (r - l) - #{ points with c >= l and d < r }.
    let mut by_r: Vec<usize> = (0..queries.len()).collect();
    by_r.sort_unstable_by_key(|&i| queries[i].1);
    let mut fenwick = FenwickTree::<i32>::with_len(n);
    let mut answers = vec![0usize; queries.len()];
    let (mut inserted, mut total) = (0usize, 0usize);
    for qi in by_r {
        let (l, r) = queries[qi];
        while inserted < r {
            if col_at[inserted] != NONE {
                fenwick.update(col_at[inserted] as usize, 1);
                total += 1;
            }
            inserted += 1;
        }
        let before_l = if l == 0 { 0 } else { fenwick.sum_first(l) as usize };
        answers[qi] = (r - l) - (total - before_l);
    }
    answers
}
