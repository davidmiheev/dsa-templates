//! Suffix array by SA-IS, LCP array, and the substring counts built on them.

fn sa_naive(s: &[u32]) -> Vec<u32> {
    let mut sa: Vec<u32> = (0..s.len() as u32).collect();
    sa.sort_by(|&a, &b| s[a as usize..].cmp(&s[b as usize..]));
    sa
}

const EMPTY: u32 = u32::MAX;

fn sa_is(s: &[u32], upper: usize) -> Vec<u32> {
    let n = s.len();
    if n < 10 {
        return sa_naive(s);
    }
    let mut ls = vec![false; n];
    for i in (0..n - 1).rev() {
        ls[i] = if s[i] == s[i + 1] { ls[i + 1] } else { s[i] < s[i + 1] };
    }
    let mut sum_l = vec![0u32; upper + 2];
    let mut sum_s = vec![0u32; upper + 2];
    for i in 0..n {
        if !ls[i] {
            sum_s[s[i] as usize] += 1;
        } else {
            sum_l[s[i] as usize + 1] += 1;
        }
    }
    for i in 0..=upper {
        sum_s[i] += sum_l[i];
        sum_l[i + 1] += sum_s[i];
    }
    let induce = |lms: &[u32]| -> Vec<u32> {
        let mut sa = vec![EMPTY; n];
        let mut buf = sum_s.clone();
        for &d in lms {
            if d as usize == n {
                continue;
            }
            let c = s[d as usize] as usize;
            sa[buf[c] as usize] = d;
            buf[c] += 1;
        }
        buf.copy_from_slice(&sum_l);
        let c = s[n - 1] as usize;
        sa[buf[c] as usize] = (n - 1) as u32;
        buf[c] += 1;
        for i in 0..n {
            let v = sa[i];
            if v != EMPTY && v >= 1 && !ls[v as usize - 1] {
                let c = s[v as usize - 1] as usize;
                sa[buf[c] as usize] = v - 1;
                buf[c] += 1;
            }
        }
        buf.copy_from_slice(&sum_l);
        for i in (0..n).rev() {
            let v = sa[i];
            if v != EMPTY && v >= 1 && ls[v as usize - 1] {
                let c = s[v as usize - 1] as usize + 1;
                buf[c] -= 1;
                sa[buf[c] as usize] = v - 1;
            }
        }
        sa
    };
    let mut lms_map = vec![EMPTY; n + 1];
    let mut m = 0usize;
    for i in 1..n {
        if !ls[i - 1] && ls[i] {
            lms_map[i] = m as u32;
            m += 1;
        }
    }
    let lms: Vec<u32> = (1..n).filter(|&i| !ls[i - 1] && ls[i]).map(|i| i as u32).collect();
    let sa = induce(&lms);
    if m == 0 {
        return sa;
    }
    let sorted_lms: Vec<u32> = sa.iter().copied().filter(|&v| lms_map[v as usize] != EMPTY).collect();
    let mut rec_s = vec![0u32; m];
    let mut rec_upper = 0usize;
    rec_s[lms_map[sorted_lms[0] as usize] as usize] = 0;
    for i in 1..m {
        let (mut l, mut r) = (sorted_lms[i - 1] as usize, sorted_lms[i] as usize);
        let (ml, mr) = (lms_map[l] as usize, lms_map[r] as usize);
        let end_l = if ml + 1 < m { lms[ml + 1] as usize } else { n };
        let end_r = if mr + 1 < m { lms[mr + 1] as usize } else { n };
        let same = if end_l - l != end_r - r {
            false
        } else {
            while l < end_l {
                if s[l] != s[r] {
                    break;
                }
                l += 1;
                r += 1;
            }
            !(l == n || s[l] != s[r])
        };
        if !same {
            rec_upper += 1;
        }
        rec_s[mr] = rec_upper as u32;
    }
    let rec_sa = sa_is(&rec_s, rec_upper);
    let sorted: Vec<u32> = rec_sa.iter().map(|&i| lms[i as usize]).collect();
    induce(&sorted)
}

/// The suffix array of `s` (entries in `0..=upper`) by SA-IS, `O(n + upper)`.
pub fn suffix_array(s: &[usize], upper: usize) -> Vec<usize> {
    assert!(s.iter().all(|&c| c <= upper) && upper < u32::MAX as usize, "character above `upper`");
    let v: Vec<u32> = s.iter().map(|&c| c as u32).collect();
    sa_is(&v, upper).into_iter().map(|x| x as usize).collect()
}

/// The suffix array of a byte string.
pub fn suffix_array_bytes(s: &[u8]) -> Vec<usize> {
    let v: Vec<u32> = s.iter().map(|&c| c as u32).collect();
    sa_is(&v, 255).into_iter().map(|x| x as usize).collect()
}

/// `lcp[i] = LCP(s[sa[i]..], s[sa[i + 1]..])` for `i < n - 1` (Kasai), `O(n)`.
pub fn lcp_array<T: Ord>(s: &[T], sa: &[usize]) -> Vec<usize> {
    let n = s.len();
    assert_eq!(sa.len(), n, "suffix array length mismatch");
    let mut rank = vec![0usize; n];
    for (i, &p) in sa.iter().enumerate() {
        rank[p] = i;
    }
    let mut lcp = vec![0usize; n.saturating_sub(1)];
    let mut h = 0usize;
    for i in 0..n {
        h = h.saturating_sub(1);
        if rank[i] == 0 {
            continue;
        }
        let j = sa[rank[i] - 1];
        while j + h < n && i + h < n && s[j + h] == s[i + h] {
            h += 1;
        }
        lcp[rank[i] - 1] = h;
    }
    lcp
}

/// The number of distinct non-empty substrings of `s`: `n (n + 1) / 2 - sum lcp`.
pub fn count_distinct_substrings(s: &[u8]) -> u64 {
    let n = s.len() as u64;
    if n == 0 {
        return 0;
    }
    let sa = suffix_array_bytes(s);
    let total: u64 = lcp_array(s, &sa).iter().map(|&x| x as u64).sum();
    n * (n + 1) / 2 - total
}

/// A longest common substring of `s` and `t` as `(a, b, c, d)` with `s[a..b] == t[c..d]`,
/// `(0, 0, 0, 0)` if it is empty. Bytes must be nonzero (0 is the separator).
pub fn longest_common_substring(s: &[u8], t: &[u8]) -> (usize, usize, usize, usize) {
    assert!(s.iter().chain(t).all(|&c| c != 0), "bytes must be nonzero");
    let mut joined = s.to_vec();
    joined.push(0);
    joined.extend_from_slice(t);
    let sa = suffix_array_bytes(&joined);
    let lcp = lcp_array(&joined, &sa);
    let ns = s.len();
    let (mut best, mut at) = (0usize, (0usize, 0usize));
    for i in 0..sa.len() - 1 {
        let (x, y) = (sa[i], sa[i + 1]);
        if (x < ns) != (y < ns) && x != ns && y != ns && lcp[i] > best {
            best = lcp[i];
            at = if x < ns { (x, y - ns - 1) } else { (y, x - ns - 1) };
        }
    }
    (at.0, at.0 + best, at.1, at.1 + best)
}
