//! Z-function, Manacher, Lyndon factorisation, wildcard matching.

use convolution::convolution_u64;

/// `z[i]` = longest common prefix of `s` and `s[i..]` (`z[0] = n`), `O(n)`.
pub fn z_algorithm<T: PartialEq>(s: &[T]) -> Vec<usize> {
    let n = s.len();
    let mut z = vec![0usize; n];
    if n == 0 {
        return z;
    }
    z[0] = n;
    let (mut l, mut r) = (0usize, 0usize);
    for i in 1..n {
        let mut k = if i < r { z[i - l].min(r - i) } else { 0 };
        while i + k < n && s[k] == s[i + k] {
            k += 1;
        }
        z[i] = k;
        if i + k > r {
            (l, r) = (i, i + k);
        }
    }
    z
}

/// The length of the longest palindrome around each of the `2n - 1` centres, in the order
/// character 0, gap 0|1, character 1, gap 1|2, ..., character `n - 1`; `O(n)`.
pub fn manacher<T: PartialEq>(s: &[T]) -> Vec<usize> {
    let n = s.len();
    if n == 0 {
        return Vec::new();
    }
    // radii over the string with separators: t[2k + 1] = s[k], t[even] = separator
    let m = 2 * n + 1;
    let mut p = vec![0usize; m];
    let (mut c, mut r) = (0usize, 0usize); // centre and right edge (exclusive) of the rightmost palindrome
    let at = |i: usize| -> Option<&T> {
        if i % 2 == 1 {
            Some(&s[i / 2])
        } else {
            None
        }
    };
    for i in 0..m {
        let mut k = if i < r { p[2 * c - i].min(r - i) } else { 0 };
        while i > k && i + k + 1 < m && at(i - k - 1) == at(i + k + 1) {
            k += 1;
        }
        p[i] = k;
        if i + k > r {
            (c, r) = (i, i + k);
        }
    }
    p[1..m - 1].to_vec()
}

/// The Lyndon factorisation by Duval's algorithm: the boundaries `0 = a_0 < a_1 < ... < a_k = n`
/// of the factors `s[a_i..a_{i+1}]`, which are Lyndon words in non-increasing order.
pub fn lyndon_factorization<T: Ord>(s: &[T]) -> Vec<usize> {
    let n = s.len();
    let mut bounds = vec![0usize];
    let mut i = 0;
    while i < n {
        let (mut j, mut k) = (i + 1, i);
        while j < n && s[k] <= s[j] {
            k = if s[k] < s[j] { i } else { k + 1 };
            j += 1;
        }
        let period = j - k;
        while i <= k {
            i += period;
            bounds.push(i);
        }
    }
    if n == 0 {
        bounds.clear();
        bounds.push(0);
    }
    bounds
}

/// For every `i`, whether `s[i..i + |t|]` matches `t`, where `*` in either string matches any
/// character (both over `a..z` and `*`). Uses `sum (s - t)^2 s t = s^3 t - 2 s^2 t^2 + s t^3`
/// over the non-wildcard positions (letters are `1..=26`, wildcard `0`), three exact
/// convolutions.
pub fn wildcard_match(s: &[u8], t: &[u8]) -> Vec<bool> {
    assert!(t.len() <= s.len() && !t.is_empty(), "need 1 <= |t| <= |s|");
    let value = |c: u8| -> u64 {
        if c == b'*' {
            0
        } else {
            assert!(c.is_ascii_lowercase(), "only lowercase letters and '*'");
            (c - b'a' + 1) as u64
        }
    };
    let a: Vec<u64> = s.iter().map(|&c| value(c)).collect();
    let b: Vec<u64> = t.iter().rev().map(|&c| value(c)).collect();
    let pow = |v: &[u64], e: u32| -> Vec<u64> { v.iter().map(|&x| x.pow(e)).collect() };
    let c1 = convolution_u64(&pow(&a, 3), &b);
    let c2 = convolution_u64(&pow(&a, 2), &pow(&b, 2));
    let c3 = convolution_u64(&a, &pow(&b, 3));
    let m = t.len();
    (0..=s.len() - m)
        .map(|i| {
            let k = i + m - 1;
            c1[k] + c3[k] == 2 * c2[k]
        })
        .collect()
}
