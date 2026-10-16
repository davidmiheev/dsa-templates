//! Stress tests of every string algorithm against brute force.
#![allow(clippy::needless_range_loop)]

use std::collections::{BTreeSet, HashMap, HashSet};
use strings::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn random_string(rng: &mut Rng, len: usize, sigma: usize) -> Vec<u8> {
    (0..len).map(|_| b'a' + rng.below(sigma) as u8).collect()
}

fn is_pal(s: &[u8]) -> bool {
    s.iter().eq(s.iter().rev())
}

fn palindromic_substrings(s: &[u8]) -> BTreeSet<Vec<u8>> {
    let mut set = BTreeSet::new();
    for i in 0..s.len() {
        for j in i + 1..=s.len() {
            if is_pal(&s[i..j]) {
                set.insert(s[i..j].to_vec());
            }
        }
    }
    set
}

#[test]
fn suffix_structures() {
    let mut st = Stress::new("strings/suffix");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = 1 + rng.below(if case % 20 == 0 { 120 } else { 14 });
        let sigma = 1 + rng.below(4);
        let s = random_string(&mut rng, n, sigma);
        let mut naive: Vec<usize> = (0..n).collect();
        naive.sort_by(|&a, &b| s[a..].cmp(&s[b..]));
        let sa = suffix_array_bytes(&s);
        assert_eq!(sa, naive, "{:?}", String::from_utf8_lossy(&s));
        let lcp = lcp_array(&s, &sa);
        for i in 0..n - 1 {
            let (a, b) = (&s[sa[i]..], &s[sa[i + 1]..]);
            assert_eq!(lcp[i], a.iter().zip(b).take_while(|(x, y)| x == y).count());
        }
        let distinct: HashSet<&[u8]> = (0..n).flat_map(|i| (i + 1..=n).map(move |j| (i, j))).map(|(i, j)| &s[i..j]).collect();
        assert_eq!(count_distinct_substrings(&s) as usize, distinct.len());
        // longest common substring
        let tl = 1 + rng.below(12);
        let t = random_string(&mut rng, tl, sigma);
        let (a, b, c, d) = longest_common_substring(&s, &t);
        let best = (0..n).flat_map(|i| (0..t.len()).map(move |j| (i, j))).map(|(i, j)| s[i..].iter().zip(&t[j..]).take_while(|(x, y)| x == y).count()).max().unwrap();
        assert_eq!(b - a, best);
        assert_eq!(d - c, best);
        assert_eq!(s[a..b], t[c..d]);
        // z, manacher, lyndon
        let z = z_algorithm(&s);
        for i in 0..n {
            assert_eq!(z[i], s.iter().zip(&s[i..]).take_while(|(x, y)| x == y).count());
        }
        let man = manacher(&s);
        assert_eq!(man.len(), 2 * n - 1);
        for (c, &len) in man.iter().enumerate() {
            let (lo, hi) = if c % 2 == 0 { (c / 2, c / 2) } else { (c / 2 + 1, c / 2 + 1) };
            let mut best = 0;
            let (mut l, mut r) = (lo as i64 - if c % 2 == 0 { 0 } else { 1 }, hi as i64);
            if c % 2 == 0 {
                best = 1;
                l -= 1;
                r += 1;
            }
            while l >= 0 && (r as usize) < n && s[l as usize] == s[r as usize] {
                best += 2;
                l -= 1;
                r += 1;
            }
            assert_eq!(len, best, "centre {c} of {:?}", String::from_utf8_lossy(&s));
        }
        let bounds = lyndon_factorization(&s);
        assert_eq!((bounds[0], *bounds.last().unwrap()), (0, n));
        let factors: Vec<&[u8]> = bounds.windows(2).map(|w| &s[w[0]..w[1]]).collect();
        for f in &factors {
            // Lyndon: strictly smaller than all proper suffixes
            assert!((1..f.len()).all(|k| f[..] < f[k..]), "factor {f:?}");
        }
        assert!(factors.windows(2).all(|w| w[0] >= w[1]));
        st.case();
        st.ops(n);
    }
    // bigger inputs agree with the doubling-free naive on a repetitive string
    let s: Vec<u8> = std::iter::repeat_n(b"abaab".iter().copied(), 300).flatten().collect();
    let sa = suffix_array_bytes(&s);
    assert!(sa.windows(2).all(|w| s[w[0]..] < s[w[1]..]));
    st.done();
}

#[test]
fn aho_corasick_and_eertree() {
    let mut st = Stress::new("strings/trie");
    for case in 0..2000 * scale() {
        let mut rng = Rng::new(10_000 + case as u64);
        let sigma = 1 + rng.below(3);
        let count = 1 + rng.below(5);
        let pats: Vec<Vec<u8>> = (0..count).map(|_| (0..1 + rng.below(6)).map(|_| rng.below(sigma) as u8).collect()).collect();
        let refs: Vec<&[u8]> = pats.iter().map(|p| p.as_slice()).collect();
        let ac = AhoCorasick::new(&refs, sigma);
        // the strings of all prefixes in creation order
        let mut names: Vec<Vec<u8>> = vec![Vec::new()];
        let mut ids: HashMap<Vec<u8>, usize> = HashMap::from([(Vec::new(), 0)]);
        for p in &pats {
            for j in 1..=p.len() {
                if !ids.contains_key(&p[..j]) {
                    ids.insert(p[..j].to_vec(), names.len());
                    names.push(p[..j].to_vec());
                }
            }
        }
        assert_eq!(ac.len(), names.len());
        for v in 1..names.len() {
            assert_eq!(ac.parent(v), ids[&names[v][..names[v].len() - 1]]);
            let link = (1..names[v].len()).find_map(|k| ids.get(&names[v][k..]).copied()).unwrap_or(0);
            assert_eq!(ac.suffix_link(v), link, "link of {:?}", names[v]);
        }
        for (i, p) in pats.iter().enumerate() {
            assert_eq!(ac.pattern_vertex(i), ids[p]);
        }
        // eertree against the set of palindromic substrings
        let n = 1 + rng.below(16);
        let s: Vec<u8> = (0..n).map(|_| rng.below(sigma.max(2)) as u8).collect();
        let mut t = Eertree::new(sigma.max(2));
        for &c in &s {
            t.push(c);
        }
        let pals = palindromic_substrings(&s.iter().map(|&c| b'a' + c).collect::<Vec<_>>());
        assert_eq!(t.distinct_palindromes(), pals.len());
        // vertex k + 1 is the k-th palindrome created: it is the first one completed in the string
        let mut order: Vec<Vec<u8>> = Vec::new();
        let mut seen = HashSet::new();
        for j in 1..=n {
            for i in 0..j {
                let sub = &s[i..j];
                if is_pal(sub) && !seen.contains(sub) {
                    seen.insert(sub.to_vec());
                    order.push(sub.to_vec());
                }
            }
        }
        // within one end position only the longest new suffix palindrome is new; sort by (end, length)
        let mut by_end: Vec<(usize, Vec<u8>)> = Vec::new();
        let mut seen2 = HashSet::new();
        for j in 1..=n {
            let mut news: Vec<Vec<u8>> = (0..j).map(|i| s[i..j].to_vec()).filter(|x| is_pal(x) && !seen2.contains(x)).collect();
            news.sort_by_key(|x| x.len());
            for x in news {
                seen2.insert(x.clone());
                by_end.push((j, x));
            }
        }
        assert_eq!(by_end.len(), t.distinct_palindromes());
        let vid: HashMap<Vec<u8>, usize> = by_end.iter().enumerate().map(|(k, (_, x))| (x.clone(), k + 2)).collect();
        for (k, (_, x)) in by_end.iter().enumerate() {
            let v = k + 2;
            assert_eq!(t.length(v), x.len());
            let par = if x.len() <= 2 { (x.len() % 2) ^ 1 } else { vid[&x[1..x.len() - 1]] };
            assert_eq!(t.parent(v), par, "parent of {x:?}");
            let link = (1..x.len()).find(|&i| is_pal(&x[i..])).map_or(1, |i| vid[&x[i..]]);
            assert_eq!(t.suffix_link(v), link, "link of {x:?}");
        }
        for j in 0..n {
            let longest = (0..=j).find(|&i| is_pal(&s[i..=j])).unwrap();
            assert_eq!(t.length(t.longest_suffix_palindrome(j)), j + 1 - longest);
        }
        st.case();
    }
    st.done();
}

#[test]
fn deque_eertree_matches_brute_force() {
    let mut st = Stress::new("strings/deque_eertree");
    for case in 0..6000 * scale() {
        let mut rng = Rng::new(20_000 + case as u64);
        let sigma = 1 + rng.below(3);
        let mut t = DequeEertree::new(sigma);
        let mut s: std::collections::VecDeque<u8> = Default::default();
        let ops = 1 + rng.below(if case % 10 == 0 { 120 } else { 30 });
        for _ in 0..ops {
            let kind = if s.is_empty() { rng.below(2) } else { rng.below(4) };
            match kind {
                0 => {
                    let c = rng.below(sigma) as u8;
                    t.push_front(c);
                    s.push_front(c);
                }
                1 => {
                    let c = rng.below(sigma) as u8;
                    t.push_back(c);
                    s.push_back(c);
                }
                2 => {
                    t.pop_front();
                    s.pop_front();
                }
                _ => {
                    t.pop_back();
                    s.pop_back();
                }
            }
            let v: Vec<u8> = s.iter().map(|&c| b'a' + c).collect();
            assert_eq!(t.distinct_palindromes(), palindromic_substrings(&v).len(), "string {:?}", String::from_utf8_lossy(&v));
            let pre = (0..=v.len()).rev().find(|&k| k == 0 || is_pal(&v[..k])).unwrap();
            let suf = (0..=v.len()).rev().find(|&k| k == 0 || is_pal(&v[v.len() - k..])).unwrap();
            assert_eq!(t.longest_prefix_palindrome(), pre, "prefix of {:?}", String::from_utf8_lossy(&v));
            assert_eq!(t.longest_suffix_palindrome(), suf, "suffix of {:?}", String::from_utf8_lossy(&v));
            assert_eq!(t.len(), v.len());
        }
        st.case();
        st.ops(ops);
    }
    st.done();
}

#[test]
fn runs_and_wildcards() {
    let mut st = Stress::new("strings/runs");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(30_000 + case as u64);
        let n = 1 + rng.below(if case % 10 == 0 { 40 } else { 14 });
        let sg = 1 + rng.below(3);
        let s = random_string(&mut rng, n, sg);
        // brute force: every (t, l, r) with minimal period t, r - l >= 2t, maximal
        let mut want = Vec::new();
        for l in 0..n {
            for r in l + 2..=n {
                let sub = &s[l..r];
                let period = (1..=sub.len()).find(|&p| (p..sub.len()).all(|i| sub[i] == sub[i - p])).unwrap();
                if r - l < 2 * period {
                    continue;
                }
                // maximal: not extendable with the same period
                let ext_l = l > 0 && s[l - 1] == s[l - 1 + period];
                let ext_r = r < n && s[r] == s[r - period];
                if !ext_l && !ext_r {
                    want.push((period, l, r));
                }
            }
        }
        want.sort_unstable();
        assert_eq!(runs(&s), want, "{:?}", String::from_utf8_lossy(&s));
        // wildcard matching
        let m = 1 + rng.below(n);
        let star = |rng: &mut Rng, len: usize| -> Vec<u8> { (0..len).map(|_| if rng.chance(1, 4) { b'*' } else { b'a' + rng.below(3) as u8 }).collect() };
        let a = star(&mut rng, n);
        let b = star(&mut rng, m);
        let got = wildcard_match(&a, &b);
        for i in 0..=n - m {
            let ok = (0..m).all(|k| a[i + k] == b'*' || b[k] == b'*' || a[i + k] == b[k]);
            assert_eq!(got[i], ok);
        }
        st.case();
    }
    st.done();
}
