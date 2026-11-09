//! Stress test: `SuffixAutomaton` against `str::contains`, a substring set, structural
//! invariants, and a suffix-array/LCP count of distinct substrings at scale.

use std::collections::HashSet;
use suffix_automaton::SuffixAutomaton;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

const ALPHABET: [char; 5] = ['a', 'b', 'c', 'é', '日'];

fn text(rng: &mut Rng, len: usize, letters: usize) -> String {
    (0..len).map(|_| ALPHABET[rng.below(letters)]).collect()
}

fn check_structure(sam: &SuffixAutomaton, n: usize, ctx: &str) {
    assert_eq!(sam.states[sam.last].len, n, "last state has the full length {ctx}");
    assert!(sam.states.len() <= 2 * n.max(2), "at most 2n states {ctx}");
    assert_eq!(sam.states[0].link, None);
    for (v, s) in sam.states.iter().enumerate().skip(1) {
        let l = s.link.expect("every non-initial state has a suffix link");
        assert!(sam.states[l].len < s.len, "link lengths strictly decrease at {v} {ctx}");
        for (_, &to) in &s.next {
            assert!(sam.states[to].len > s.len || sam.states[to].len >= 1, "sane transition at {v} {ctx}");
        }
    }
}

#[test]
fn matches_naive_substring_queries() {
    let mut st = Stress::new("suffix_automaton/naive");
    for case in 0..800 * scale() {
        let mut rng = Rng::new(case as u64);
        let letters = rng.range(1, 5) as usize;
        let len = rng.below(35);
        let s = text(&mut rng, len, letters);
        let sam = SuffixAutomaton::from_str(&s);
        let chars: Vec<char> = s.chars().collect();
        st.case();
        check_structure(&sam, chars.len(), &format!("case {case} s={s:?}"));

        let mut subs: HashSet<Vec<char>> = HashSet::new();
        for i in 0..chars.len() {
            for j in i + 1..=chars.len() {
                subs.insert(chars[i..j].to_vec());
            }
        }
        assert_eq!(sam.distinct_substring_count(), subs.len() as u64, "distinct count s={s:?}");
        for sub in &subs {
            assert!(sam.contains(&sub.iter().collect::<String>()), "substring of {s:?}");
        }
        assert!(sam.contains(""), "empty pattern");
        for _ in 0..25 {
            let plen = rng.below(8);
            let p = text(&mut rng, plen, (letters + case % 2).min(ALPHABET.len())); // sometimes uses an unseen letter
            assert_eq!(sam.contains(&p), s.contains(&p), "contains({p:?}) in {s:?}");
        }
        // Incremental construction equals from_str.
        let mut inc = SuffixAutomaton::new();
        for &c in &chars {
            inc.extend(c);
        }
        assert_eq!(inc.distinct_substring_count(), sam.distinct_substring_count());
        st.ops(subs.len() + 27);
    }
    st.done();
}

/// Distinct substrings = n(n+1)/2 - sum(LCP) over a suffix array built by prefix doubling.
fn distinct_via_suffix_array(s: &[u8]) -> u64 {
    let n = s.len();
    let mut sa: Vec<usize> = (0..n).collect();
    let mut rank: Vec<i64> = s.iter().map(|&c| c as i64).collect();
    let mut k = 1;
    let mut tmp = vec![0i64; n];
    while k < n {
        let key = |i: usize, rank: &Vec<i64>| (rank[i], if i + k < n { rank[i + k] } else { -1 });
        sa.sort_by_key(|&i| key(i, &rank));
        tmp[sa[0]] = 0;
        for w in 1..n {
            tmp[sa[w]] = tmp[sa[w - 1]] + (key(sa[w], &rank) != key(sa[w - 1], &rank)) as i64;
        }
        rank.copy_from_slice(&tmp);
        if rank[sa[n - 1]] as usize == n - 1 {
            break;
        }
        k *= 2;
    }
    // Kasai's LCP
    let mut inv = vec![0usize; n];
    for (i, &p) in sa.iter().enumerate() {
        inv[p] = i;
    }
    let (mut h, mut lcp_sum) = (0usize, 0u64);
    for i in 0..n {
        if inv[i] > 0 {
            let j = sa[inv[i] - 1];
            while i + h < n && j + h < n && s[i + h] == s[j + h] {
                h += 1;
            }
            lcp_sum += h as u64;
            h = h.saturating_sub(1);
        } else {
            h = 0;
        }
    }
    (n as u64 * (n as u64 + 1)) / 2 - lcp_sum
}

#[test]
fn large_distinct_count_matches_suffix_array() {
    let mut st = Stress::new("suffix_automaton/large");
    let mut rng = Rng::new(23);
    for (n, letters) in [(20_000usize, 2usize), (20_000, 3), (20_000, 26), (50_000, 2)] {
        let s: String = (0..n).map(|_| (b'a' + rng.below(letters) as u8) as char).collect();
        let sam = SuffixAutomaton::from_str(&s);
        st.case();
        assert_eq!(sam.distinct_substring_count(), distinct_via_suffix_array(s.as_bytes()), "n={n} letters={letters}");
        assert!(sam.states.len() <= 2 * n);
        for _ in 0..2_000 {
            let i = rng.below(n);
            let j = (i + 1 + rng.below(40)).min(n);
            assert!(sam.contains(&s[i..j]), "a substring of s must be found");
            st.ops(1);
        }
        st.ops(n);
    }
    // periodic and constant strings have closed forms
    let a = SuffixAutomaton::from_str(&"a".repeat(10_000));
    assert_eq!(a.distinct_substring_count(), 10_000);
    st.done();
}
