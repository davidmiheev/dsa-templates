//! Stress test: existence against the brute-force definition, and validity of every output.

use st_numbering::st_numbering;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn valid(n: usize, edges: &[(usize, usize)], s: usize, t: usize, p: &[usize]) -> bool {
    let mut seen = vec![false; n];
    for &x in p {
        if x >= n || seen[x] {
            return false;
        }
        seen[x] = true;
    }
    if p[s] != 0 || p[t] != n - 1 {
        return false;
    }
    // every vertex is reachable from s and reaches t along increasing numbers
    let mut lower = vec![false; n];
    let mut higher = vec![false; n];
    for &(u, v) in edges {
        let (a, b) = if p[u] < p[v] { (u, v) } else { (v, u) };
        higher[a] = true;
        lower[b] = true;
    }
    (0..n).all(|v| v == s || v == t || (lower[v] && higher[v]))
}

/// Brute force: try every permutation with p[s] = 0, p[t] = n - 1.
fn brute(n: usize, edges: &[(usize, usize)], s: usize, t: usize) -> bool {
    let mid: Vec<usize> = (0..n).filter(|&v| v != s && v != t).collect();
    let mut perm: Vec<usize> = (1..n - 1).collect();
    loop {
        let mut p = vec![0usize; n];
        p[t] = n - 1;
        for (i, &v) in mid.iter().enumerate() {
            p[v] = perm[i];
        }
        if valid(n, edges, s, t, &p) {
            return true;
        }
        // next permutation
        let mut i = perm.len();
        while i >= 2 && perm[i - 2] >= perm[i - 1] {
            i -= 1;
        }
        if i < 2 {
            return false;
        }
        let mut j = perm.len() - 1;
        while perm[j] <= perm[i - 2] {
            j -= 1;
        }
        perm.swap(i - 2, j);
        perm[i - 1..].reverse();
    }
}

#[test]
fn existence_and_validity() {
    let mut st = Stress::new("st_numbering");
    let (mut yes, mut no) = (0, 0);
    for case in 0..5000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(2, 7) as usize;
        let m = rng.range(0, 12) as usize;
        let edges: Vec<(usize, usize)> = (0..m)
            .map(|_| (rng.below(n), rng.below(n)))
            .filter(|&(a, b)| a != b)
            .collect();
        let s = rng.below(n);
        let t = (s + 1 + rng.below(n - 1)) % n;
        let p = st_numbering(n, &edges, s, t);
        // For n == 2 the definition needs an actual s-t edge.
        let expect = if n == 2 { !edges.is_empty() } else { brute(n, &edges, s, t) };
        if p.is_empty() {
            no += 1;
            assert!(!expect, "case {case}: {edges:?} s={s} t={t} should have a numbering");
        } else {
            yes += 1;
            assert!(expect, "case {case}: {edges:?} s={s} t={t} must be rejected");
            assert!(valid(n, &edges, s, t, &p), "case {case}: {edges:?} s={s} t={t} p={p:?}");
        }
        st.case();
    }
    assert!(yes > 100 && no > 100, "{yes}/{no}");
    st.done();
}
