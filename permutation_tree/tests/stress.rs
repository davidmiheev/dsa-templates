//! Stress test: `PermutationTree` against the definition (enumerate common intervals, keep the
//! strong ones, classify linear / prime by checking runs of consecutive children).

use permutation_tree::{Kind, PermutationTree};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn random_perm(rng: &mut Rng, n: usize) -> Vec<usize> {
    let mut p: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        p.swap(i, rng.below(i + 1));
    }
    p
}

fn is_common(p: &[usize], l: usize, r: usize) -> bool {
    let (mn, mx) = p[l..=r].iter().fold((usize::MAX, 0), |(a, b), &x| (a.min(x), b.max(x)));
    mx - mn == r - l
}

/// (l, r, linear?) of every strong interval, parent index in the sorted list.
fn brute(p: &[usize]) -> Vec<(usize, usize, bool)> {
    let n = p.len();
    let common: Vec<(usize, usize)> = (0..n).flat_map(|l| (l..n).map(move |r| (l, r))).filter(|&(l, r)| is_common(p, l, r)).collect();
    let overlaps = |a: (usize, usize), b: (usize, usize)| (a.0 < b.0 && b.0 <= a.1 && a.1 < b.1) || (b.0 < a.0 && a.0 <= b.1 && b.1 < a.1);
    let strong: Vec<(usize, usize)> = common.iter().copied().filter(|&a| !common.iter().any(|&b| overlaps(a, b))).collect();
    // children of an interval: maximal strong intervals strictly inside, in order
    strong
        .iter()
        .map(|&(l, r)| {
            let mut kids: Vec<(usize, usize)> = strong.iter().copied().filter(|&(a, b)| (a, b) != (l, r) && l <= a && b <= r).collect();
            kids.retain(|&(a, b)| !strong.iter().any(|&(c, d)| (c, d) != (a, b) && (c, d) != (l, r) && l <= c && d <= r && c <= a && b <= d));
            kids.sort_unstable();
            // linear iff every run of >= 2 consecutive children (and < all) is a common interval
            let linear = if kids.is_empty() {
                true
            } else {
                (0..kids.len()).all(|i| (i + 1..kids.len()).all(|j| is_common(p, kids[i].0, kids[j].1)))
            };
            (l, r, linear)
        })
        .collect()
}

#[test]
fn matches_definition_on_small_permutations() {
    let mut st = Stress::new("permutation_tree/naive");
    for case in 0..3000 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, if case % 4 == 0 { 13 } else { 9 }) as usize;
        let p = if case % 3 == 0 {
            // Make structured permutations (blocks), which exercise linear nodes.
            let mut q: Vec<usize> = Vec::new();
            let mut base = 0;
            while base < n {
                let len = 1 + rng.below(n - base);
                let mut block: Vec<usize> = (base..base + len).collect();
                if rng.chance(1, 2) {
                    block.reverse();
                }
                q.extend(block);
                base += len;
            }
            q
        } else {
            random_perm(&mut rng, n)
        };
        let t = PermutationTree::new(&p);
        let mut got: Vec<(usize, usize, bool)> = (0..t.len()).map(|v| (t.left[v], t.right[v], t.is_linear(v))).collect();
        got.sort_unstable();
        let mut want = brute(&p);
        want.sort_unstable();
        assert_eq!(got, want, "p={p:?}");
        // Structure: the root spans everything, children tile the parent, leaves have no children.
        assert_eq!((t.left[t.root], t.right[t.root]), (0, n - 1), "p={p:?}");
        for v in 0..t.len() {
            let kids = t.children(v);
            if t.kind[v] == Kind::Leaf {
                assert!(kids.is_empty());
                assert_eq!(t.left[v], t.right[v]);
            } else {
                assert!(kids.len() >= 2, "inner node with < 2 children, p={p:?}");
                assert_eq!(t.left[kids[0]], t.left[v]);
                assert_eq!(t.right[*kids.last().unwrap()], t.right[v]);
                for w in kids.windows(2) {
                    assert_eq!(t.right[w[0]] + 1, t.left[w[1]], "children tile, p={p:?}");
                }
                for &c in kids {
                    assert_eq!(t.parent[c], v);
                }
            }
        }
        st.case();
        st.ops(1);
    }
    st.done();
}

#[test]
fn large_permutations() {
    let mut st = Stress::new("permutation_tree/large");
    let n = 500_000;
    let mut rng = Rng::new(6);
    st.case();
    for p in [random_perm(&mut rng, n), (0..n).collect(), (0..n).rev().collect(), {
        // blocks of mixed directions
        let mut q = Vec::with_capacity(n);
        let mut base = 0;
        while base < n {
            let len = (1 + rng.below(50)).min(n - base);
            let mut b: Vec<usize> = (base..base + len).collect();
            if rng.chance(1, 2) {
                b.reverse();
            }
            q.extend(b);
            base += len;
        }
        q
    }] {
        let t = PermutationTree::new(&p);
        assert_eq!((t.left[t.root], t.right[t.root]), (0, n - 1));
        assert!(t.len() < 2 * n);
    }
    st.ops(4 * n);
    st.done();
}

#[test]
#[should_panic(expected = "not a permutation")]
fn rejects_non_permutations() {
    PermutationTree::new(&[0, 0, 1]);
}
