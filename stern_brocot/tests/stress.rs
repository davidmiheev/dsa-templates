//! Stress test against a literal walk of the tree and brute-force best approximations.

use stern_brocot::*;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// Walks `moves` from the root, returning every node on the way (including the root).
fn walk(moves: &[bool]) -> Vec<(u64, u64)> {
    let (mut lo, mut hi) = ((0u64, 1u64), (1u64, 0u64));
    let mut nodes = vec![(1, 1)];
    for &right in moves {
        let mid = (lo.0 + hi.0, lo.1 + hi.1);
        if right {
            lo = mid;
        } else {
            hi = mid;
        }
        nodes.push((lo.0 + hi.0, lo.1 + hi.1));
    }
    nodes
}

fn expand(runs: &[Run]) -> Vec<bool> {
    runs.iter().flat_map(|&(c, n)| std::iter::repeat_n(c, n as usize)).collect()
}

#[test]
fn tree_operations_match_walks() {
    let mut st = Stress::new("stern_brocot/tree");
    for case in 0..20_000 * scale() {
        let mut rng = Rng::new(case as u64);
        let len = rng.below(25);
        let moves: Vec<bool> = (0..len).map(|_| rng.chance(1, 2)).collect();
        let nodes = walk(&moves);
        let u = *nodes.last().unwrap();
        let runs = encode_path(u.0, u.1);
        assert_eq!(expand(&runs), moves);
        assert!(runs.windows(2).all(|w| w[0].0 != w[1].0));
        assert_eq!(decode_path(&runs), u);
        // ancestors
        for (k, &node) in nodes.iter().enumerate() {
            assert_eq!(ancestor(k as u64, u), node);
        }
        assert_eq!(ancestor(len as u64 + 1, u), (0, 0));
        // lca with another node
        let len2 = rng.below(25);
        let moves2: Vec<bool> = (0..len2).map(|_| rng.chance(1, 2)).collect();
        let nodes2 = walk(&moves2);
        let v = *nodes2.last().unwrap();
        let common = moves.iter().zip(&moves2).take_while(|(a, b)| a == b).count();
        assert_eq!(lca(u, v), nodes[common]);
        // subtree range: the descendants are exactly the fractions strictly between the ends
        let (lo, hi) = descendant_range(u);
        assert!(lo.0 * u.1 < u.0 * lo.1, "lo < u");
        assert!(u.0 * hi.1 < hi.0 * u.1, "u < hi");
        assert_eq!(hi.0 * lo.1 - lo.0 * hi.1, 1, "neighbours in the Farey sense");
        st.case();
        st.ops(len);
    }
    st.done();
}

#[test]
fn approximations_match_brute_force() {
    let mut st = Stress::new("stern_brocot/approx");
    for case in 0..20_000 * scale() {
        let mut rng = Rng::new(10_000 + case as u64);
        let n = rng.range(1, 25) as u64;
        let (x, y) = (rng.range(1, 80) as u64, rng.range(1, 80) as u64);
        let mut lower = (0u64, 1u64);
        let mut upper = (1u64, 0u64);
        for p in 1..=n {
            for q in 1..=n {
                // p/q <= x/y
                if p * y <= x * q && p * lower.1 > lower.0 * q {
                    lower = (p, q);
                }
                if p * y >= x * q && (upper.1 == 0 || p * upper.1 < upper.0 * q) {
                    upper = (p, q);
                }
            }
        }
        let norm = |f: (u64, u64)| {
            let g = gcd(f.0, f.1);
            (f.0 / g, f.1 / g)
        };
        let (lo, hi) = best_approximations(n, x, y);
        assert_eq!((norm(lo), norm(hi)), (norm(lower), norm(upper)), "n={n} x={x} y={y}");
        st.case();
    }
    // large values
    let (lo, hi) = best_approximations(1_000_000_000, 1_000_000_000, 999_999_999);
    assert_eq!((lo, hi), ((1_000_000_000, 999_999_999), (1_000_000_000, 999_999_999)));
    st.done();
}
