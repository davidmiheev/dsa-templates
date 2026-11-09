//! Stress test: `Bst` against `BTreeSet`.

use bst::Bst;
use std::collections::BTreeSet;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_btreeset() {
    let mut st = Stress::new("bst/naive");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let span = *[4i64, 20, 100, 1_000_000_000_000].get(case % 4).unwrap();
        let mut bst = Bst::new();
        let mut set = BTreeSet::new();
        st.case();
        assert_eq!((bst.min(), bst.max(), bst.inorder()), (None, None, vec![]), "empty tree");
        for step in 0..120 {
            let key = rng.range(-span, span);
            if rng.chance(2, 3) {
                assert_eq!(bst.insert(key), set.insert(key), "insert({key}) case {case} step {step}");
            } else {
                assert_eq!(bst.contains(key), set.contains(&key), "contains({key}) case {case} step {step}");
            }
            assert_eq!(bst.min(), set.iter().next().copied(), "min case {case} step {step}");
            assert_eq!(bst.max(), set.iter().next_back().copied(), "max case {case} step {step}");
            st.ops(1);
        }
        assert_eq!(bst.inorder(), set.iter().copied().collect::<Vec<_>>(), "inorder case {case}");
        st.ops(1);
    }
    st.done();
}

#[test]
fn degenerate_sorted_insertions() {
    // Sorted input makes this non-balancing tree a linked list; recursion in the
    // crate must still cope with a moderately deep tree.
    let mut st = Stress::new("bst/sorted");
    let n = 3_000;
    let mut bst = Bst::new();
    st.case();
    for k in 0..n {
        assert!(bst.insert(k));
    }
    for k in 0..n {
        assert!(!bst.insert(k), "duplicate insert must report false");
        assert!(bst.contains(k));
    }
    assert!(!bst.contains(n) && !bst.contains(-1));
    assert_eq!(bst.inorder(), (0..n).collect::<Vec<_>>());
    st.ops(3 * n as usize);
    st.done();
}
