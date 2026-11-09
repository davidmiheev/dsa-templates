//! Stress test: `BinaryTrie` against a plain multiset (sorted `Vec`).

use binary_trie::BinaryTrie;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn run<const BITS: usize>(st: &mut Stress, case: usize, universe: u64, steps: usize) {
    let mut rng = Rng::new(case as u64 * 31 + BITS as u64);
    let mut trie = BinaryTrie::<BITS>::new();
    let mut model: Vec<u64> = vec![];
    st.case();
    for step in 0..steps {
        let x = rng.next_u64() % universe;
        let ctx = format!("BITS={BITS} case {case} step {step} x={x}");
        match rng.below(5) {
            0 | 1 => {
                trie.insert(x);
                model.push(x);
            }
            2 => {
                let at = model.iter().position(|&k| k == x);
                assert_eq!(trie.remove(x), at.is_some(), "remove {ctx}");
                if let Some(i) = at {
                    model.swap_remove(i);
                }
            }
            3 => {
                assert_eq!(trie.count(x), model.iter().filter(|&&k| k == x).count(), "count {ctx}");
                assert_eq!(trie.contains(x), model.contains(&x), "contains {ctx}");
            }
            _ => {
                if !model.is_empty() {
                    assert_eq!(trie.min_xor(x), model.iter().map(|&k| k ^ x).min().unwrap(), "min_xor {ctx}");
                    assert_eq!(trie.max_xor(x), model.iter().map(|&k| k ^ x).max().unwrap(), "max_xor {ctx}");
                }
            }
        }
        assert_eq!((trie.len(), trie.is_empty()), (model.len(), model.is_empty()), "len {ctx}");
        st.ops(1);
    }
}

#[test]
fn matches_multiset() {
    let mut st = Stress::new("binary_trie/naive");
    for case in 0..400 * scale() {
        run::<4>(&mut st, case, 1 << 4, 80); // dense: duplicates and removals collide a lot
        run::<30>(&mut st, case, 1 << 30, 80);
        run::<64>(&mut st, case, u64::MAX, 80);
    }
    st.done();
}

#[test]
fn keys_must_fit_and_queries_need_keys() {
    let mut t = BinaryTrie::<8>::new();
    assert!(catch_unwind(AssertUnwindSafe(|| t.min_xor(0))).is_err()); // empty
    assert!(catch_unwind(AssertUnwindSafe(|| t.insert(256))).is_err()); // does not fit
    assert_eq!((t.count(300), t.contains(1000)), (0, false)); // out-of-range lookups are just absent
    t.insert(255);
    assert_eq!((t.min_xor(255), t.max_xor(0), t.max_xor(255)), (0, 255, 0));
}

#[test]
fn large_matches_sorted_model() {
    // 300k operations over 30-bit keys; min/max xor verified against a linear scan on a sample.
    let mut st = Stress::new("binary_trie/large");
    let mut rng = Rng::new(5);
    let mut trie = BinaryTrie::<30>::new();
    let mut keys: Vec<u64> = vec![];
    st.case();
    for i in 0..300_000 * scale() {
        let x = rng.next_u64() % (1 << 30);
        if rng.chance(2, 3) || keys.is_empty() {
            trie.insert(x);
            keys.push(x);
        } else {
            let j = rng.below(keys.len());
            let k = keys.swap_remove(j);
            assert!(trie.remove(k));
        }
        if i % 1000 == 0 && !keys.is_empty() {
            let q = rng.next_u64() % (1 << 30);
            assert_eq!(trie.min_xor(q), keys.iter().map(|&k| k ^ q).min().unwrap());
            assert_eq!(trie.max_xor(q), keys.iter().map(|&k| k ^ q).max().unwrap());
        }
        st.ops(1);
    }
    assert_eq!(trie.len(), keys.len());
    st.done();
}
