//! Stress test: treap-based `SortedList` against a sorted `Vec`.

use sorted_list::SortedList;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn check_all(list: &SortedList, model: &[i64], ctx: &str) {
    assert_eq!(list.len(), model.len(), "len {ctx}");
    assert_eq!(list.is_empty(), model.is_empty(), "is_empty {ctx}");
    for (i, &v) in model.iter().enumerate() {
        assert_eq!(list.get(i), Some(v), "get({i}) {ctx}");
    }
    assert_eq!(list.get(model.len()), None, "get(len) {ctx}");
    assert_eq!(list.get(model.len() + 5), None, "get(len+5) {ctx}");
}

#[test]
fn matches_sorted_vec() {
    let mut st = Stress::new("sorted_list/naive");
    for case in 0..300 * scale() {
        let mut rng = Rng::new(case as u64);
        let span = *[3i64, 8, 30, 1_000_000_000].get(case % 4).unwrap();
        let mut list = SortedList::new();
        let mut model: Vec<i64> = Vec::new();
        st.case();
        for step in 0..150 {
            let key = rng.range(-span, span);
            let ctx = format!("case {case} step {step} key {key}");
            match rng.below(5) {
                0 | 1 => {
                    list.insert(key);
                    let at = model.partition_point(|&x| x < key);
                    model.insert(at, key);
                }
                2 => {
                    list.remove(key); // removing an absent key must be a no-op
                    if let Ok(at) = model.binary_search(&key) {
                        // remove the first occurrence (all equal keys are interchangeable)
                        let first = model.partition_point(|&x| x < key);
                        assert_eq!(model[at], key);
                        model.remove(first);
                    }
                }
                3 => {
                    assert_eq!(list.bisect_left(key), model.partition_point(|&x| x < key), "bisect_left {ctx}");
                    assert_eq!(list.bisect_right(key), model.partition_point(|&x| x <= key), "bisect_right {ctx}");
                }
                _ => {
                    let i = rng.below(model.len() + 2);
                    assert_eq!(list.get(i), model.get(i).copied(), "get({i}) {ctx}");
                }
            }
            assert_eq!(list.len(), model.len(), "len {ctx}");
            st.ops(1);
        }
        check_all(&list, &model, &format!("final case {case}"));
        st.ops(model.len() + 3);
    }
    st.done();
}

#[test]
fn drain_and_refill_reuses_nodes() {
    let mut st = Stress::new("sorted_list/drain");
    let mut rng = Rng::new(42);
    let mut list = SortedList::new();
    let mut model: Vec<i64> = Vec::new();
    st.case();
    for round in 0..20 * scale() {
        for _ in 0..500 {
            let k = rng.range(-50, 50);
            list.insert(k);
            let at = model.partition_point(|&x| x < k);
            model.insert(at, k);
        }
        while !model.is_empty() {
            let k = model[rng.below(model.len())];
            list.remove(k);
            let first = model.partition_point(|&x| x < k);
            model.remove(first);
            st.ops(1);
        }
        check_all(&list, &model, &format!("drained round {round}"));
        assert!(list.is_empty());
    }
    st.done();
}

#[test]
fn large_random_matches_sorted_vec() {
    let mut st = Stress::new("sorted_list/large");
    let mut rng = Rng::new(11);
    let mut list = SortedList::new();
    let mut model: Vec<i64> = Vec::new();
    st.case();
    for _ in 0..40_000 * scale() {
        let key = rng.range(-1_000_000, 1_000_000);
        match rng.below(4) {
            0 | 1 => {
                list.insert(key);
                let at = model.partition_point(|&x| x < key);
                model.insert(at, key);
            }
            2 => {
                list.remove(key);
                let first = model.partition_point(|&x| x < key);
                if model.get(first) == Some(&key) {
                    model.remove(first);
                }
            }
            _ => {
                assert_eq!(list.bisect_left(key), model.partition_point(|&x| x < key));
                assert_eq!(list.bisect_right(key), model.partition_point(|&x| x <= key));
                let i = rng.below(model.len() + 1);
                assert_eq!(list.get(i), model.get(i).copied());
            }
        }
        st.ops(1);
    }
    assert_eq!(list.len(), model.len());
    st.done();
}

#[test]
fn extreme_keys_insert_and_remove() {
    // `remove` used to split at `key + 1`, which overflows for i64::MAX (a panic in debug
    // builds, a silent no-op in release builds).
    let mut list = SortedList::new();
    let mut model = vec![];
    for k in [i64::MAX, i64::MIN, 0, i64::MAX, i64::MIN, i64::MAX - 1, i64::MIN + 1] {
        list.insert(k);
        let at = model.partition_point(|&x| x < k);
        model.insert(at, k);
    }
    check_all(&list, &model, "after inserts");
    for k in [i64::MAX, i64::MIN, i64::MAX, i64::MAX, 7, i64::MIN + 1] {
        list.remove(k);
        let first = model.partition_point(|&x| x < k);
        if model.get(first) == Some(&k) {
            model.remove(first);
        }
        check_all(&list, &model, &format!("after remove({k})"));
        assert_eq!(list.bisect_left(k), model.partition_point(|&x| x < k));
        assert_eq!(list.bisect_right(k), model.partition_point(|&x| x <= k));
    }
}
