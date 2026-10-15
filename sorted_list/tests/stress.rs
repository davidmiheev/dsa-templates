//! Stress test: bucketed `SortedList` against a sorted `Vec`. Sizes go well past the bucket
//! capacity (1024), so bucket splits and removals of emptied buckets are exercised.

use sorted_list::SortedList;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

fn check_all(list: &SortedList<i64>, model: &[i64], ctx: &str) {
    assert_eq!(list.len(), model.len(), "len {ctx}");
    assert_eq!(list.is_empty(), model.is_empty(), "is_empty {ctx}");
    for (i, &v) in model.iter().enumerate() {
        assert_eq!(list.get(i), v, "get({i}) {ctx}");
    }
    assert!(list.iter().eq(model.iter().copied()), "iter {ctx}");
    assert!(catch_unwind(AssertUnwindSafe(|| list.get(model.len()))).is_err(), "get(len) must panic {ctx}");
}

fn model_insert(model: &mut Vec<i64>, k: i64) {
    let at = model.partition_point(|&x| x < k);
    model.insert(at, k);
}

fn model_remove(model: &mut Vec<i64>, k: i64) -> bool {
    let first = model.partition_point(|&x| x < k);
    if model.get(first) == Some(&k) {
        model.remove(first);
        true
    } else {
        false
    }
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
            match rng.below(6) {
                0 | 1 => {
                    list.insert(key);
                    model_insert(&mut model, key);
                }
                2 => assert_eq!(list.remove(key), model_remove(&mut model, key), "remove {ctx}"),
                3 => {
                    assert_eq!(list.bisect_left(key), model.partition_point(|&x| x < key), "bisect_left {ctx}");
                    assert_eq!(list.bisect_right(key), model.partition_point(|&x| x <= key), "bisect_right {ctx}");
                }
                4 => assert_eq!(list.contains(key), model.binary_search(&key).is_ok(), "contains {ctx}"),
                _ => {
                    if !model.is_empty() {
                        let i = rng.below(model.len());
                        assert_eq!(list.get(i), model[i], "get({i}) {ctx}");
                    }
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
fn from_iterator_matches_sorted_input() {
    let mut st = Stress::new("sorted_list/from_iter");
    for case in 0..100 * scale() {
        let mut rng = Rng::new(7_000 + case as u64);
        let n = rng.below(5_000);
        let items: Vec<i64> = (0..n).map(|_| rng.range(-1000, 1000)).collect();
        let mut model = items.clone();
        model.sort_unstable();
        let mut list: SortedList<i64> = items.into_iter().collect();
        check_all(&list, &model, &format!("from_iter case {case}"));
        // The structure must stay valid when mutated after bulk construction.
        for _ in 0..200 {
            let k = rng.range(-1000, 1000);
            if rng.chance(1, 2) {
                list.insert(k);
                model_insert(&mut model, k);
            } else {
                assert_eq!(list.remove(k), model_remove(&mut model, k));
            }
        }
        check_all(&list, &model, &format!("mutated case {case}"));
        st.case();
        st.ops(n + 200);
    }
    st.done();
}

#[test]
fn drain_and_refill_crosses_bucket_boundaries() {
    let mut st = Stress::new("sorted_list/drain");
    let mut rng = Rng::new(42);
    let mut list = SortedList::new();
    let mut model: Vec<i64> = Vec::new();
    st.case();
    for round in 0..6 * scale() {
        // Sorted, reverse-sorted and random insertion orders all force different split patterns.
        for i in 0..4_000 {
            let k = match round % 3 {
                0 => i,
                1 => 4_000 - i,
                _ => rng.range(-50, 50),
            };
            list.insert(k);
            model_insert(&mut model, k);
        }
        while !model.is_empty() {
            let k = model[rng.below(model.len())];
            assert!(list.remove(k));
            model_remove(&mut model, k);
            st.ops(1);
        }
        check_all(&list, &model, &format!("drained round {round}"));
        assert!(list.is_empty() && !list.contains(0) && list.bisect_left(5) == 0);
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
                model_insert(&mut model, key);
            }
            2 => assert_eq!(list.remove(key), model_remove(&mut model, key)),
            _ => {
                assert_eq!(list.bisect_left(key), model.partition_point(|&x| x < key));
                assert_eq!(list.bisect_right(key), model.partition_point(|&x| x <= key));
                if !model.is_empty() {
                    let i = rng.below(model.len());
                    assert_eq!(list.get(i), model[i]);
                }
            }
        }
        st.ops(1);
    }
    assert_eq!(list.len(), model.len());
    st.done();
}

#[test]
fn extreme_keys_insert_and_remove() {
    // Keys at the ends of the i64 range must work (a former `key + 1` split overflowed).
    let mut list = SortedList::new();
    let mut model = vec![];
    for k in [i64::MAX, i64::MIN, 0, i64::MAX, i64::MIN, i64::MAX - 1, i64::MIN + 1] {
        list.insert(k);
        model_insert(&mut model, k);
    }
    check_all(&list, &model, "after inserts");
    for k in [i64::MAX, i64::MIN, i64::MAX, i64::MAX, 7, i64::MIN + 1] {
        assert_eq!(list.remove(k), model_remove(&mut model, k));
        check_all(&list, &model, &format!("after remove({k})"));
        assert_eq!(list.bisect_left(k), model.partition_point(|&x| x < k));
        assert_eq!(list.bisect_right(k), model.partition_point(|&x| x <= k));
    }
}

#[test]
fn works_for_other_types() {
    let mut words: SortedList<&str> = ["pear", "apple", "fig"].into_iter().collect();
    words.insert("kiwi");
    assert_eq!(words.iter().collect::<Vec<_>>(), ["apple", "fig", "kiwi", "pear"]);
    assert_eq!((words.bisect_left("g"), words.get(3)), (2, "pear"));
}
