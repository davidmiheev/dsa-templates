//! Stress test: FastMap against std::collections::HashMap under random operation sequences.

use hash_map::FastMap;
use std::collections::HashMap;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_std_hash_map() {
    let mut st = Stress::new("hash_map");
    for case in 0..300 * scale() {
        let mut rng = Rng::new(case as u64);
        let universe = [4usize, 40, 1000][rng.below(3)];
        let big = rng.chance(1, 2);
        let key = |x: usize| if big { (x as u64).wrapping_mul(0x1234_5678_9ABC_DEF1) } else { x as u64 };
        let (mut a, mut b) = (FastMap::new(), HashMap::new());
        for _ in 0..rng.range(1, 600) {
            let k = key(rng.below(universe));
            match rng.below(5) {
                0 | 1 => {
                    let v = rng.next_u64();
                    a.insert(k, v);
                    b.insert(k, v);
                }
                2 => assert_eq!(a.remove(k), b.remove(&k).is_some()),
                3 => {
                    *a.entry(k, 7) += 1;
                    *b.entry(k).or_insert(7) += 1;
                }
                _ => {
                    assert_eq!(a.contains(k), b.contains_key(&k));
                    assert_eq!(*a.get_or(k, &u64::MAX), b.get(&k).copied().unwrap_or(u64::MAX));
                }
            }
            assert_eq!(a.len(), b.len());
            st.ops(1);
        }
        let mut items: Vec<(u64, u64)> = a.iter().map(|(k, v)| (k, *v)).collect();
        items.sort();
        let mut expect: Vec<(u64, u64)> = b.into_iter().collect();
        expect.sort();
        assert_eq!(items, expect, "case {case}");
        st.case();
    }
    st.done();
}

#[test]
#[should_panic]
fn missing_key_panics() {
    let m: FastMap<u32> = FastMap::new();
    m.get(5);
}
