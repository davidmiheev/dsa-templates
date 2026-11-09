//! Stress test: Deque against std::collections::VecDeque.

use deque::Deque;
use std::collections::VecDeque;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_vec_deque() {
    let mut st = Stress::new("deque");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let (mut a, mut b): (Deque<u32>, VecDeque<u32>) = (Deque::with_capacity(rng.below(5)), VecDeque::new());
        for _ in 0..rng.range(1, 300) {
            match rng.below(6) {
                0 => {
                    let x = rng.next_u64() as u32;
                    a.push_front(x);
                    b.push_front(x);
                }
                1 | 2 => {
                    let x = rng.next_u64() as u32;
                    a.push_back(x);
                    b.push_back(x);
                }
                3 if !b.is_empty() => assert_eq!(a.pop_front(), b.pop_front().unwrap()),
                4 if !b.is_empty() => assert_eq!(a.pop_back(), b.pop_back().unwrap()),
                _ if !b.is_empty() => {
                    let i = rng.below(b.len());
                    assert_eq!(a[i], b[i]);
                    a[i] ^= 1;
                    b[i] ^= 1;
                }
                _ => {}
            }
            assert_eq!(a.len(), b.len());
            st.ops(1);
        }
        assert!(a.iter().eq(b.iter().copied()));
        st.case();
    }
    st.done();
}

#[test]
#[should_panic]
fn pop_from_empty_panics() {
    Deque::<u8>::new().pop_back();
}
