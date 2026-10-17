//! Stress test: `BinaryTree` traversals against iterative (explicit-stack) references.

use binary_tree::BinaryTree;
use std::collections::VecDeque;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

/// `slot[i]` is `Some(value)` if node `i` has been set. Only indices `1..=n` are visited.
fn inorder(slot: &[Option<u32>], n: usize) -> Vec<u32> {
    let (mut out, mut stack, mut cur) = (vec![], vec![], 1usize);
    while cur <= n || !stack.is_empty() {
        while cur <= n {
            stack.push(cur);
            cur *= 2;
        }
        let v = stack.pop().unwrap();
        out.extend(slot[v]);
        cur = v * 2 + 1;
    }
    out
}

fn preorder(slot: &[Option<u32>], n: usize) -> Vec<u32> {
    let (mut out, mut stack) = (vec![], vec![1usize]);
    while let Some(v) = stack.pop() {
        if v > n {
            continue;
        }
        out.extend(slot[v]);
        stack.push(v * 2 + 1);
        stack.push(v * 2);
    }
    out
}

fn postorder(slot: &[Option<u32>], n: usize) -> Vec<u32> {
    // reverse of (root, right, left)
    let (mut out, mut stack) = (vec![], vec![1usize]);
    while let Some(v) = stack.pop() {
        if v > n {
            continue;
        }
        out.extend(slot[v]);
        stack.push(v * 2);
        stack.push(v * 2 + 1);
    }
    out.reverse();
    out
}

fn level_order(slot: &[Option<u32>], n: usize) -> Vec<u32> {
    let (mut out, mut q) = (vec![], VecDeque::from([1usize]));
    while let Some(v) = q.pop_front() {
        if v > n {
            continue;
        }
        out.extend(slot[v]);
        q.push_back(v * 2);
        q.push_back(v * 2 + 1);
    }
    out
}

#[test]
fn traversals_match_iterative_references() {
    let mut st = Stress::new("binary_tree/naive");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let size = rng.range(1, 80) as usize;
        let mut tree = BinaryTree::new(size);
        let mut slot: Vec<Option<u32>> = vec![None; size + 1];
        st.case();
        for _ in 0..rng.below(2 * size + 1) {
            let (i, v) = (rng.range(1, size as i64) as usize, rng.next_u64() as u32);
            tree.set(i, v);
            slot[i] = Some(v);
        }
        for i in 1..=size {
            assert_eq!(tree.get(i), slot[i].as_ref(), "get({i}) case {case}");
        }
        for n in [size, rng.range(0, size as i64) as usize, 0] {
            assert_eq!(tree.inorder(n), inorder(&slot, n), "inorder n={n} case {case}");
            assert_eq!(tree.preorder(n), preorder(&slot, n), "preorder n={n} case {case}");
            assert_eq!(tree.postorder(n), postorder(&slot, n), "postorder n={n} case {case}");
            assert_eq!(tree.level_order(n), level_order(&slot, n), "level_order n={n} case {case}");
            st.ops(4);
        }
    }
    st.done();
}

#[test]
fn index_zero_and_out_of_range_panic() {
    let mut t = BinaryTree::new(4);
    assert!(std::panic::catch_unwind(|| BinaryTree::<u8>::new(4).get(0).copied()).is_err());
    assert!(std::panic::catch_unwind(|| BinaryTree::<u8>::new(4).get(5).copied()).is_err());
    assert!(std::panic::catch_unwind(move || t.set(0, 1u8)).is_err());
    let mut t = BinaryTree::new(4);
    assert!(std::panic::catch_unwind(move || t.set(5, 1u8)).is_err());
}

#[test]
fn deep_complete_tree() {
    // 2^16 - 1 nodes: every traversal visits each node once and pre/level order start at the root.
    let mut st = Stress::new("binary_tree/large");
    let n = (1usize << 16) - 1;
    let mut tree = BinaryTree::new(n);
    for i in 1..=n {
        tree.set(i, i as u32);
    }
    st.case();
    let (pre, ino, post, lvl) = (tree.preorder(n), tree.inorder(n), tree.postorder(n), tree.level_order(n));
    assert_eq!(lvl, (1..=n as u32).collect::<Vec<_>>());
    assert_eq!((pre[0], *post.last().unwrap()), (1, 1));
    for v in [&pre, &ino, &post] {
        let mut s = v.clone();
        s.sort_unstable();
        assert_eq!(s, (1..=n as u32).collect::<Vec<_>>(), "each node exactly once");
    }
    st.ops(4 * n);
    st.done();
}
