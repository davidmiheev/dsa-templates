use segment_tree::SegmentTree;

/// Run the segment tree driver.
///
/// Builds `[1, 3, 5, 7, 9, 11]`, queries the sum of `[1, 4)` (15), updates index 2
/// to 6, and re-queries (16).
pub fn run() {
    let arr = vec![1, 3, 5, 7, 9, 11];

    // Create a segment tree for summation.
    let mut seg_tree = SegmentTree::new(&arr, 0, |a, b| a + b);

    // Query the sum of a half-open range.
    // The sum of elements in [1, 4) is 3 + 5 + 7 = 15.
    println!("Sum of range [1, 4): {}", seg_tree.query(1, 4));

    // Update a value.
    // Update the element at index 2 to be 6.
    seg_tree.update(2, 6);

    // The new array is [1, 3, 6, 7, 9, 11].
    // The new sum of elements in [1, 4) is 3 + 6 + 7 = 16.
    println!("Sum of range [1, 4) after update: {}", seg_tree.query(1, 4));
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_prints_expected_sums() {
        // Just smoke-test that the driver runs without panicking;
        // it prints to stdout, so we only assert no panic.
        run();
    }

    #[test]
    fn sum_query_after_update() {
        let arr = vec![1, 3, 5, 7, 9, 11];
        let mut seg_tree = SegmentTree::new(&arr, 0, |a, b| a + b);
        assert_eq!(seg_tree.query(1, 4), 15);
        seg_tree.update(2, 6);
        assert_eq!(seg_tree.query(1, 4), 16);
    }

    #[test]
    fn min_query() {
        let arr = vec![5, 2, 8, 1, 9, 3];
        let seg_tree = SegmentTree::new(&arr, i64::MAX, |a, b| a.min(b));
        assert_eq!(seg_tree.query(0, 6), 1);
        assert_eq!(seg_tree.query(2, 4), 1);
        assert_eq!(seg_tree.query(2, 3), 8);
    }

    #[test]
    fn empty_tree_and_empty_range() {
        let st = SegmentTree::new(&[] as &[i64], 0, |a, b| a + b);
        assert!(st.is_empty());
        assert_eq!(st.all_prod(), 0);
        let st = SegmentTree::new(&[4, 5], 0, |a, b| a + b);
        assert_eq!(st.query(1, 1), 0);
    }

    #[test]
    #[should_panic(expected = "invalid range")]
    fn out_of_range_panics() {
        SegmentTree::new(&[1, 2, 3], 0, |a, b| a + b).query(2, 4);
    }
}