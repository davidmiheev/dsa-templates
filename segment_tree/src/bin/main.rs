use segment_tree::SegmentTree;

/// Run the segment tree driver.
///
/// The original example built an array `[1, 3, 5, 7, 9, 11]`,
/// queried the sum of `[1, 3]` (expecting 15), updated index 2 to 6,
/// and re-queried (expecting 16).
pub fn run() {
    let arr = vec![1, 3, 5, 7, 9, 11];

    // Create a segment tree for summation.
    let mut seg_tree = SegmentTree::new(&arr, 0, |a, b| a + b);

    // Query the sum of a range.
    // The sum of elements from index 1 to 3 (inclusive) is 3 + 5 + 7 = 15.
    println!("Sum of range [1, 3]: {}", seg_tree.query(1, 3));

    // Update a value.
    // Update the element at index 2 to be 6.
    seg_tree.update(2, 6);

    // The new array is [1, 3, 6, 7, 9, 11].
    // The new sum of elements from index 1 to 3 is 3 + 6 + 7 = 16.
    println!("Sum of range [1, 3] after update: {}", seg_tree.query(1, 3));
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
        assert_eq!(seg_tree.query(1, 3), 15);
        seg_tree.update(2, 6);
        assert_eq!(seg_tree.query(1, 3), 16);
    }

    #[test]
    fn min_query() {
        let arr = vec![5, 2, 8, 1, 9, 3];
        let seg_tree = SegmentTree::new(&arr, i64::MAX, |a, b| a.min(b));
        assert_eq!(seg_tree.query(0, 5), 1);
        assert_eq!(seg_tree.query(2, 4), 1);
    }
}