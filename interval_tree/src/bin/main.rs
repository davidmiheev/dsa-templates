use interval_tree::{Interval, IntervalTree};

/// Run the interval-tree example.
///
/// Insert six intervals, then query an overlapping interval and print the
/// payload of every overlap.
pub fn run() {
    println!("interval_tree");

    let intervals = [
        Interval::new(15, 20, "a"),
        Interval::new(10, 30, "b"),
        Interval::new(17, 19, "c"),
        Interval::new(5,  20, "d"),
        Interval::new(12, 15, "e"),
        Interval::new(30, 40, "f"),
    ];

    let tree = IntervalTree::new(&intervals);

    // Query: anything overlapping [14, 17).
    let hits = tree.query(14, 17);
    let names: Vec<String> = hits.iter().map(|iv| iv.value.to_string()).collect();
    println!("{}", names.join(" "));
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small_tree() -> IntervalTree<&'static str> {
        let intervals = [
            Interval::new(15, 20, "a"),
            Interval::new(10, 30, "b"),
            Interval::new(17, 19, "c"),
            Interval::new(5,  20, "d"),
            Interval::new(12, 15, "e"),
            Interval::new(30, 40, "f"),
        ];
        IntervalTree::new(&intervals)
    }

    #[test]
    fn query_finds_overlapping_intervals() {
        let tree = small_tree();
        let overlaps = tree.query(14, 17);
        let names: Vec<&str> = overlaps.iter().map(|iv| iv.value).collect();
        // For half-open intervals, `[a,b)` overlaps `[qlo,qhi)` iff
        // `a < qhi && qlo < b`. Applied to query [14, 17):
        //   a [15,20): 15<17 ✓ 14<20 ✓ → overlaps
        //   b [10,30): ✓ ✓ → overlaps
        //   c [17,19): 17<17 ✗ → no overlap (boundary)
        //   d [5,20):  ✓ ✓ → overlaps
        //   e [12,15): 12<17 ✓ 14<15 ✓ → overlaps
        //   f [30,40): 30<17 ✗ → no overlap
        assert!(names.contains(&"a"), "missing a: {names:?}");
        assert!(names.contains(&"b"), "missing b: {names:?}");
        assert!(!names.contains(&"c"), "unexpectedly included c: {names:?}");
        assert!(names.contains(&"d"), "missing d: {names:?}");
        assert!(names.contains(&"e"), "missing e: {names:?}");
        assert!(!names.contains(&"f"), "unexpectedly included f: {names:?}");
    }

    #[test]
    fn empty_query_at_no_overlap() {
        let tree = small_tree();
        let overlaps = tree.query(100, 200);
        assert!(overlaps.is_empty());
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
