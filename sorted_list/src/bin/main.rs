use sorted_list::SortedList;

/// Driver for the sorted list example.
pub fn run() {
    println!("Sorted List (Treap)");

    let mut sl = SortedList::new();
    sl.insert(5);
    sl.insert(3);
    sl.insert(8);
    sl.insert(1);
    sl.insert(9);
    println!("len: {}", sl.len());

    println!("bisect_left(5): {}", sl.bisect_left(5));
    println!("bisect_right(5): {}", sl.bisect_right(5));

    for i in 0..sl.len() {
        print!("{:?} ", sl.get(i).unwrap());
    }
    println!();

    sl.remove(5);
    println!("After removing 5:");
    for i in 0..sl.len() {
        print!("{:?} ", sl.get(i).unwrap());
    }
    println!();
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_smoke() {
        run();
    }

    #[test]
    fn insert_and_iterate_sorted() {
        let mut sl = SortedList::new();
        for &k in &[5, 3, 8, 1, 9] {
            sl.insert(k);
        }
        assert_eq!(sl.len(), 5);
        let collected: Vec<i64> = (0..sl.len()).map(|i| sl.get(i).unwrap()).collect();
        assert_eq!(collected, vec![1, 3, 5, 8, 9]);
    }

    #[test]
    fn bisect_left_right_consistent_with_iteration() {
        // The treap stores keys by `Ord`, so for any `key`,
        // `bisect_left(key)` must equal the count of iterated keys < key.
        let mut sl = SortedList::new();
        for &k in &[5, 3, 8, 1, 9, 3, 7] {
            sl.insert(k);
        }
        let collected: Vec<i64> = (0..sl.len()).map(|i| sl.get(i).unwrap()).collect();
        for probe in &[0i64, 1, 4, 7, 10] {
            let expected_left = collected.iter().filter(|&&x| x < *probe).count();
            let expected_right = collected.iter().filter(|&&x| x <= *probe).count();
            assert_eq!(
                sl.bisect_left(*probe),
                expected_left,
                "bisect_left({}) mismatch",
                probe
            );
            assert_eq!(
                sl.bisect_right(*probe),
                expected_right,
                "bisect_right({}) mismatch",
                probe
            );
        }
    }

    #[test]
    fn remove_and_duplicates() {
        let mut sl = SortedList::new();
        for &k in &[5, 3, 8, 3, 5] {
            sl.insert(k);
        }
        assert_eq!(sl.len(), 5);
        sl.remove(5);
        // one 5 should remain
        assert_eq!(sl.len(), 4);
        let collected: Vec<i64> = (0..sl.len()).map(|i| sl.get(i).unwrap()).collect();
        assert_eq!(collected, vec![3, 3, 5, 8]);
    }

    #[test]
    fn is_empty() {
        let mut sl = SortedList::new();
        assert!(sl.is_empty());
        sl.insert(1);
        assert!(!sl.is_empty());
    }
}