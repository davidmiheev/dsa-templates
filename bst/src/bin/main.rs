use bst::Bst;

/// Run the BST example: insert a handful of keys, then print them in
/// sorted order along with min/max membership.
pub fn run() {
    println!("bst");

    let mut tree = Bst::new();
    for k in [5, 3, 7, 1, 4, 6, 8] {
        tree.insert(k);
    }

    // Inserting a duplicate must report `false` and not change structure.
    let added = tree.insert(4);

    let sorted = tree.inorder()
        .into_iter()
        .map(|k| k.to_string())
        .collect::<Vec<_>>()
        .join(" ");

    let contains_4 = tree.contains(4);
    let min = tree.min().map(|k| k.to_string()).unwrap_or_default();
    let max = tree.max().map(|k| k.to_string()).unwrap_or_default();

    println!("{}", sorted);
    println!("{} {} {} {}", added, contains_4, min, max);
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_returns_new_for_unique_and_false_for_duplicate() {
        let mut tree = Bst::new();
        assert!(tree.insert(3));
        assert!(tree.insert(1));
        assert!(tree.insert(5));
        assert!(!tree.insert(3)); // duplicate
    }

    #[test]
    fn inorder_yields_sorted_keys() {
        let mut tree = Bst::new();
        for k in [5, 3, 7, 1, 4, 6, 8] {
            tree.insert(k);
        }
        assert_eq!(tree.inorder(), vec![1, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn contains_distinguishes_present_and_absent() {
        let mut tree = Bst::new();
        for k in [2, 1, 3] {
            tree.insert(k);
        }
        assert!(tree.contains(2));
        assert!(!tree.contains(42));
    }

    #[test]
    fn min_and_max_track_extrema() {
        let mut tree = Bst::new();
        assert_eq!(tree.min(), None);
        for k in [5, 3, 7, 1, 4, 6, 8] {
            tree.insert(k);
        }
        assert_eq!(tree.min(), Some(1));
        assert_eq!(tree.max(), Some(8));
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
