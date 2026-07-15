use multiset::MultiSet;

/// Driver for the multiset example.
pub fn run() {
    let mut ms = MultiSet::from_iter(vec![1, 2, 2, 3, 4, 4, 4]);

    ms.insert(1);
    ms.insert(2);
    ms.insert(2);
    ms.insert(3);

    println!("Multiset: {:?}", ms);
    println!("Count of 2: {}", ms.count(&2));

    if let Some(first) = ms.first() {
        println!("First element: {}", first);
    }

    if let Some(last) = ms.last() {
        println!("Last element: {}", last);
    }

    ms.remove(&2);
    println!("After removing one occurrence of 2: {:?}", ms);

    ms.remove_all(&4);
    println!("After removing all occurrences of 4: {:?}", ms);

    if let Some(first) = ms.first() {
        println!("First element: {}", first);
    }

    if let Some(last) = ms.last() {
        println!("Last element: {}", last);
    }
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
    fn insert_remove_count() {
        let mut ms = MultiSet::from_iter(vec![1, 2, 2, 3, 4, 4, 4]);
        assert_eq!(ms.count(&2), 2);
        assert_eq!(ms.count(&4), 3);
        ms.insert(2);
        assert_eq!(ms.count(&2), 3);
        assert!(ms.remove(&2));
        assert_eq!(ms.count(&2), 2);
        assert!(ms.remove_all(&4));
        assert_eq!(ms.count(&4), 0);
        assert!(!ms.contains(&4));
    }

    #[test]
    fn first_last_len() {
        let ms = MultiSet::from_iter(vec![3, 1, 4, 1, 5, 9, 2, 6]);
        assert_eq!(ms.first(), Some(1));
        assert_eq!(ms.last(), Some(9));
        assert_eq!(ms.len(), 8);
        assert!(!ms.is_empty());
    }
}