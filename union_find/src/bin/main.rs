use union_find::UnionFind;

/// Driver for the Union-Find example.
///
/// Builds a Union-Find over `0..=10`, runs four `union` calls (sum as the
/// combine op), and prints representative / character / size / connectivity.
pub fn run() {
    // Example usage
    let mut uf = UnionFind::new(10, |a, b| a + b, 0);
    uf.union(1, 2, 5);
    uf.union(2, 3, 10);
    uf.union(4, 5, 15);
    uf.union(6, 5, 20);
    println!("Root of 1: {}", uf.get_root(1));
    println!("Character of 1: {}", uf.get_character(1));
    println!("Size of component containing 1: {}", uf.get_size(1));
    println!("Are 1 and 3 connected? {}", uf.are_connected(1, 3));
    println!("Number of components: {}", uf.get_components());
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
    fn union_connects_components() {
        let mut uf = UnionFind::new(5, |a, b| a + b, 0);
        assert_eq!(uf.get_components(), 5);
        assert!(uf.union(1, 2, 0));
        assert!(uf.are_connected(1, 2));
        assert_eq!(uf.get_components(), 4);
        assert!(!uf.union(1, 2, 0));
        assert_eq!(uf.get_components(), 4);
    }

    #[test]
    fn character_combines_via_op() {
        // XOR combines characters.
        let mut uf = UnionFind::new(4, |a, b| a ^ b, 0);
        uf.union(1, 2, 7);
        assert_eq!(uf.get_character(1), 7);
        uf.union(2, 3, 3);
        assert_eq!(uf.get_character(1), 7 ^ 3);
    }
}