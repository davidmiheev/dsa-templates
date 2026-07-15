use binary_tree::BinaryTree;

/// Run the binary-tree example.
///
/// Builds the tree
///
/// ```text
///         1
///        / \
///       2   3
///      / \
///     4   5
/// ```
///
/// and prints the four traversal orderings.
pub fn run() {
    println!("binary_tree");

    // Indices in the heap convention: root = 1, left child = 2i, right child = 2i+1.
    let mut tree: BinaryTree<i32> = BinaryTree::new(7);
    tree.set(1, 1);
    tree.set(2, 2);
    tree.set(3, 3);
    tree.set(4, 4);
    tree.set(5, 5);

    let fmt = |v: &[i32]| -> String {
        v.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(" ")
    };

    println!("{}", fmt(&tree.preorder(7)));
    println!("{}", fmt(&tree.inorder(7)));
    println!("{}", fmt(&tree.postorder(7)));
    println!("{}", fmt(&tree.level_order(7)));
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_tree() -> BinaryTree<i32> {
        let mut tree: BinaryTree<i32> = BinaryTree::new(7);
        tree.set(1, 1);
        tree.set(2, 2);
        tree.set(3, 3);
        tree.set(4, 4);
        tree.set(5, 5);
        tree
    }

    #[test]
    fn preorder_root_left_right() {
        let tree = sample_tree();
        assert_eq!(tree.preorder(7), vec![1, 2, 4, 5, 3]);
    }

    #[test]
    fn inorder_left_root_right() {
        let tree = sample_tree();
        assert_eq!(tree.inorder(7), vec![4, 2, 5, 1, 3]);
    }

    #[test]
    fn postorder_left_right_root() {
        let tree = sample_tree();
        assert_eq!(tree.postorder(7), vec![4, 5, 2, 3, 1]);
    }

    #[test]
    fn level_order_top_down() {
        let tree = sample_tree();
        assert_eq!(tree.level_order(7), vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
