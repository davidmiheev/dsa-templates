use centroid_decomposition::{CentroidTree, ContourSum};

/// Run the centroid decomposition example.
pub fn run() {
    println!("centroid_decomposition");

    // Path 0 -- 1 -- 2 -- 3 -- 4 -- 5 -- 6
    let edges: Vec<(usize, usize)> = (1..7).map(|i| (i - 1, i)).collect();
    let tree = CentroidTree::from_edges(7, &edges);

    println!("Centroid tree parents / depths:");
    for v in 0..tree.len() {
        let p = if tree.parent[v] == CentroidTree::NONE { "root".to_string() } else { tree.parent[v].to_string() };
        println!("v={} parent={} depth={}", v, p, tree.depth[v]);
    }

    let mut c = ContourSum::<i64>::new(tree);
    for v in 0..7 {
        c.add(v, 1 << v);
    }
    println!("sum of values at distance 1..3 from vertex 2: {}", c.sum(2, 1, 3));
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn path_centroid_is_middle() {
        let tree = CentroidTree::from_edges(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
        assert_eq!(tree.parent[2], CentroidTree::NONE);
        assert_eq!(tree.root, 2);
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
