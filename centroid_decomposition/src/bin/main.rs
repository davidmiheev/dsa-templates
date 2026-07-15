use centroid_decomposition::build_centroid_tree;

/// Run the centroid decomposition example.
pub fn run() {
    println!("centroid_decomposition");

    // Path 0 -- 1 -- 2 -- 3 -- 4 -- 5 -- 6
    let adj = vec![
        vec![1],
        vec![0, 2],
        vec![1, 3],
        vec![2, 4],
        vec![3, 5],
        vec![4, 6],
        vec![5],
    ];
    let (parent, depth) = build_centroid_tree(&adj);

    println!("Centroid tree parents / depths:");
    for v in 0..adj.len() {
        let p = match parent[v] {
            Some(p) => p.to_string(),
            None => "root".to_string(),
        };
        println!("v={} parent={} depth={}", v, p, depth[v]);
    }
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn path_centroid_is_middle() {
        // 0 -- 1 -- 2 -- 3 -- 4
        let adj = vec![vec![1], vec![0, 2], vec![1, 3], vec![2, 4], vec![3]];
        let (parent, _depth) = build_centroid_tree(&adj);
        // Middle vertex 2 has no parent.
        assert_eq!(parent[2], None);
    }

    #[test]
    fn star_centroid_is_center() {
        // 0 connected to 1, 2, 3, 4
        let adj = vec![
            vec![1, 2, 3, 4],
            vec![0],
            vec![0],
            vec![0],
            vec![0],
        ];
        let (parent, depth) = build_centroid_tree(&adj);
        // Vertex 0 is the centroid root.
        assert_eq!(parent[0], None);
        assert_eq!(depth[0], 0);
    }

    #[test]
    fn depth_is_bounded_by_log_n() {
        // Build a binary tree of depth 3 (15 vertices).
        let n = 15;
        let mut adj = vec![vec![]; n];
        for i in 1..n {
            adj[i].push((i - 1) / 2);
            adj[(i - 1) / 2].push(i);
        }
        let (_, depth) = build_centroid_tree(&adj);
        // log2(15) ~= 3.9; allow up to 6 for safety across shapes.
        assert!(*depth.iter().max().unwrap() <= 6);
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
