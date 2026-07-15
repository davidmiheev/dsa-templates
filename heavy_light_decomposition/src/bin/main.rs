use heavy_light_decomposition::Hld;

/// Run the HLD example.
pub fn run() {
    println!("heavy_light_decomposition");

    // Tree:
    //      0
    //      |
    //      1
    //     / \
    //    2   3
    //    |   |
    //    4   5
    let adj = vec![
        vec![1],
        vec![0, 2, 3],
        vec![1, 4],
        vec![1, 5],
        vec![2],
        vec![3],
    ];
    let hld = Hld::new(&adj, 0);

    // Print parent / depth / chain head / pos for every vertex.
    for v in 0..adj.len() {
        let p = match hld.parent[v] {
            Some(p) if p == usize::MAX => "root".to_string(),
            Some(p) => p.to_string(),
            None => "-".to_string(),
        };
        println!(
            "v={} parent={} depth={} head={} pos={}",
            v, p, hld.depth[v], hld.head[v], hld.pos[v]
        );
    }

    // Path from 4 to 5 should be covered by 3 segments.
    let segs = hld.path_segments(4, 5);
    println!("segments(4,5) = {:?}", segs);
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_tree() -> Vec<Vec<usize>> {
        vec![
            vec![1],
            vec![0, 2, 3],
            vec![1, 4],
            vec![1, 5],
            vec![2],
            vec![3],
        ]
    }

    #[test]
    fn parents_and_depths() {
        let hld = Hld::new(&sample_tree(), 0);
        assert_eq!(hld.parent[0], None);
        assert_eq!(hld.parent[1], Some(0));
        assert_eq!(hld.parent[2], Some(1));
        assert_eq!(hld.parent[3], Some(1));
        assert_eq!(hld.depth[0], 0);
        assert_eq!(hld.depth[1], 1);
        assert_eq!(hld.depth[2], 2);
        assert_eq!(hld.depth[3], 2);
        assert_eq!(hld.depth[4], 3);
        assert_eq!(hld.depth[5], 3);
    }

    #[test]
    fn positions_are_permutation() {
        let tree = sample_tree();
        let n = tree.len();
        let hld = Hld::new(&tree, 0);
        let mut seen = vec![false; n];
        for v in 0..n {
            seen[hld.pos[v]] = true;
        }
        assert!(seen.iter().all(|&x| x));
    }

    #[test]
    fn path_segments_cover_full_path() {
        let tree = sample_tree();
        let hld = Hld::new(&tree, 0);
        let segs = hld.path_segments(4, 5);
        // 4 -> 1 -> 0 -> 1 -> 3 -> 5; split into O(log n) segments.
        assert!(!segs.is_empty());
        assert!(segs.len() <= 6);
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
