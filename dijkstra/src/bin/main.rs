use dijkstra::dijkstra;

/// Run the Dijkstra example.
pub fn run() {
    println!("dijkstra");

    // Directed weighted graph:
    //   0 -> 1 (4), 0 -> 2 (1)
    //   1 -> 3 (1)
    //   2 -> 1 (2), 2 -> 3 (5)
    //   3 -> (none)
    let graph: Vec<Vec<(usize, u64)>> = vec![
        vec![(1, 4), (2, 1)],
        vec![(3, 1)],
        vec![(1, 2), (3, 5)],
        vec![],
    ];

    let dist = dijkstra(&graph, 0);
    println!("{}", dist.get(&0).unwrap()); // 0
    println!("{}", dist.get(&1).unwrap()); // 3
    println!("{}", dist.get(&2).unwrap()); // 1
    println!("{}", dist.get(&3).unwrap()); // 4
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shortest_paths_classic_example() {
        let graph: Vec<Vec<(usize, u64)>> = vec![
            vec![(1, 4), (2, 1)],
            vec![(3, 1)],
            vec![(1, 2), (3, 5)],
            vec![],
        ];
        let dist = dijkstra(&graph, 0);
        assert_eq!(dist[&0], 0);
        assert_eq!(dist[&1], 3); // 0 -> 2 -> 1 (cheaper than 0 -> 1)
        assert_eq!(dist[&2], 1);
        assert_eq!(dist[&3], 4); // 0 -> 2 -> 1 -> 3
    }

    #[test]
    fn unreachable_vertex_absent() {
        let graph: Vec<Vec<(usize, u64)>> = vec![vec![(1, 1)], vec![], vec![]];
        let dist = dijkstra(&graph, 0);
        assert!(!dist.contains_key(&2));
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
