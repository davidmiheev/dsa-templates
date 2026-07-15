use bfs::{bfs_distances, bfs_order};

/// Run the BFS example.
pub fn run() {
    println!("bfs");

    // Graph:
    //      0 -- 1 -- 3
    //      |    |
    //      2    4
    let graph: Vec<Vec<usize>> = vec![
        vec![1, 2],       // 0
        vec![0, 3, 4],    // 1
        vec![0],          // 2
        vec![1],          // 3
        vec![1],          // 4
    ];

    let dist = bfs_distances(&graph, 0);
    let ordered: Vec<String> = (0..5)
        .map(|v| match dist.get(&v) {
            Some(d) => format!("{}:{}", v, d),
            None => format!("{}:inf", v),
        })
        .collect();
    println!("{}", ordered.join(" "));

    let order = bfs_order(&graph, 0);
    let order_str: Vec<String> = order.iter().map(|v| v.to_string()).collect();
    println!("{}", order_str.join(" "));
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distances_from_source() {
        let graph: Vec<Vec<usize>> = vec![
            vec![1, 2],       // 0
            vec![0, 3, 4],    // 1
            vec![0],          // 2
            vec![1],          // 3
            vec![1],          // 4
        ];
        let dist = bfs_distances(&graph, 0);
        assert_eq!(dist[&0], 0);
        assert_eq!(dist[&1], 1);
        assert_eq!(dist[&2], 1);
        assert_eq!(dist[&3], 2);
        assert_eq!(dist[&4], 2);
        assert_eq!(dist.len(), 5);
    }

    #[test]
    fn disconnected_vertex_unreachable() {
        // 0 -- 1   and isolated 2
        let graph: Vec<Vec<usize>> = vec![vec![1], vec![0], vec![]];
        let dist = bfs_distances(&graph, 0);
        assert!(!dist.contains_key(&2));
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
