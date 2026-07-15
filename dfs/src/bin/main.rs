use dfs::{dfs_order, has_cycle_undirected};

/// Run the DFS example.
pub fn run() {
    println!("dfs");

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

    let order = dfs_order(&graph, 0);
    let s: Vec<String> = order.iter().map(|v| v.to_string()).collect();
    println!("{}", s.join(" "));

    // Cycle present: 0 - 1 - 2 - 0
    let cyclic: Vec<Vec<usize>> = vec![vec![1, 2], vec![0, 2], vec![0, 1]];
    // Acyclic: 0 - 1 - 2
    let acyclic: Vec<Vec<usize>> = vec![vec![1], vec![0, 2], vec![1]];
    println!("{}", has_cycle_undirected(&cyclic));
    println!("{}", has_cycle_undirected(&acyclic));
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dfs_visits_all_reachable() {
        let graph: Vec<Vec<usize>> = vec![vec![1, 2], vec![0, 3, 4], vec![0], vec![1], vec![1]];
        let order = dfs_order(&graph, 0);
        assert_eq!(order.len(), 5);
        assert_eq!(order[0], 0);
    }

    #[test]
    fn cycle_detection() {
        let cyclic: Vec<Vec<usize>> = vec![vec![1, 2], vec![0, 2], vec![0, 1]];
        let acyclic: Vec<Vec<usize>> = vec![vec![1], vec![0, 2], vec![1]];
        assert!(has_cycle_undirected(&cyclic));
        assert!(!has_cycle_undirected(&acyclic));
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
