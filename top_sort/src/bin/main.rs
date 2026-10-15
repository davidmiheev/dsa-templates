use top_sort::{is_acyclic, topological_sort};

/// Run the topological sort example.
pub fn run() {
    println!("top_sort");

    // DAG: 5 -> 3, 5 -> 2, 4 -> 2, 4 -> 1, 3 -> 1, 2 -> 0, 2 -> 3, 1 -> 0
    let graph: Vec<Vec<usize>> = vec![
        vec![],        // 0
        vec![0],       // 1
        vec![0, 3],    // 2
        vec![1],       // 3
        vec![2, 1],    // 4
        vec![3, 2],    // 5
    ];

    assert!(is_acyclic(&graph));
    let order = topological_sort(&graph);
    let s: Vec<String> = order.iter().map(|v| v.to_string()).collect();
    println!("{}", s.join(" "));

    // Cyclic: 0 -> 1 -> 0 (`topological_sort` would panic, so test it first)
    let cyclic: Vec<Vec<usize>> = vec![vec![1], vec![0], vec![]];
    println!("{}", if is_acyclic(&cyclic) { format!("{:?}", topological_sort(&cyclic)) } else { "cycle".to_string() });
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_dag_is_sorted() {
        // 0 <- 1 <- 2 <- 3 (i.e., 3 -> 2 -> 1 -> 0)
        let graph: Vec<Vec<usize>> = vec![vec![], vec![0], vec![1], vec![2]];
        let order = topological_sort(&graph);
        assert_eq!(order, vec![3, 2, 1, 0]);
    }

    #[test]
    #[should_panic(expected = "cycle")]
    fn cyclic_panics() {
        let graph: Vec<Vec<usize>> = vec![vec![1], vec![0]];
        topological_sort(&graph);
    }

    #[test]
    fn is_acyclic_detects_cycles() {
        assert!(is_acyclic(&[vec![], vec![0]]));
        assert!(!is_acyclic(&[vec![1], vec![0]]));
        assert!(!is_acyclic(&[vec![0]]));
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
