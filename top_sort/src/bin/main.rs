use top_sort::topological_sort;

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

    match topological_sort(&graph) {
        Some(order) => {
            let s: Vec<String> = order.iter().map(|v| v.to_string()).collect();
            println!("{}", s.join(" "));
        }
        None => println!("cycle"),
    }

    // Cyclic: 0 -> 1 -> 0
    let cyclic: Vec<Vec<usize>> = vec![vec![1], vec![0], vec![]];
    println!("{}", match topological_sort(&cyclic) {
        Some(o) => format!("{:?}", o),
        None => "cycle".to_string(),
    });
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
        let order = topological_sort(&graph).unwrap();
        assert_eq!(order, vec![3, 2, 1, 0]);
    }

    #[test]
    fn cyclic_returns_none() {
        let graph: Vec<Vec<usize>> = vec![vec![1], vec![0]];
        assert!(topological_sort(&graph).is_none());
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
