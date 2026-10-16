use csr_graph::Csr;
use dijkstra::dijkstra;

/// Run the Dijkstra example.
pub fn run() {
    println!("dijkstra");
    // 0 -> 1 (4), 0 -> 2 (1), 2 -> 1 (2), 1 -> 3 (5)
    let edges = [(0, 1), (0, 2), (2, 1), (1, 3)];
    let g = Csr::directed(5, &edges);
    let sp = dijkstra(&g, &[4u64, 1, 2, 5], 0);
    for v in 0..5 {
        if sp.reachable(v) {
            println!("dist[{}] = {} via {:?}", v, sp.dist[v], sp.path_vertices(v));
        } else {
            println!("dist[{}] = unreachable", v);
        }
    }
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_graph() {
        let g = Csr::directed(3, &[(0, 1), (1, 2), (0, 2)]);
        let sp = dijkstra(&g, &[1u32, 1, 5], 0);
        assert_eq!(sp.dist, vec![0, 1, 2]);
        assert_eq!(sp.path_edges(2), vec![0, 1]);
    }

    #[test]
    #[should_panic(expected = "not reachable")]
    fn unreachable_path_panics() {
        let g = Csr::directed(2, &[]);
        dijkstra(&g, &[] as &[u32], 0).path_edges(1);
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
