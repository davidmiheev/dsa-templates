use heavy_light_decomposition::Hld;

/// Run the HLD example.
pub fn run() {
    println!("heavy_light_decomposition");

    // 0 - 1 - {2 - 4, 3 - 5}
    let hld = Hld::from_edges(6, &[(0, 1), (1, 2), (1, 3), (2, 4), (3, 5)], 0);

    for v in 0..hld.len() {
        let p = if hld.parent[v] == Hld::NONE { "root".to_string() } else { hld.parent[v].to_string() };
        println!("v={} parent={} depth={} head={} pos={}", v, p, hld.depth[v], hld.head[v], hld.pos[v]);
    }

    println!("segments(4,5) = {:?}", hld.path_segments(4, 5));
    println!("lca(4,5) = {}, dist(4,5) = {}", hld.lca(4, 5), hld.dist(4, 5));
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_tree_basics() {
        let hld = Hld::from_edges(6, &[(0, 1), (1, 2), (1, 3), (2, 4), (3, 5)], 0);
        assert_eq!(hld.parent[0], Hld::NONE);
        assert_eq!(hld.parent[4], 2);
        assert_eq!(hld.lca(4, 5), 1);
        assert_eq!(hld.dist(4, 5), 4);
        assert_eq!(hld.jump(4, 2), 1);
        assert_eq!(hld.path_vertex(4, 5, 3), 3);
    }

    #[test]
    #[should_panic(expected = "ancestor")]
    fn jump_too_far_panics() {
        Hld::from_edges(2, &[(0, 1)], 0).jump(1, 2);
    }

    #[test]
    #[should_panic(expected = "cycle")]
    fn cycle_is_rejected() {
        Hld::new(&[vec![1, 2], vec![0, 2], vec![0, 1]], 0);
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
