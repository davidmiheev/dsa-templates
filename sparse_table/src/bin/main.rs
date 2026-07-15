use sparse_table::{gcd_sparse_table, max_sparse_table, min_sparse_table};

/// Run the Sparse Table example.
pub fn run() {
    println!("sparse_table");

    let a = vec![1, 3, 2, 5, 4, 0, 7, 6];

    let min_st = min_sparse_table(&a);
    let max_st = max_sparse_table(&a);
    let gcd_st = gcd_sparse_table(&a);

    // Range [1, 5) = [3, 2, 5, 4]
    let mn = min_st.query(1, 5).unwrap();
    let mx = max_st.query(1, 5).unwrap();
    let g = gcd_st.query(1, 5).unwrap();
    println!("min={} max={} gcd={}", mn, mx, g);
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn min_query_basic() {
        let a = vec![1, 3, 2, 5, 4];
        let st = min_sparse_table(&a);
        assert_eq!(st.query(0, 5), Some(1));
        assert_eq!(st.query(1, 4), Some(2));
        assert_eq!(st.query(3, 5), Some(4));
    }

    #[test]
    fn max_query_basic() {
        let a = vec![1, 3, 2, 5, 4];
        let st = max_sparse_table(&a);
        assert_eq!(st.query(0, 5), Some(5));
        assert_eq!(st.query(1, 4), Some(5));
    }

    #[test]
    fn gcd_query_basic() {
        let a = vec![12, 18, 24, 36];
        let st = gcd_sparse_table(&a);
        assert_eq!(st.query(0, 4), Some(6));
        assert_eq!(st.query(1, 3), Some(6));
    }

    #[test]
    fn out_of_bounds_returns_none() {
        let a = vec![1, 2, 3];
        let st = min_sparse_table(&a);
        assert_eq!(st.query(0, 0), None);
        assert_eq!(st.query(2, 10), None);
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
