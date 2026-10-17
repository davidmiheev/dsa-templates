use sparse_table::{gcd_sparse_table, max_sparse_table, min_sparse_table};

/// Run the Sparse Table example.
pub fn run() {
    println!("sparse_table");

    let a = vec![1, 3, 2, 5, 4, 0, 7, 6];

    let min_st = min_sparse_table(&a);
    let max_st = max_sparse_table(&a);
    let gcd_st = gcd_sparse_table(&a);

    // Range [1, 5) = [3, 2, 5, 4]
    let mn = min_st.query(1, 5);
    let mx = max_st.query(1, 5);
    let g = gcd_st.query(1, 5);
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
        assert_eq!(st.query(0, 5), 1);
        assert_eq!(st.query(1, 4), 2);
        assert_eq!(st.query(3, 5), 4);
    }

    #[test]
    fn max_query_basic() {
        let a = vec![1, 3, 2, 5, 4];
        let st = max_sparse_table(&a);
        assert_eq!(st.query(0, 5), 5);
        assert_eq!(st.query(1, 4), 5);
    }

    #[test]
    fn gcd_query_basic() {
        let a = vec![12, 18, 24, 36];
        let st = gcd_sparse_table(&a);
        assert_eq!(st.query(0, 4), 6);
        assert_eq!(st.query(1, 3), 6);
    }

    #[test]
    #[should_panic(expected = "invalid range")]
    fn empty_range_panics() {
        min_sparse_table(&[1, 2, 3]).query(1, 1);
    }

    #[test]
    #[should_panic(expected = "invalid range")]
    fn out_of_bounds_panics() {
        min_sparse_table(&[1, 2, 3]).query(2, 10);
    }

    #[test]
    fn run_smoke() {
        run();
    }
}
