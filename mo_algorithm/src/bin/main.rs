use mo_algorithm::{count_distinct_in_range, Query};

/// Run the Mo's algorithm example.
pub fn run() {
    println!("mo_algorithm");

    // Array: [1, 1, 2, 1, 3, 2, 3]
    let a = vec![1, 1, 2, 1, 3, 2, 3];
    let queries = vec![
        Query { lo: 0, hi: 5 }, // {1, 1, 2, 1, 3} -> 3 distinct
        Query { lo: 2, hi: 7 }, // {2, 1, 3, 2, 3} -> 3 distinct
        Query { lo: 0, hi: 7 }, // {1, 1, 2, 1, 3, 2, 3} -> 3 distinct
    ];

    let answers = count_distinct_in_range(&a, &queries);
    let answer_str: Vec<String> = answers.iter().map(|v| v.to_string()).collect();
    println!("{}", answer_str.join(" "));
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_smoke() {
        run();
    }

    #[test]
    fn distinct_counts_basic() {
        let a = vec![1, 1, 2, 1, 3];
        let queries = vec![Query { lo: 0, hi: 3 }, Query { lo: 2, hi: 5 }];
        let answers = count_distinct_in_range(&a, &queries);
        assert_eq!(answers, vec![2, 3]);
    }
}
