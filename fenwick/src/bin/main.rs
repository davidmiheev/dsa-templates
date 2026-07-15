use fenwick::FenwickTree;

/// Run the Fenwick Tree example.
pub fn run() {
    println!("fenwick");

    let nums = [1, 3, 5, 7, 9, 11];
    let mut fenwick = FenwickTree::new(&nums);

    println!("{}", fenwick.sum_range(0, 5));
    fenwick.update(1, -2);
    println!("{}", fenwick.sum_range(0, 5));
}

fn main() {
    run();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sum_range_after_update() {
        let nums = [1, 3, 5, 7, 9, 11];
        let mut ft = FenwickTree::new(&nums);
        assert_eq!(ft.sum_range(0, 5), 36);
        ft.update(1, -2);
        assert_eq!(ft.sum_range(0, 5), 34);
    }

    #[test]
    #[should_panic(expected = "Index out of bounds")]
    fn update_out_of_bounds_panics() {
        let nums = [1, 2, 3];
        let mut ft = FenwickTree::new(&nums);
        ft.update(99, 1);
    }

    #[test]
    fn pref_matches_naive() {
        let nums = [1, 3, 5, 7, 9, 11];
        let ft = FenwickTree::new(&nums);
        for i in 0..nums.len() {
            let naive: i64 = nums[..=i].iter().sum();
            assert_eq!(ft.pref(i), naive);
        }
    }
}