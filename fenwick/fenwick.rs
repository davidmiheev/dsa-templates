use std::ops::{AddAssign, Sub};

/// Fenwick Tree (Binary Indexed Tree)
///
/// Fenwick Tree is a data structure that can efficiently update elements
/// and calculate prefix sums in a table of numbers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FenwickTree<T> {
    size: usize,
    tree: Vec<T>,
}

impl<T> FenwickTree<T>
where
    T: Copy + Default + AddAssign + Sub<Output = T>,
{
    pub fn new(nums: &[T]) -> Self {
        let size = nums.len();
        let mut tree = vec![T::default(); size + 1];

        let mut ft = Self { size, tree };

        for (i, &val) in nums.iter().enumerate() {
            ft.update(i, val);
        }
        ft
    }

    pub fn update(&mut self, index: usize, delta: T) {
        if index >= self.size {
            panic!("Index out of bounds");
        }

        let mut i = index + 1;
        loop {
            if i > self.size { break }
            self.tree[i] += delta;
            i += i & (!i + 1); // equivalent to i & -i
        }
    }

    pub fn pref(&self, index: usize) -> T {
        if index >= self.size {
            panic!("Index out of bounds");
        }
        let mut result = T::default();
        let mut i = index + 1;
        loop {
            if i <= 0 { break }
            result += self.tree[i];
            i -= i & (!i + 1); // equivalent to i & -i
        }
        result
    }

    pub fn sum_range(&self, left: usize, right: usize) -> T {
        if right < left || right >= self.size {
            panic!("Invalid range");
        }
        match left == 0 {
            true => self.pref(right),
            false => self.pref(right) - self.pref(left - 1)
        }
    }
}

pub fn main() {
    println!("fenwick");

    let nums = [1, 3, 5, 7, 9, 11];
    let mut fenwick = FenwickTree::new(&nums);
    
    println!("{}", fenwick.sum_range(0, 5));
    fenwick.update(1, -2);
    println!("{}", fenwick.sum_range(0, 5));

    println!("{:?}", fenwick.tree);
}