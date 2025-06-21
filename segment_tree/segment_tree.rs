/// A generic segment tree.
///
/// `T` is the type of the data in the array.
/// `F` is the type of the merge function.
pub struct SegmentTree<T, F> {
    size: usize,
    tree: Vec<T>,
    default: T,
    merge: F,
}

impl<T, F> SegmentTree<T, F>
where
    T: Clone + Copy,
    F: Fn(T, T) -> T,
{
    pub fn new(arr: &[T], default: T, merge: F) -> Self {
        let size = arr.len();
        let mut tree = vec![default; 4 * size];
        let mut seg_tree = SegmentTree {
            size,
            tree,
            default,
            merge,
        };
        seg_tree.build(arr, 0, 0, size - 1);
        seg_tree
    }

    fn build(&mut self, arr: &[T], node: usize, start: usize, end: usize) {
        if start == end {
            self.tree[node] = arr[start];
            return;
        }
        let mid = (start + end) / 2;
        self.build(arr, 2 * node + 1, start, mid);
        self.build(arr, 2 * node + 2, mid + 1, end);
        self.tree[node] = (self.merge)(self.tree[2 * node + 1], self.tree[2 * node + 2]);
    }

    pub fn query(&self, l: usize, r: usize) -> T {
        self.query_recursive(0, 0, self.size - 1, l, r)
    }

    fn query_recursive(&self, node: usize, start: usize, end: usize, l: usize, r: usize) -> T {
        if r < start || end < l {
            return self.default;
        }
        if l <= start && end <= r {
            return self.tree[node];
        }
        let mid = (start + end) / 2;
        let p1 = self.query_recursive(2 * node + 1, start, mid, l, r);
        let p2 = self.query_recursive(2 * node + 2, mid + 1, end, l, r);
        (self.merge)(p1, p2)
    }

    pub fn update(&mut self, idx: usize, val: T) {
        self.update_recursive(0, 0, self.size - 1, idx, val);
    }

    fn update_recursive(&mut self, node: usize, start: usize, end: usize, idx: usize, val: T) {
        if start == end {
            self.tree[node] = val;
            return;
        }
        let mid = (start + end) / 2;
        if start <= idx && idx <= mid {
            self.update_recursive(2 * node + 1, start, mid, idx, val);
        } else {
            self.update_recursive(2 * node + 2, mid + 1, end, idx, val);
        }
        self.tree[node] = (self.merge)(self.tree[2 * node + 1], self.tree[2 * node + 2]);
    }
}

fn main() {
    let arr = vec![1, 3, 5, 7, 9, 11];
    
    // Create a segment tree for summation.
    let mut seg_tree = SegmentTree::new(&arr, 0, |a, b| a + b);

    // Query the sum of a range.
    // The sum of elements from index 1 to 3 (inclusive) is 3 + 5 + 7 = 15.
    println!("Sum of range [1, 3]: {}", seg_tree.query(1, 3)); 

    // Update a value.
    // Update the element at index 2 to be 6.
    seg_tree.update(2, 6);

    // The new array is [1, 3, 6, 7, 9, 11].
    // The new sum of elements from index 1 to 3 is 3 + 6 + 7 = 16.
    println!("Sum of range [1, 3] after update: {}", seg_tree.query(1, 3));
}