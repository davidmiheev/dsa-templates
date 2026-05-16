use std::io;
use std::cmp::{min, max};
use std::iter::FromIterator;
use std::fmt::{Display, Debug};
use std::collections::{VecDeque, BTreeMap, BTreeSet, BinaryHeap, HashSet};
use std::cmp::Ordering;
//use itertools::Itertools;

const INF: i64 = std::i64::MAX;
const MOD: i64 = 998244353;

fn read_line_as<T: std::str::FromStr>() -> T {
    let mut line = String::new();
    io::stdin().read_line(&mut line).unwrap();
    line.trim().parse().ok().unwrap()
}
 
fn read_line_as_vec<T: std::str::FromStr>() -> Vec<T> {
    let mut line = String::new();
    io::stdin().read_line(&mut line).unwrap();
    line.split_whitespace()
        .map(|s| s.parse().ok().unwrap())
        .collect()
}
 
fn print_vec<T: Display>(vec: &[T]) {
    for (i, item) in vec.iter().enumerate() {
        print!("{}", item);
        if i < vec.len() - 1 { print!(" ") }
        
    }
    println!();
}
 
fn gcd(a: i64, b: i64) -> i64 {
    match b {
        0 => a,
        _ => gcd(b, a % b)
    }
}
 

fn power(mut base: i128, mut exp: i128) -> i128 {

    let mut res: i128 = 1;
    base %= MOD as i128;

    while exp > 0 {
        if exp % 2 == 1 {
            res = (res * base) % MOD as i128;
        }

        exp >>= 1;

        base = (base * base) % MOD as i128;
    }

    res
}

#[inline(always)]
fn solve() {
    let n = read_line_as::<usize>();
    println!("{}", n);
}


fn main() {
    let t: i32 = 1; //. read_line_as();
    for _ in 0..t {
        solve()
    }
}


struct UnionFind {
    root: Vec<usize>,
    size: Vec<usize>,
    components: usize,
}


impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            root: (0..=n).collect(),
            size: vec![1; n + 1],
            components: n,
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.root[x] == x {
            x
        } else {
            let r = self.root[x];
            self.root[x] = self.find(r);
            self.root[x]
        }
    }

    fn union(&mut self, x: usize, y: usize) -> bool {
        let root_x = self.find(x);
        let root_y = self.find(y);
        if root_x == root_y {
            return false;
        }

        if self.size[root_x] > self.size[root_y] {
            self.size[root_x] += self.size[root_y];
            self.root[root_y] = root_x;
        } else {
            self.size[root_y] += self.size[root_x];
            self.root[root_x] = root_y;
        }

        self.components -= 1;
        true
    }

    fn are_connected(&mut self, x: usize, y: usize) -> bool {
        self.find(x) == self.find(y)
    }
}

