use std::io::{self, Read, Write, BufWriter};
use std::cmp::{min, max};
use std::iter::FromIterator;
use std::fmt::{Display, Debug};
use std::collections::{VecDeque, BTreeMap, BTreeSet, BinaryHeap, HashSet};
use std::cmp::Ordering;

const INF: i64 = std::i64::MAX;
const MOD: i64 = 998244353;


fn print_vec<T: Display>(vec: &[T], sep: &str) {
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    if let Some((first, rest)) = vec.split_first() {
        write!(out, "{}", first).unwrap();
        for item in rest {
            write!(out, "{}{}", sep, item).unwrap();
        }
    }
    writeln!(out).unwrap();
}


fn gcd(a: i64, b: i64) -> i64 {
    match b {
        0 => a,
        _ => gcd(b, a % b)
    }
}

fn solve(cin: &mut FastScan) {
    let n: usize = cin.next();
    let mut arr: Vec<i64> = cin.read_vec(n);
    arr.sort();
    print_vec(&arr, " ");
}

fn main() {
    let mut cin = FastScan::new();
    let t = 1; // cin.next::<usize>();
    for _ in 0..t {
        solve(&mut cin);
    }
}

struct FastScan<'a> {
    iter: std::str::SplitWhitespace<'a>,
}

impl<'a> FastScan<'a> {
    fn new() -> Self {
        let mut input_buffer = String::new();
        let _ = io::stdin().read_to_string(&mut input_buffer);
        let static_str: &'static str = Box::leak(input_buffer.into_boxed_str());
        FastScan {
            iter: static_str.split_whitespace(),
        }
    }

    fn next<T: std::str::FromStr>(&mut self) -> T {
        self.iter.next().and_then(|s| s.parse().ok()).unwrap()
    }

    fn read_vec<T: std::str::FromStr>(&mut self, size: usize) -> Vec<T> {
        (0..size).map(|_| self.next()).collect()
    }
}