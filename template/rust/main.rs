use std::io::{self, Read, BufRead, Write, BufWriter, StdinLock, StdoutLock};
use std::cmp::{min, max};
use std::iter::FromIterator;
use std::fmt::{Display, Debug};
use std::collections::{VecDeque, BTreeMap, BTreeSet, BinaryHeap, HashSet};
use std::cmp::Ordering;
use ac_library::{ModInt998244353 as mint, Modulus, StaticModInt};

const INF: i64 = std::i64::MAX / 2;
const MOD: i64 = 998244353;
// Set to true for interactive problems: FastScan then reads stdin lazily
// (line by line) and Out flushes after every println/print_vec.
const INTERACTIVE: bool = false;


fn gcd(a: i64, b: i64) -> i64 {
    match b {
        0 => a,
        _ => gcd(b, a % b)
    }
}

fn solve(cin: &mut FastScan, out: &mut Out) {
    let n: usize = cin.next();
    let q: usize = cin.next();

    let a: Vec<i64> = cin.read_vec(n);
    let b: Vec<i64> = cin.read_vec(n);

    let hub_dist: Vec<i64> = {
        let a3: Vec<i64> = a.iter().cycle().take(3 * n).cloned().collect();
        let b3: Vec<i64> = b.iter().cycle().take(3 * n).cloned().collect();

        let mut p = vec![0i64; 3 * n + 1];
        for k in 0..3 * n {
            p[k + 1] = p[k] + a3[k];
        }

        let mut fwd = vec![0i64; 3 * n];
        {
            let mut dq: VecDeque<usize> = VecDeque::new();
            for pos in 0..3 * n {
                let key = |j: usize| b3[j] - p[j];
                while let Some(&back) = dq.back() {
                    if key(back) >= key(pos) {
                        dq.pop_back();
                    } else {
                        break;
                    }
                }
                dq.push_back(pos);
                while *dq.front().unwrap() + n <= pos {
                    dq.pop_front();
                }
                let best = *dq.front().unwrap();
                fwd[pos] = p[pos] + (b3[best] - p[best]);
            }
        }

        let mut bwd = vec![0i64; 3 * n];
        {
            let mut dq: VecDeque<usize> = VecDeque::new();
            for pos in (0..3 * n).rev() {
                let key = |j: usize| b3[j] + p[j];
                while let Some(&back) = dq.back() {
                    if key(back) >= key(pos) {
                        dq.pop_back();
                    } else {
                        break;
                    }
                }
                dq.push_back(pos);
                while *dq.front().unwrap() >= pos + n {
                    dq.pop_front();
                }
                let best = *dq.front().unwrap();
                bwd[pos] = (b3[best] + p[best]) - p[pos];
            }
        }

        (0..n).map(|i| fwd[i + n].min(bwd[i + n])).collect()
    };

    let mut pre = vec![0i64; n + 1];
    for i in 0..n {
        pre[i + 1] = pre[i] + a[i];
    }
    let total: i64 = pre[n];

    let mut ans = vec![];
    for _ in 0..q {
        let l: usize = cin.next::<usize>() - 1;
        let r: usize = cin.next::<usize>() - 1;
        let (lo, hi) = if l < r { (l, r) } else { (r, l) };

        if hi == n {
            ans.push(hub_dist[lo]);
        } else {
            let cw = pre[hi] - pre[lo];
            let ccw = total - cw;
            let via_hub = hub_dist[lo] + hub_dist[hi];
            ans.push(cw.min(ccw).min(via_hub));
        }
    }

    out.print_vec(&ans, "\n");
}


fn main() {
    let mut cin = FastScan::new(INTERACTIVE);
    let mut out = Out::new(INTERACTIVE);
    // Interactive quick-start (INTERACTIVE = true):
    //   out.println(format_args!("? {} {}", a, b)); // written and flushed at once
    //   let reply: i64 = cin.next();                // blocks only until the reply line
    let t = 1; // cin.next::<usize>();
    for _ in 0..t {
        solve(&mut cin, &mut out);
    }
}

/// Whitespace-separated token reader over stdin.
/// - `interactive = false`: reads all of stdin up front (fastest; needs EOF, so not
///   usable when a judge answers your queries).
/// - `interactive = true`: reads one line at a time, so it only blocks until the next
///   reply line arrives. Tokens may share a line or span lines in both modes.
struct FastScan {
    stdin: StdinLock<'static>,
    interactive: bool,
    buf: String,
    pos: usize,
}

impl FastScan {
    fn new(interactive: bool) -> Self {
        let mut s = FastScan { stdin: io::stdin().lock(), interactive, buf: String::new(), pos: 0 };
        if !interactive {
            s.stdin.read_to_string(&mut s.buf).unwrap();
        }
        s
    }

    fn token(&mut self) -> &str {
        loop {
            let bytes = self.buf.as_bytes();
            let mut i = self.pos;
            while i < bytes.len() && bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            if i < bytes.len() {
                let start = i;
                while i < bytes.len() && !bytes[i].is_ascii_whitespace() {
                    i += 1;
                }
                self.pos = i;
                return &self.buf[start..i];
            }
            self.buf.clear();
            self.pos = 0;
            if !self.interactive || self.stdin.read_line(&mut self.buf).unwrap() == 0 {
                panic!("Unexpected end of input: check number of testcases or input format");
            }
        }
    }

    fn next<T: std::str::FromStr>(&mut self) -> T {
        let t = self.token();
        t.parse().ok().unwrap_or_else(|| panic!("Cannot parse {:?} as {}", t, std::any::type_name::<T>()))
    }

    fn read_vec<T: std::str::FromStr>(&mut self, size: usize) -> Vec<T> {
        (0..size).map(|_| self.next()).collect()
    }
}


/// Buffered stdout shared by the whole program (create once in `main`, pass as `&mut Out`).
/// One lock + one big buffer instead of a fresh lock/BufWriter per call, so printing
/// inside a loop (e.g. once per query) stays cheap. Flushed on drop.
/// With `interactive = true`, `println` and `print_vec` also flush, so every line
/// reaches the judge immediately; `print` never flushes (finish the line with
/// `println`, or call `flush()`).
struct Out {
    w: BufWriter<StdoutLock<'static>>,
    interactive: bool,
}

impl Out {
    fn new(interactive: bool) -> Self {
        Out { w: BufWriter::with_capacity(1 << 17, io::stdout().lock()), interactive }
    }

    /// `print!` analog: writes `x` without a trailing newline.
    fn print<T: Display>(&mut self, x: T) {
        write!(self.w, "{}", x).unwrap();
    }

    /// `println!` analog: writes `x` followed by a newline.
    fn println<T: Display>(&mut self, x: T) {
        writeln!(self.w, "{}", x).unwrap();
        self.end_line();
    }

    /// Writes all items of `v` joined by `sep`, followed by a newline
    /// (an empty slice prints just the newline).
    fn print_vec<T: Display>(&mut self, v: &[T], sep: &str) {
        let mut it = v.iter();
        if let Some(first) = it.next() {
            write!(self.w, "{}", first).unwrap();
            for x in it {
                write!(self.w, "{}{}", sep, x).unwrap();
            }
        }
        writeln!(self.w).unwrap();
        self.end_line();
    }

    fn end_line(&mut self) {
        if self.interactive {
            self.w.flush().unwrap();
        }
    }

    fn flush(&mut self) {
        self.w.flush().unwrap();
    }
}

impl Drop for Out {
    fn drop(&mut self) {
        let _ = self.w.flush();
    }
}


pub struct Math<M: Modulus> {
    fact: Vec<StaticModInt<M>>,
    inv_fact: Vec<StaticModInt<M>>,
}

impl<M: Modulus> Math<M> {
    pub fn new(max_n: usize) -> Self {
        let mut fact = vec![StaticModInt::<M>::new(1); max_n + 1];
        let mut inv_fact = vec![StaticModInt::<M>::new(1); max_n + 1];

        for i in 1..=max_n {
            fact[i] = fact[i - 1] * StaticModInt::<M>::new(i);
        }

        inv_fact[max_n] = fact[max_n].inv();
        for i in (1..=max_n).rev() {
            inv_fact[i - 1] = inv_fact[i] * StaticModInt::<M>::new(i);
        }

        Self { fact, inv_fact }
    }

    pub fn ncr(&self, n: usize, r: usize) -> StaticModInt<M> {
        if r > n {
            return StaticModInt::<M>::new(0);
        }
        self.fact[n] * self.inv_fact[r] * self.inv_fact[n - r]
    }

    pub fn fact(&self, n: usize) -> StaticModInt<M> {
        self.fact[n]
    }

    pub fn inv_fact(&self, n: usize) -> StaticModInt<M> {
        self.inv_fact[n]
    }
}


#[derive(Debug)]
struct Node {
    key: i64,
    priority: u32,
    left: usize,
    right: usize,
    size: usize,
}

impl Node {
    fn new(key: i64, seed: &mut u32) -> Self {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 17;
        *seed ^= *seed << 5;
        Self {
            key,
            priority: *seed,
            left: 0,
            right: 0,
            size: 1,
        }
    }

    fn update(idx: usize, nodes: &mut [Node]) {
        if idx == 0 {
            return;
        }
        let l = nodes[idx - 1].left;
        let r = nodes[idx - 1].right;
        let l_size = if l == 0 { 0 } else { nodes[l - 1].size };
        let r_size = if r == 0 { 0 } else { nodes[r - 1].size };
        nodes[idx - 1].size = 1 + l_size + r_size;
    }
}

pub struct SortedList {
    nodes: Vec<Node>,
    root: usize,
    seed: u32,
    free_stack: Vec<usize>,
}

impl SortedList {
    pub fn new() -> Self {
        Self {
            nodes: Vec::with_capacity(1024 * 8),
            root: 0,
            seed: 12345,
            free_stack: Vec::new(),
        }
    }

    fn split(&mut self, node_idx: usize, key: i64, inclusive: bool) -> (usize, usize) {
        if node_idx == 0 {
            return (0, 0);
        }

        let n_key = self.nodes[node_idx - 1].key;
        let go_left = if inclusive { n_key >= key } else { n_key > key };

        if go_left {
            let left_child = self.nodes[node_idx - 1].left;
            let (l, r) = self.split(left_child, key, inclusive);
            self.nodes[node_idx - 1].left = r;
            Node::update(node_idx, &mut self.nodes);
            (l, node_idx)
        } else {
            let right_child = self.nodes[node_idx - 1].right;
            let (l, r) = self.split(right_child, key, inclusive);
            self.nodes[node_idx - 1].right = l;
            Node::update(node_idx, &mut self.nodes);
            (node_idx, r)
        }
    }

    fn merge(&mut self, l: usize, r: usize) -> usize {
        if l == 0 || r == 0 {
            return if l == 0 { r } else { l };
        }

        if self.nodes[l - 1].priority > self.nodes[r - 1].priority {
            let l_right = self.nodes[l - 1].right;
            self.nodes[l - 1].right = self.merge(l_right, r);
            Node::update(l, &mut self.nodes);
            l
        } else {
            let r_left = self.nodes[r - 1].left;
            self.nodes[r - 1].left = self.merge(l, r_left);
            Node::update(r, &mut self.nodes);
            r
        }
    }

    pub fn insert(&mut self, key: i64) {
        let (l, r) = self.split(self.root, key, true);
        let node_idx = if let Some(idx) = self.free_stack.pop() {
            self.nodes[idx - 1] = Node::new(key, &mut self.seed);
            idx
        } else {
            self.nodes.push(Node::new(key, &mut self.seed));
            self.nodes.len()
        };
        let merged_left = self.merge(l, node_idx);
        self.root = self.merge(merged_left, r);
    }


    pub fn remove(&mut self, key: i64) {
        let (l, mid_r) = self.split(self.root, key, true);
        let (mid, r) = self.split(mid_r, key + 1, true);
        if mid != 0 {
            let m_left = self.nodes[mid - 1].left;
            let m_right = self.nodes[mid - 1].right;
            let new_mid = self.merge(m_left, m_right);
            self.free_stack.push(mid);
            let merged_left = self.merge(l, new_mid);
            self.root = self.merge(merged_left, r);
        } else {
            self.root = self.merge(l, r);
        }
    }


    pub fn bisect_left(&self, key: i64) -> usize {
        let mut curr = self.root;
        let mut rank = 0;
        while curr != 0 {
            let node = &self.nodes[curr - 1];
            if node.key >= key {
                curr = node.left;
            } else {
                let l_size = if node.left == 0 { 0 } else { self.nodes[node.left - 1].size };
                rank += 1 + l_size;
                curr = node.right;
            }
        }
        rank
    }


    pub fn bisect_right(&self, key: i64) -> usize {
        let mut curr = self.root;
        let mut rank = 0;
        while curr != 0 {
            let node = &self.nodes[curr - 1];
            if node.key > key {
                curr = node.left;
            } else {
                let l_size = if node.left == 0 { 0 } else { self.nodes[node.left - 1].size };
                rank += 1 + l_size;
                curr = node.right;
            }
        }
        rank
    }

    pub fn get(&self, index: usize) -> Option<i64> {
        let mut curr = self.root;
        let mut i = index;
        while curr != 0 {
            let node = &self.nodes[curr - 1];
            let left_size = if node.left == 0 { 0 } else { self.nodes[node.left - 1].size };
            match i.cmp(&left_size) {
                Ordering::Less => curr = node.left,
                Ordering::Equal => return Some(node.key),
                Ordering::Greater => {
                    i -= left_size + 1;
                    curr = node.right;
                }
            }
        }
        None
    }


    pub fn len(&self) -> usize {
        if self.root == 0 {
            0
        } else {
            self.nodes[self.root - 1].size
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
