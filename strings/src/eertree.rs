//! The palindromic tree (eertree).

const NONE: u32 = u32::MAX;

/// The eertree of a growing string over `0..sigma`: one vertex per distinct non-empty palindromic
/// substring plus the two roots, built online in amortised `O(1)` per character.
///
/// Vertex ids: `0` is the odd root (length -1), `1` the even root (length 0), the vertex created
/// `k`-th is `k + 1`.
pub struct Eertree {
    sigma: usize,
    s: Vec<u8>,
    len: Vec<i32>,
    link: Vec<u32>,
    parent: Vec<u32>,
    next: Vec<u32>,
    /// the longest palindromic suffix after each prefix
    suffix_of_prefix: Vec<u32>,
    last: u32,
}

impl Eertree {
    pub fn new(sigma: usize) -> Self {
        let mut t = Eertree { sigma, s: Vec::new(), len: Vec::new(), link: Vec::new(), parent: Vec::new(), next: Vec::new(), suffix_of_prefix: Vec::new(), last: 1 };
        t.new_vertex(-1, 0, 0);
        t.new_vertex(0, 0, 0);
        t
    }

    fn new_vertex(&mut self, len: i32, link: u32, parent: u32) -> u32 {
        self.len.push(len);
        self.link.push(link);
        self.parent.push(parent);
        self.next.extend(std::iter::repeat_n(NONE, self.sigma));
        (self.len.len() - 1) as u32
    }

    fn can_extend(&self, v: u32, c: u8) -> bool {
        let n = self.s.len() as i32;
        let l = self.len[v as usize];
        l == -1 || (n - l > 0 && self.s[(n - l - 1) as usize] == c)
    }

    /// Appends the symbol `c`.
    pub fn push(&mut self, c: u8) {
        assert!((c as usize) < self.sigma, "symbol out of range");
        let mut cur = self.last;
        while !self.can_extend(cur, c) {
            cur = self.link[cur as usize];
        }
        self.s.push(c);
        let slot = cur as usize * self.sigma + c as usize;
        if self.next[slot] == NONE {
            let new_len = self.len[cur as usize] + 2;
            let link = if new_len == 1 {
                1
            } else {
                let mut p = self.link[cur as usize];
                // the string already contains c at the end, so test against s[..n-1]
                loop {
                    let n = self.s.len() as i32;
                    let l = self.len[p as usize];
                    if n - l - 2 >= 0 && self.s[(n - l - 2) as usize] == c {
                        break;
                    }
                    p = self.link[p as usize];
                }
                self.next[p as usize * self.sigma + c as usize]
            };
            let v = self.new_vertex(new_len, link, cur);
            self.next[slot] = v;
        }
        self.last = self.next[slot];
        self.suffix_of_prefix.push(self.last);
    }

    /// The number of distinct non-empty palindromic substrings so far.
    pub fn distinct_palindromes(&self) -> usize {
        self.len.len() - 2
    }

    /// The length of the palindrome of vertex `v`.
    pub fn length(&self, v: usize) -> usize {
        assert!(v >= 1, "roots have no length");
        self.len[v].max(0) as usize
    }

    /// The vertex obtained by removing the first and last symbol (`0` = odd root for single
    /// symbols, `1` = even root for pairs).
    pub fn parent(&self, v: usize) -> usize {
        self.parent[v] as usize
    }

    /// The vertex of the longest proper palindromic suffix (the even root if there is none).
    pub fn suffix_link(&self, v: usize) -> usize {
        self.link[v] as usize
    }

    /// The vertex of the longest palindromic suffix of the prefix of length `i + 1`.
    pub fn longest_suffix_palindrome(&self, i: usize) -> usize {
        self.suffix_of_prefix[i] as usize
    }
}
