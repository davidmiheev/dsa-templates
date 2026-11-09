//! The double-ended palindromic tree.

use std::collections::VecDeque;

const NONE: u32 = u32::MAX;
const ODD: u32 = 0;
const EVEN: u32 = 1;

/// A palindromic tree of a string that supports `push_front`, `push_back`, `pop_front` and
/// `pop_back`, each in `O(log n)` worst case per step through quick links (series links).
///
/// The structure follows the surface-occurrence technique of Rubinchik and Shur: every position
/// stores the *surface* palindrome that starts (`pre_surf`) and ends (`suf_surf`) there, i.e. an
/// occurrence that is not a proper prefix / suffix of another recorded occurrence. A vertex is
/// removed when it has no surface occurrence left and is the suffix link of no other vertex.
pub struct DequeEertree {
    sigma: usize,
    deq: VecDeque<Slot>,
    // vertex arrays
    len: Vec<i32>,
    parent: Vec<u32>,
    link: Vec<u32>,
    quick: Vec<u32>,
    next: Vec<u32>,
    cnt: Vec<u32>,
    link_cnt: Vec<u32>,
    free: Vec<u32>,
    vertices: usize,
}

#[derive(Clone, Copy)]
struct Slot {
    ch: u8,
    pre_surf: u32,
    suf_surf: u32,
}

impl DequeEertree {
    pub fn new(sigma: usize) -> Self {
        let mut t = DequeEertree {
            sigma,
            deq: VecDeque::new(),
            len: Vec::new(),
            parent: Vec::new(),
            link: Vec::new(),
            quick: Vec::new(),
            next: Vec::new(),
            cnt: Vec::new(),
            link_cnt: Vec::new(),
            free: Vec::new(),
            vertices: 0,
        };
        for len in [-1, 0] {
            let id = t.alloc();
            t.len[id as usize] = len;
            t.parent[id as usize] = ODD;
            t.link[id as usize] = ODD;
            t.quick[id as usize] = ODD;
        }
        t
    }

    fn alloc(&mut self) -> u32 {
        if let Some(id) = self.free.pop() {
            let i = id as usize;
            self.next[i * self.sigma..(i + 1) * self.sigma].fill(NONE);
            self.cnt[i] = 0;
            self.link_cnt[i] = 0;
            return id;
        }
        self.len.push(0);
        self.parent.push(0);
        self.link.push(0);
        self.quick.push(0);
        self.cnt.push(0);
        self.link_cnt.push(0);
        self.next.extend(std::iter::repeat_n(NONE, self.sigma));
        (self.len.len() - 1) as u32
    }

    fn child(&self, v: u32, c: u8) -> u32 {
        self.next[v as usize * self.sigma + c as usize]
    }

    /// Character `i` of the current string.
    fn ch(&self, i: usize) -> u8 {
        self.deq[i].ch
    }

    /// The first palindrome on the link chain of `p` (starting with `p`) that can be extended by
    /// `c` at the back of a string of length `n`.
    fn back_appendable(&self, c: u8, mut p: u32) -> u32 {
        let n = self.deq.len();
        loop {
            let lp = self.len[p as usize];
            if lp == -1 || ((lp as usize) < n && self.ch(n - lp as usize - 1) == c) {
                return p;
            }
            let q = self.link[p as usize];
            let lq = self.len[q as usize];
            if lq == -1 || self.ch(n - lq as usize - 1) == c {
                return q;
            }
            p = self.quick[p as usize];
        }
    }

    fn front_appendable(&self, c: u8, mut p: u32) -> u32 {
        let n = self.deq.len();
        loop {
            let lp = self.len[p as usize];
            if lp == -1 || ((lp as usize) < n && self.ch(lp as usize) == c) {
                return p;
            }
            let q = self.link[p as usize];
            let lq = self.len[q as usize];
            if lq == -1 || self.ch(lq as usize) == c {
                return q;
            }
            p = self.quick[p as usize];
        }
    }

    pub fn push_back(&mut self, c: u8) {
        assert!((c as usize) < self.sigma, "symbol out of range");
        let par = if self.deq.is_empty() { ODD } else { self.back_appendable(c, self.deq.back().expect("non-empty").suf_surf) };
        let slot = Slot { ch: c, pre_surf: EVEN, suf_surf: EVEN };
        let existing = self.child(par, c);
        let (v, w);
        if existing == NONE {
            let id = self.alloc();
            let iv = id as usize;
            self.parent[iv] = par;
            self.len[iv] = self.len[par as usize] + 2;
            let mut ww = EVEN;
            if par != ODD {
                let x = self.back_appendable(c, self.link[par as usize]);
                ww = self.child(x, c);
            }
            self.link[iv] = ww;
            self.link_cnt[ww as usize] += 1;
            self.deq.push_back(slot);
            let n = self.deq.len();
            let lw = self.len[ww as usize] as usize;
            let lw_link = self.link[ww as usize];
            self.quick[iv] = if lw_link != ODD && self.ch(n - lw - 1) == self.ch(n - self.len[lw_link as usize] as usize - 1) { self.quick[ww as usize] } else { lw_link };
            self.next[par as usize * self.sigma + c as usize] = id;
            self.vertices += 1;
            v = id;
            w = ww;
        } else {
            self.deq.push_back(slot);
            v = existing;
            w = self.link[v as usize];
        }
        let n = self.deq.len();
        let lv = self.len[v as usize] as usize;
        let lw = self.len[w as usize];
        self.deq[n - 1].suf_surf = v;
        self.deq[n - lv].pre_surf = v;
        if lw >= 1 {
            let idx = n - lv + lw as usize - 1;
            if self.deq[idx].suf_surf == w {
                self.deq[idx].suf_surf = EVEN;
            }
        }
        self.cnt[v as usize] += 1;
    }

    pub fn push_front(&mut self, c: u8) {
        assert!((c as usize) < self.sigma, "symbol out of range");
        let par = if self.deq.is_empty() { ODD } else { self.front_appendable(c, self.deq.front().expect("non-empty").pre_surf) };
        let slot = Slot { ch: c, pre_surf: EVEN, suf_surf: EVEN };
        let existing = self.child(par, c);
        let (v, w);
        if existing == NONE {
            let id = self.alloc();
            let iv = id as usize;
            self.parent[iv] = par;
            self.len[iv] = self.len[par as usize] + 2;
            let mut ww = EVEN;
            if par != ODD {
                let x = self.front_appendable(c, self.link[par as usize]);
                ww = self.child(x, c);
            }
            self.link[iv] = ww;
            self.link_cnt[ww as usize] += 1;
            self.deq.push_front(slot);
            let lw = self.len[ww as usize] as usize;
            let lw_link = self.link[ww as usize];
            self.quick[iv] = if lw_link != ODD && self.ch(lw) == self.ch(self.len[lw_link as usize] as usize) { self.quick[ww as usize] } else { lw_link };
            self.next[par as usize * self.sigma + c as usize] = id;
            self.vertices += 1;
            v = id;
            w = ww;
        } else {
            self.deq.push_front(slot);
            v = existing;
            w = self.link[v as usize];
        }
        let lv = self.len[v as usize] as usize;
        let lw = self.len[w as usize];
        self.deq[0].pre_surf = v;
        self.deq[lv - 1].suf_surf = v;
        if lw >= 1 {
            let idx = lv - lw as usize;
            if self.deq[idx].pre_surf == w {
                self.deq[idx].pre_surf = EVEN;
            }
        }
        self.cnt[v as usize] += 1;
    }

    fn release(&mut self, v: u32, c: u8) {
        self.cnt[v as usize] -= 1;
        if self.link_cnt[v as usize] == 0 && self.cnt[v as usize] == 0 {
            let par = self.parent[v as usize];
            self.next[par as usize * self.sigma + c as usize] = NONE;
            let w = self.link[v as usize];
            self.link_cnt[w as usize] -= 1;
            self.free.push(v);
            self.vertices -= 1;
        }
    }

    pub fn pop_back(&mut self) {
        let last = *self.deq.back().expect("the string is empty");
        let (v, c) = (last.suf_surf, last.ch);
        let w = self.link[v as usize];
        let n = self.deq.len();
        let (lv, lw) = (self.len[v as usize] as usize, self.len[w as usize] as usize);
        if lv >= 2 && lw > 0 && self.len[self.deq[n - lv + lw - 1].suf_surf as usize] < lw as i32 {
            self.deq[n - lv + lw - 1].suf_surf = w;
            self.deq[n - lv].pre_surf = w;
        } else {
            self.deq[n - lv].pre_surf = EVEN;
        }
        self.release(v, c);
        self.deq.pop_back();
    }

    pub fn pop_front(&mut self) {
        let first = *self.deq.front().expect("the string is empty");
        let (v, c) = (first.pre_surf, first.ch);
        let w = self.link[v as usize];
        let (lv, lw) = (self.len[v as usize] as usize, self.len[w as usize] as usize);
        if lv >= 2 && lw > 0 && self.len[self.deq[lv - lw].pre_surf as usize] < lw as i32 {
            self.deq[lv - lw].pre_surf = w;
            self.deq[lv - 1].suf_surf = w;
        } else {
            self.deq[lv - 1].suf_surf = EVEN;
        }
        self.release(v, c);
        self.deq.pop_front();
    }

    /// The length of the current string.
    pub fn len(&self) -> usize {
        self.deq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.deq.is_empty()
    }

    /// The number of distinct non-empty palindromic substrings.
    pub fn distinct_palindromes(&self) -> usize {
        self.vertices
    }

    /// The length of the longest palindromic prefix (0 for the empty string).
    pub fn longest_prefix_palindrome(&self) -> usize {
        self.deq.front().map_or(0, |s| self.len[s.pre_surf as usize] as usize)
    }

    /// The length of the longest palindromic suffix (0 for the empty string).
    pub fn longest_suffix_palindrome(&self) -> usize {
        self.deq.back().map_or(0, |s| self.len[s.suf_surf as usize] as usize)
    }
}
