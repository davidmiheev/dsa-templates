//! A trie of patterns with suffix links.

/// The trie of a set of patterns over `0..sigma` with the Aho–Corasick suffix links.
///
/// Vertices are numbered in the order of creation while inserting the patterns one after another,
/// so vertex `v` corresponds to the lexicographically first `(pattern, prefix length)` that
/// spells it.
pub struct AhoCorasick {
    sigma: usize,
    parent: Vec<u32>,
    link: Vec<u32>,
    /// `goto[v * sigma + c]`: the automaton transition (after the links are built)
    goto: Vec<u32>,
    end: Vec<u32>,
}

const NONE: u32 = u32::MAX;

impl AhoCorasick {
    /// Builds the trie of `patterns`, whose symbols must be below `sigma`.
    pub fn new(patterns: &[&[u8]], sigma: usize) -> Self {
        let total: usize = patterns.iter().map(|p| p.len()).sum();
        let mut ac = AhoCorasick {
            sigma,
            parent: Vec::with_capacity(total + 1),
            link: Vec::with_capacity(total + 1),
            goto: Vec::with_capacity((total + 1) * sigma),
            end: Vec::with_capacity(patterns.len()),
        };
        ac.parent.push(0);
        ac.goto.extend(std::iter::repeat_n(NONE, sigma));
        for p in patterns {
            let mut v = 0usize;
            for &c in p.iter() {
                assert!((c as usize) < sigma, "symbol out of range");
                let nxt = ac.goto[v * sigma + c as usize];
                v = if nxt == NONE {
                    let id = ac.parent.len();
                    ac.parent.push(v as u32);
                    ac.goto.extend(std::iter::repeat_n(NONE, sigma));
                    ac.goto[v * sigma + c as usize] = id as u32;
                    id
                } else {
                    nxt as usize
                };
            }
            ac.end.push(v as u32);
        }
        ac.build_links();
        ac
    }

    fn build_links(&mut self) {
        let (n, sigma) = (self.parent.len(), self.sigma);
        self.link = vec![0; n];
        let mut queue = std::collections::VecDeque::new();
        for c in 0..sigma {
            let v = self.goto[c];
            if v == NONE {
                self.goto[c] = 0;
            } else {
                queue.push_back(v as usize);
            }
        }
        while let Some(v) = queue.pop_front() {
            let link_v = self.link[v] as usize;
            for c in 0..sigma {
                let u = self.goto[v * sigma + c];
                if u == NONE {
                    self.goto[v * sigma + c] = self.goto[link_v * sigma + c];
                } else {
                    self.link[u as usize] = self.goto[link_v * sigma + c];
                    queue.push_back(u as usize);
                }
            }
        }
    }

    /// The number of vertices (including the root 0).
    pub fn len(&self) -> usize {
        self.parent.len()
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    /// The parent of vertex `v >= 1`.
    pub fn parent(&self, v: usize) -> usize {
        assert!(v >= 1 && v < self.len(), "no parent");
        self.parent[v] as usize
    }

    /// The suffix link of `v >= 1`: the vertex of the longest proper suffix of `str(v)` that is a
    /// prefix of some pattern (0 for none).
    pub fn suffix_link(&self, v: usize) -> usize {
        assert!(v >= 1 && v < self.len(), "no suffix link");
        self.link[v] as usize
    }

    /// The vertex spelling pattern `i`.
    pub fn pattern_vertex(&self, i: usize) -> usize {
        self.end[i] as usize
    }

    /// The automaton transition from `v` on symbol `c`.
    pub fn next(&self, v: usize, c: u8) -> usize {
        self.goto[v * self.sigma + c as usize] as usize
    }
}
