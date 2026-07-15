//! Suffix Automaton (SAM).
//!
//! A suffix automaton is a directed acyclic automaton that recognises
//! exactly the set of substrings of a given string. It has at most
//! `2n - 1` states for an input of length `n`, and every state
//! corresponds to a distinct set of end positions in the original
//! string — hence to a distinct set of substrings.
//!
//! Once built, you can answer:
//!
//! * Is `p` a substring of `s`? (O(|p|))
//! * How many distinct substrings does `s` have? (sum over states of
//!   `len(state) - len(link(state))`)
//! * Number of occurrences of a pattern (when `state.occ` is computed)
//! * Longest common substring between two strings
//!
//! # Example
//!
//! ```
//! use suffix_automaton::SuffixAutomaton;
//!
//! let mut sam = SuffixAutomaton::new();
//! for ch in "ababa".chars() {
//!     sam.extend(ch);
//! }
//! assert!(sam.contains("aba"));
//! assert!(sam.contains("bab"));
//! assert!(!sam.contains("bb"));
//! ```

/// State of a suffix automaton.
#[derive(Debug, Clone)]
pub struct State {
    /// Length of the longest substring represented by this state.
    pub len: usize,
    /// Suffix link: points to the state representing the longest proper
    /// suffix of this state's substrings.
    pub link: Option<usize>,
    /// Transitions keyed by character.
    pub next: std::collections::BTreeMap<char, usize>,
}

/// Suffix Automaton over a sequence of `char`s.
#[derive(Debug, Clone)]
pub struct SuffixAutomaton {
    /// All states, indexed by id. State 0 is the initial state.
    pub states: Vec<State>,
    /// Id of the state representing the entire current string.
    pub last: usize,
}

impl SuffixAutomaton {
    /// Build an empty suffix automaton.
    pub fn new() -> Self {
        SuffixAutomaton {
            states: vec![State {
                len: 0,
                link: None,
                next: std::collections::BTreeMap::new(),
            }],
            last: 0,
        }
    }

    /// Extend the automaton with a single character.
    pub fn extend(&mut self, ch: char) {
        let cur = self.states.len();
        self.states.push(State {
            len: self.states[self.last].len + 1,
            link: None,
            next: std::collections::BTreeMap::new(),
        });

        let mut p = self.last;
        while p != usize::MAX && !self.states[p].next.contains_key(&ch) {
            self.states[p].next.insert(ch, cur);
            p = self.states[p].link.unwrap_or(usize::MAX);
        }

        if p == usize::MAX {
            self.states[cur].link = Some(0);
        } else {
            let q = self.states[p].next[&ch];
            if self.states[p].len + 1 == self.states[q].len {
                self.states[cur].link = Some(q);
            } else {
                // Clone state q into a new state `clone`.
                let clone = self.states.len();
                self.states.push(State {
                    len: self.states[p].len + 1,
                    link: self.states[q].link,
                    next: self.states[q].next.clone(),
                });
                while p != usize::MAX && self.states[p].next.get(&ch) == Some(&q) {
                    self.states[p].next.insert(ch, clone);
                    p = self.states[p].link.unwrap_or(usize::MAX);
                }
                self.states[q].link = Some(clone);
                self.states[cur].link = Some(clone);
            }
        }

        self.last = cur;
    }

    /// Build the automaton from an entire string in one call.
    pub fn from_str(s: &str) -> Self {
        let mut sam = SuffixAutomaton::new();
        for ch in s.chars() {
            sam.extend(ch);
        }
        sam
    }

    /// Returns `true` iff `pattern` is a substring of the built string.
    pub fn contains(&self, pattern: &str) -> bool {
        let mut v = 0usize;
        for ch in pattern.chars() {
            match self.states[v].next.get(&ch) {
                Some(&u) => v = u,
                None => return false,
            }
        }
        true
    }

    /// Number of distinct substrings of the built string.
    pub fn distinct_substring_count(&self) -> u64 {
        // Σ (len(v) - len(link(v))) over all states v != 0.
        let mut total: u64 = 0;
        for (i, s) in self.states.iter().enumerate() {
            if i == 0 {
                continue;
            }
            let link_len = s.link.map(|l| self.states[l].len).unwrap_or(0);
            total += (s.len - link_len) as u64;
        }
        total
    }
}

impl Default for SuffixAutomaton {
    fn default() -> Self {
        Self::new()
    }
}
