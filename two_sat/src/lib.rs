//! 2-SAT in `O(n + m)`: the implication graph is decomposed into strongly connected components
//! ([`scc`]); the formula is satisfiable iff no variable shares a component with its negation, and
//! a variable is true iff its positive literal is later in the topological order than its negation.
//!
//! # Example
//!
//! ```
//! use two_sat::TwoSat;
//!
//! let mut sat = TwoSat::new(2);
//! sat.add_clause(0, true, 1, true); // x0 or x1
//! sat.add_clause(0, false, 1, false); // not x0 or not x1
//! let (ok, x) = sat.solve();
//! assert!(ok && x[0] != x[1]);
//! sat.add_clause(0, true, 0, true); // x0
//! sat.add_clause(1, true, 1, true); // x1
//! assert!(!sat.solve().0);
//! ```

use scc::Scc;

/// A 2-SAT instance over variables `0..n`.
pub struct TwoSat {
    n: usize,
    arcs: Vec<(usize, usize)>,
}

impl TwoSat {
    pub fn new(n: usize) -> Self {
        TwoSat { n, arcs: Vec::new() }
    }

    /// The clause `(x_i == f) or (x_j == g)`.
    pub fn add_clause(&mut self, i: usize, f: bool, j: usize, g: bool) {
        assert!(i < self.n && j < self.n, "variable out of range");
        // literal of x_i == f has id 2 i + f; its negation is id ^ 1
        let (a, b) = (2 * i + f as usize, 2 * j + g as usize);
        self.arcs.push((a ^ 1, b));
        self.arcs.push((b ^ 1, a));
    }

    /// `(satisfiable, assignment)`; the assignment is empty when unsatisfiable.
    pub fn solve(&self) -> (bool, Vec<bool>) {
        let scc = Scc::from_edges(2 * self.n, &self.arcs);
        let mut assignment = Vec::with_capacity(self.n);
        for i in 0..self.n {
            let (neg, pos) = (scc.comp[2 * i], scc.comp[2 * i + 1]);
            if neg == pos {
                return (false, Vec::new());
            }
            assignment.push(pos > neg);
        }
        (true, assignment)
    }
}
