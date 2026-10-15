//! Union-Find with potentials over an arbitrary [`Group`] (weighted union-find).
//!
//! Every vertex `v` has an unknown group element `a_v`. The structure stores relations of the
//! form `a_u = a_v ⊗ x` and answers whether two vertices are related and, if so, the relation
//! `r` with `a_u = a_v ⊗ r`. The group may be non-commutative (matrices, permutations,
//! affine maps); for a commutative one such as `(Z_p, +)` this is the classic "weighted
//! union-find".
//!
//! Internally `diff[v]` is the element with `a_v = a_parent(v) ⊗ diff[v]`; path compression
//! multiplies the differences along the path in order, and union by size keeps paths
//! `O(log n)`, so every operation is `O(log n)` and recursion depth is at most `log2 n`.
//!
//! # Example
//!
//! ```
//! use algebra::Group;
//! use union_find::PotentialUnionFind;
//!
//! /// Integers under addition.
//! enum Add {}
//! impl Monoid for Add {
//!     type S = i64;
//!     fn identity() -> i64 { 0 }
//!     fn op(a: &i64, b: &i64) -> i64 { a + b }
//! }
//! impl Group for Add {
//!     fn inv(x: &i64) -> i64 { -x }
//! }
//!
//! let mut uf = PotentialUnionFind::<Add>::new(4);
//! assert!(uf.unite(0, 1, 5)); // a_0 = a_1 + 5
//! assert!(uf.unite(1, 2, 2)); // a_1 = a_2 + 2
//! assert_eq!(uf.relation(0, 2), 7); // a_0 = a_2 + 7
//! assert!(!uf.unite(0, 2, 8)); // contradicts the above
//! assert!(!uf.same(0, 3));
//! ```

use algebra::Group;

/// Union-find with group-valued potentials.
pub struct PotentialUnionFind<G: Group> {
    /// `< 0`: root with size `-value`; otherwise the parent.
    parent_or_size: Vec<i32>,
    /// `a_v = a_parent(v) ⊗ diff[v]`; the identity for roots.
    diff: Vec<G::S>,
}

impl<G: Group> PotentialUnionFind<G>
where
    G::S: PartialEq,
{
    /// `n` vertices with no relations.
    pub fn new(n: usize) -> Self {
        assert!(n < i32::MAX as usize, "too many elements: {n}");
        Self { parent_or_size: vec![-1; n], diff: vec![G::identity(); n] }
    }

    /// Root of `x`; afterwards `diff[x]` is relative to that root.
    fn find(&mut self, x: usize) -> usize {
        if self.parent_or_size[x] < 0 {
            return x;
        }
        let p = self.parent_or_size[x] as usize;
        let root = self.find(p);
        if p != root {
            // a_x = a_p ⊗ d_x and a_p = a_root ⊗ d_p  =>  a_x = a_root ⊗ (d_p ⊗ d_x)
            self.diff[x] = G::op(&self.diff[p], &self.diff[x]);
            self.parent_or_size[x] = root as i32;
        }
        root
    }

    /// Whether `u` and `v` are related.
    pub fn same(&mut self, u: usize, v: usize) -> bool {
        self.find(u) == self.find(v)
    }

    /// Size of the component of `x`.
    pub fn size(&mut self, x: usize) -> usize {
        let r = self.find(x);
        (-self.parent_or_size[r]) as usize
    }

    /// The element `r` with `a_u = a_v ⊗ r`.
    ///
    /// # Panics
    /// Panics if `u` and `v` are not related (check with [`Self::same`]).
    pub fn relation(&mut self, u: usize, v: usize) -> G::S {
        assert!(self.same(u, v), "vertices {u} and {v} are not related");
        // a_u = a_r ⊗ d_u, a_v = a_r ⊗ d_v  =>  a_u = a_v ⊗ (d_v^-1 ⊗ d_u)
        G::op(&G::inv(&self.diff[v]), &self.diff[u])
    }

    /// Record `a_u = a_v ⊗ x`. Returns `true` if the information is consistent with everything
    /// recorded so far (it is then stored), `false` if it contradicts it (nothing changes).
    pub fn unite(&mut self, u: usize, v: usize, x: G::S) -> bool {
        let (ru, rv) = (self.find(u), self.find(v));
        if ru == rv {
            return self.diff[u] == G::op(&self.diff[v], &x);
        }
        // a_ru = a_rv ⊗ (d_v ⊗ x ⊗ d_u^-1)
        let ru_in_rv = G::op(&G::op(&self.diff[v], &x), &G::inv(&self.diff[u]));
        if self.parent_or_size[ru] <= self.parent_or_size[rv] {
            // ru is the larger component: hang rv below it with the inverse relation.
            self.parent_or_size[ru] += self.parent_or_size[rv];
            self.parent_or_size[rv] = ru as i32;
            self.diff[rv] = G::inv(&ru_in_rv);
        } else {
            self.parent_or_size[rv] += self.parent_or_size[ru];
            self.parent_or_size[ru] = rv as i32;
            self.diff[ru] = ru_in_rv;
        }
        true
    }
}
