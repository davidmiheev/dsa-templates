//! Link-cut tree over a forest whose vertices carry values of a [`Monoid`].
//!
//! * [`link`](LinkCutTree::link) / [`cut`](LinkCutTree::cut) edges in `O(log n)` amortized,
//! * [`connected`](LinkCutTree::connected), [`lca`](LinkCutTree::lca),
//!   [`find_root`](LinkCutTree::find_root),
//! * [`set`](LinkCutTree::set) / [`get`](LinkCutTree::get) a vertex value,
//! * [`path_prod`](LinkCutTree::path_prod)`(u, v)`: the fold of the values on the path from
//!   `u` to `v`, **in path order** — the monoid need not be commutative.
//!
//! Operations that need a precondition (`link` of two different trees, `cut` of an existing
//! edge, `path_prod` between connected vertices) assert it instead of returning an `Option`.
//!
//! # Example
//!
//! ```
//! use algebra::Monoid;
//! use link_cut_tree::LinkCutTree;
//!
//! enum Concat {}
//! impl Monoid for Concat {
//!     type S = String;
//!     fn identity() -> String { String::new() }
//!     fn op(a: &String, b: &String) -> String { format!("{a}{b}") }
//! }
//!
//! let mut t = LinkCutTree::<Concat>::new(&["a", "b", "c", "d"].map(String::from));
//! t.link(0, 1);
//! t.link(1, 2);
//! t.link(1, 3);
//! assert_eq!(t.path_prod(0, 2), "abc");
//! assert_eq!(t.path_prod(2, 3), "cbd");
//! t.cut(1, 2);
//! assert!(!t.connected(0, 2));
//! assert_eq!(t.path_prod(3, 0), "dba");
//! ```

use algebra::Monoid;

const NIL: u32 = 0;

struct Node<S> {
    child: [u32; 2],
    parent: u32,
    reversed: bool,
    value: S,
    /// Fold of the splay subtree in left-to-right order, and in right-to-left order.
    forward: S,
    backward: S,
}

/// See the [module documentation](self).
pub struct LinkCutTree<M: Monoid> {
    /// Node `v + 1` is vertex `v`; node 0 is the empty tree.
    nodes: Vec<Node<M::S>>,
}

impl<M: Monoid> LinkCutTree<M> {
    /// A forest of `values.len()` isolated vertices carrying `values`.
    pub fn new(values: &[M::S]) -> Self {
        let mut nodes = Vec::with_capacity(values.len() + 1);
        let id = M::identity();
        nodes.push(Node { child: [NIL; 2], parent: NIL, reversed: false, value: id.clone(), forward: id.clone(), backward: id });
        for v in values {
            nodes.push(Node { child: [NIL; 2], parent: NIL, reversed: false, value: v.clone(), forward: v.clone(), backward: v.clone() });
        }
        LinkCutTree { nodes }
    }

    /// Number of vertices.
    pub fn len(&self) -> usize {
        self.nodes.len() - 1
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn check(&self, v: usize) -> u32 {
        assert!(v < self.len(), "vertex {v} out of range for {} vertices", self.len());
        v as u32 + 1
    }

    fn is_root(&self, x: u32) -> bool {
        let p = self.nodes[x as usize].parent;
        p == NIL || (self.nodes[p as usize].child[0] != x && self.nodes[p as usize].child[1] != x)
    }

    fn pull(&mut self, x: u32) {
        let [l, r] = self.nodes[x as usize].child;
        let forward = M::op(&M::op(&self.nodes[l as usize].forward, &self.nodes[x as usize].value), &self.nodes[r as usize].forward);
        let backward = M::op(&M::op(&self.nodes[r as usize].backward, &self.nodes[x as usize].value), &self.nodes[l as usize].backward);
        let node = &mut self.nodes[x as usize];
        node.forward = forward;
        node.backward = backward;
    }

    fn flip(&mut self, x: u32) {
        if x != NIL {
            let node = &mut self.nodes[x as usize];
            node.child.swap(0, 1);
            std::mem::swap(&mut node.forward, &mut node.backward);
            node.reversed ^= true;
        }
    }

    fn push(&mut self, x: u32) {
        if self.nodes[x as usize].reversed {
            let [l, r] = self.nodes[x as usize].child;
            self.flip(l);
            self.flip(r);
            self.nodes[x as usize].reversed = false;
        }
    }

    fn rotate(&mut self, x: u32) {
        let p = self.nodes[x as usize].parent;
        let g = self.nodes[p as usize].parent;
        let dir = (self.nodes[p as usize].child[1] == x) as usize;
        if !self.is_root(p) {
            let gd = (self.nodes[g as usize].child[1] == p) as usize;
            self.nodes[g as usize].child[gd] = x;
        }
        self.nodes[x as usize].parent = g;
        let inner = self.nodes[x as usize].child[dir ^ 1];
        self.nodes[p as usize].child[dir] = inner;
        if inner != NIL {
            self.nodes[inner as usize].parent = p;
        }
        self.nodes[x as usize].child[dir ^ 1] = p;
        self.nodes[p as usize].parent = x;
        self.pull(p);
    }

    fn splay(&mut self, x: u32) {
        // Push reversals from the splay root down to x.
        let mut stack = vec![x];
        let mut y = x;
        while !self.is_root(y) {
            y = self.nodes[y as usize].parent;
            stack.push(y);
        }
        while let Some(y) = stack.pop() {
            self.push(y);
        }
        while !self.is_root(x) {
            let p = self.nodes[x as usize].parent;
            if !self.is_root(p) {
                let g = self.nodes[p as usize].parent;
                let same = (self.nodes[g as usize].child[0] == p) == (self.nodes[p as usize].child[0] == x);
                self.rotate(if same { p } else { x });
            }
            self.rotate(x);
        }
        self.pull(x);
    }

    /// Make the root-to-`x` path preferred and splay `x` to the top of it. Returns the last
    /// node whose preferred path was changed (for `access(u); access(v)` that is the LCA).
    fn access(&mut self, x: u32) -> u32 {
        let mut last = NIL;
        let mut y = x;
        let mut top = x;
        while y != NIL {
            self.splay(y);
            self.nodes[y as usize].child[1] = last;
            self.pull(y);
            top = y;
            last = y;
            y = self.nodes[y as usize].parent;
        }
        self.splay(x);
        top
    }

    fn make_root(&mut self, x: u32) {
        self.access(x);
        self.flip(x);
    }

    fn find_root_node(&mut self, x: u32) -> u32 {
        self.access(x);
        let mut y = x;
        loop {
            self.push(y);
            let l = self.nodes[y as usize].child[0];
            if l == NIL {
                break;
            }
            y = l;
        }
        self.splay(y);
        y
    }

    /// The root of the tree containing `v` in its current orientation (the orientation changes
    /// with [`link`](Self::link), [`cut`](Self::cut) and [`path_prod`](Self::path_prod)); two
    /// vertices are in the same tree exactly when this agrees.
    pub fn find_root(&mut self, v: usize) -> usize {
        let x = self.check(v);
        (self.find_root_node(x) - 1) as usize
    }

    /// Whether `u` and `v` are in the same tree.
    pub fn connected(&mut self, u: usize, v: usize) -> bool {
        self.find_root(u) == self.find_root(v)
    }

    /// Add the edge `(u, v)`.
    ///
    /// # Panics
    /// Panics if `u` and `v` are already connected.
    pub fn link(&mut self, u: usize, v: usize) {
        assert!(!self.connected(u, v), "link({u}, {v}): the vertices are already connected");
        let (x, y) = (self.check(u), self.check(v));
        self.make_root(x);
        self.nodes[x as usize].parent = y;
    }

    /// Remove the edge `(u, v)`.
    ///
    /// # Panics
    /// Panics if there is no such edge.
    pub fn cut(&mut self, u: usize, v: usize) {
        let (x, y) = (self.check(u), self.check(v));
        self.make_root(x);
        self.access(y);
        // The path x..y must consist of exactly these two vertices.
        self.push(y);
        let l = self.nodes[y as usize].child[0];
        assert!(l == x && x != y, "cut({u}, {v}): no such edge");
        self.push(x);
        assert!(self.nodes[x as usize].child == [NIL; 2], "cut({u}, {v}): no such edge");
        self.nodes[y as usize].child[0] = NIL;
        self.nodes[x as usize].parent = NIL;
        self.pull(y);
    }

    /// The lowest common ancestor of `u` and `v` when the tree is rooted at `root`.
    ///
    /// # Panics
    /// Panics unless `u`, `v` and `root` are in one tree.
    pub fn lca(&mut self, root: usize, u: usize, v: usize) -> usize {
        assert!(self.connected(root, u) && self.connected(root, v), "lca: the vertices are not in one tree");
        let r = self.check(root);
        self.make_root(r);
        let (x, y) = (self.check(u), self.check(v));
        self.access(x);
        (self.access(y) - 1) as usize
    }

    /// Replace the value of `v`.
    pub fn set(&mut self, v: usize, value: M::S) {
        let x = self.check(v);
        self.access(x);
        self.nodes[x as usize].value = value;
        self.pull(x);
    }

    /// The value of `v`.
    pub fn get(&mut self, v: usize) -> M::S {
        let x = self.check(v);
        self.access(x);
        self.nodes[x as usize].value.clone()
    }

    /// The fold of the values on the path from `u` to `v`, in path order (`u` first).
    ///
    /// # Panics
    /// Panics if `u` and `v` are not connected.
    pub fn path_prod(&mut self, u: usize, v: usize) -> M::S {
        assert!(self.connected(u, v), "path_prod({u}, {v}): the vertices are not connected");
        let (x, y) = (self.check(u), self.check(v));
        self.make_root(x);
        self.access(y);
        self.nodes[y as usize].forward.clone()
    }
}
