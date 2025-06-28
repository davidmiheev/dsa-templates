Union Find (Disjoint Set)
=========================

Union Find is a data structure that keeps track of elements which are partitioned into disjoint sets.
It provides operations to add new sets, merge existing sets, and find the set to which a particular element belongs.

Rust Implementation
-------------------

.. code-block:: rust

    struct UnionFind<F> 
    where F: Fn(usize, usize) -> usize 
    {
        root: Vec<usize>,
        size: Vec<usize>,
        components: usize, 
        character: Vec<usize>,
        op: F,
    }
    
    impl<F> UnionFind<F> 
    where F: Fn(usize, usize) -> usize
    {
        pub fn new(n: usize, op: F, init: usize) -> Self {
            Self {
                root: (0..=n).collect(),
                size: vec![1; n + 1],
                components: n,
                character: vec![init; n + 1],
                op,
            }
        }
    
        pub fn find(&mut self, x: usize) -> usize {
            if self.root[x] == x { return x }
            self.root[x] = self.find(self.root[x]);
            self.root[x]
        }
    
        pub fn union(&mut self, x: usize, y: usize, w: usize) -> bool {
            let root_x = self.find(x);
            let root_y = self.find(y);
    
            if root_x == root_y {
                self.character[root_x] = (self.op)(self.character[root_x], w);
                return false;
            }
    
            if self.size[root_x] > self.size[root_y] {
                self.size[root_x] += self.size[root_y];
                self.character[root_x] = (self.op)((self.op)(self.character[root_x], self.character[root_y]), w);
                self.root[root_y] = root_x;
            } else {
                self.size[root_y] += self.size[root_x];
                self.character[root_y] = (self.op)((self.op)(self.character[root_x], self.character[root_y]), w);
                self.root[root_x] = root_y;
            }
            self.components -= 1;
            true
        }

        pub fn are_connected(&mut self, x: usize, y: usize) -> bool {
            self.find(x) == self.find(y)
        }

        pub fn get_components(&self) -> usize {
            self.components
        }

        pub fn get_root(&mut self, x: usize) -> usize {
            self.find(x)
        }

        pub fn get_character(&mut self, x: usize) -> usize {
            let root_x = self.find(x);
            self.character[root_x]
        }

        pub fn get_size(&mut self, x: usize) -> usize {
            let root_x = self.find(x);
            self.size[root_x]
        }
    }



Python Implementation
---------------------
.. automodule:: union_find.union_find
    :members:
    :undoc-members:
    :show-inheritance:
    :inherited-members:


