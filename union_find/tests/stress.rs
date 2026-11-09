//! Stress test: `UnionFind` against a naive model that relabels whole components.
//!
//! The naive model follows the documented semantics: `union(x, y, w)` combines
//! characters as `op(op(char[x's component], char[y's component]), w)` (always in
//! that argument order, since `op` need not be commutative), or
//! `op(char, w)` when `x` and `y` are already connected.

use union_find::UnionFind;

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

type Op = fn(usize, usize) -> usize;

fn add(a: usize, b: usize) -> usize {
    a.wrapping_add(b)
}
fn xor(a: usize, b: usize) -> usize {
    a ^ b
}
fn max(a: usize, b: usize) -> usize {
    a.max(b)
}
/// Order-sensitive on purpose: catches swapped `(x, y)` arguments.
fn poly(a: usize, b: usize) -> usize {
    a.wrapping_mul(31).wrapping_add(b)
}
const OPS: [Op; 4] = [add, xor, max, poly];

struct Naive {
    comp: Vec<usize>,
    chr: Vec<usize>,
    size: Vec<usize>,
    components: usize,
    op: Op,
}

impl Naive {
    fn new(n: usize, op: Op, init: usize) -> Self {
        Naive {
            comp: (0..=n).collect(),
            chr: vec![init; n + 1],
            size: vec![1; n + 1],
            components: n,
            op,
        }
    }

    fn union(&mut self, x: usize, y: usize, w: usize) -> bool {
        let (cx, cy) = (self.comp[x], self.comp[y]);
        if cx == cy {
            self.chr[cx] = (self.op)(self.chr[cx], w);
            return false;
        }
        let combined = (self.op)((self.op)(self.chr[cx], self.chr[cy]), w);
        for c in self.comp.iter_mut() {
            if *c == cy {
                *c = cx;
            }
        }
        self.chr[cx] = combined;
        self.size[cx] += self.size[cy];
        self.components -= 1;
        true
    }
}

#[test]
fn matches_naive_model() {
    let mut st = Stress::new("union_find/naive");
    for case in 0..400 * scale() {
        let mut rng = Rng::new(case as u64);
        let n = rng.range(1, 40) as usize;
        let op = rng.pick(&OPS);
        let init = rng.pick(&[0usize, 1, 7, usize::MAX]);
        let mut uf = UnionFind::new(n, op, init);
        let mut model = Naive::new(n, op, init);
        st.case();

        for step in 0..200 {
            let x = rng.below(n + 1);
            let y = rng.below(n + 1);
            let ctx = format!("case {case} step {step} n={n} x={x} y={y}");
            match rng.below(6) {
                0..=2 => {
                    let w = rng.pick(&[0usize, 1, 5, 1 << 40, usize::MAX]);
                    assert_eq!(uf.union(x, y, w), model.union(x, y, w), "union {ctx}");
                }
                3 => assert_eq!(
                    uf.are_connected(x, y),
                    model.comp[x] == model.comp[y],
                    "are_connected {ctx}"
                ),
                4 => {
                    assert_eq!(uf.get_size(x), model.size[model.comp[x]], "size {ctx}");
                    assert_eq!(uf.get_character(x), model.chr[model.comp[x]], "character {ctx}");
                }
                _ => {
                    let (rx, ry) = (uf.get_root(x), uf.get_root(y));
                    assert_eq!(rx == ry, model.comp[x] == model.comp[y], "root equality {ctx}");
                    assert_eq!(uf.find(rx), rx, "root must be its own representative {ctx}");
                    assert_eq!(uf.get_root(x), rx, "get_root must be stable {ctx}");
                }
            }
            assert_eq!(uf.get_components(), model.components, "components {ctx}");
            st.ops(1);
        }
        // Final full sweep.
        for x in 0..=n {
            assert_eq!(uf.get_size(x), model.size[model.comp[x]], "final size case {case}");
            assert_eq!(uf.get_character(x), model.chr[model.comp[x]], "final char case {case}");
        }
        st.ops(n + 1);
    }
    st.done();
}

/// Independent reference: plain DSU with per-root sums (op = wrapping add, so order is irrelevant).
struct RefDsu {
    p: Vec<usize>,
    size: Vec<usize>,
    sum: Vec<usize>,
}

impl RefDsu {
    fn find(&mut self, mut x: usize) -> usize {
        while self.p[x] != x {
            self.p[x] = self.p[self.p[x]];
            x = self.p[x];
        }
        x
    }
}

#[test]
fn large_random_matches_reference_dsu() {
    let mut st = Stress::new("union_find/large");
    for case in 0..3 * scale() {
        let mut rng = Rng::new(1000 + case as u64);
        let n = 200_000;
        let mut uf = UnionFind::new(n, add, 0);
        let mut r = RefDsu { p: (0..=n).collect(), size: vec![1; n + 1], sum: vec![0; n + 1] };
        let mut merges = 0;
        st.case();
        for _ in 0..300_000 {
            let (x, y, w) = (rng.below(n + 1), rng.below(n + 1), rng.below(1000));
            let (rx, ry) = (r.find(x), r.find(y));
            let expect = rx != ry;
            if expect {
                let (big, small) = if r.size[rx] >= r.size[ry] { (rx, ry) } else { (ry, rx) };
                r.p[small] = big;
                r.size[big] += r.size[small];
                r.sum[big] = r.sum[rx] + r.sum[ry] + w;
                merges += 1;
            } else {
                r.sum[rx] += w;
            }
            assert_eq!(uf.union(x, y, w), expect, "case {case}");
        }
        assert_eq!(uf.get_components(), n - merges);
        for _ in 0..200_000 {
            let (x, y) = (rng.below(n + 1), rng.below(n + 1));
            let (rx, ry) = (r.find(x), r.find(y));
            assert_eq!(uf.are_connected(x, y), rx == ry);
            assert_eq!(uf.get_size(x), r.size[rx]);
            assert_eq!(uf.get_character(x), r.sum[rx]);
        }
        st.ops(500_000);
    }
    st.done();
}

#[test]
fn adversarial_chains_stay_shallow() {
    // Chains and reverse chains are the classic worst cases for naive linking;
    // union by size must keep `find` fast and non-recursive-overflowing.
    let mut st = Stress::new("union_find/chains");
    let n = 1_000_000;
    for order in 0..3 {
        let mut uf = UnionFind::new(n, add, 0);
        st.case();
        for i in 0..n {
            let (x, y) = match order {
                0 => (i, i + 1),
                1 => (n - i, n - i - 1),
                _ => (0, i + 1),
            };
            uf.union(x, y, 1);
        }
        let root = uf.get_root(0);
        assert_eq!(uf.get_size(0), n + 1, "order {order}");
        for x in (0..=n).step_by(997) {
            assert_eq!(uf.get_root(x), root, "order {order}");
        }
        st.ops(n + n / 997);
    }
    st.done();
}
