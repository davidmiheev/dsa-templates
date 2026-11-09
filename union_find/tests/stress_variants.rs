//! Stress tests for the union-find variants: `PotentialUnionFind` (commutative and
//! non-commutative groups), `RollbackUnionFind` and `RangeParallelUnionFind`.

use algebra::{Group, Monoid};
use std::collections::VecDeque;
use union_find::{PotentialUnionFind, RangeParallelUnionFind, RollbackUnionFind};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

// ---- groups ---------------------------------------------------------------------------------
/// Z_7 under addition (commutative).
enum AddMod7 {}
impl Monoid for AddMod7 {
    type S = u8;
    fn identity() -> u8 {
        0
    }
    fn op(a: &u8, b: &u8) -> u8 {
        (a + b) % 7
    }
}
impl Group for AddMod7 {
    fn inv(x: &u8) -> u8 {
        (7 - x) % 7
    }
}

/// S_4: permutations of four points, `(a ⊗ b)[i] = a[b[i]]` (non-commutative).
enum Perm4 {}
impl Monoid for Perm4 {
    type S = [u8; 4];
    fn identity() -> [u8; 4] {
        [0, 1, 2, 3]
    }
    fn op(a: &[u8; 4], b: &[u8; 4]) -> [u8; 4] {
        [a[b[0] as usize], a[b[1] as usize], a[b[2] as usize], a[b[3] as usize]]
    }
}
impl Group for Perm4 {
    fn inv(x: &[u8; 4]) -> [u8; 4] {
        let mut r = [0u8; 4];
        for i in 0..4 {
            r[x[i] as usize] = i as u8;
        }
        r
    }
}

fn random_perm(rng: &mut Rng) -> [u8; 4] {
    let mut p = [0u8, 1, 2, 3];
    for i in (1..4).rev() {
        p.swap(i, rng.below(i + 1));
    }
    p
}

/// Naive model: the accepted relations as an edge list, queried by BFS.
struct Naive<G: Group> {
    n: usize,
    edges: Vec<(usize, usize, G::S)>, // a_u = a_v ⊗ x
}

impl<G: Group> Naive<G>
where
    G::S: PartialEq + Copy,
{
    /// `rel[w]` with `a_w = a_src ⊗ rel[w]` for every vertex reachable from `src`.
    fn bfs(&self, src: usize) -> Vec<Option<G::S>> {
        let mut rel: Vec<Option<G::S>> = vec![None; self.n];
        rel[src] = Some(G::identity());
        let mut q = VecDeque::from([src]);
        while let Some(w) = q.pop_front() {
            let rw = rel[w].unwrap();
            for &(u, v, x) in &self.edges {
                // a_u = a_v ⊗ x  <=>  a_v = a_u ⊗ x^-1
                if v == w && rel[u].is_none() {
                    // a_u = a_v ⊗ x = a_src ⊗ rw ⊗ x
                    rel[u] = Some(G::op(&rw, &x));
                    q.push_back(u);
                }
                if u == w && rel[v].is_none() {
                    rel[v] = Some(G::op(&rw, &G::inv(&x)));
                    q.push_back(v);
                }
            }
        }
        rel
    }

    /// The relation `r` with `a_u = a_v ⊗ r`, if `u` and `v` are connected.
    fn relation(&self, u: usize, v: usize) -> Option<G::S> {
        // a_v = a_u ⊗ rel[v]  =>  a_u = a_v ⊗ rel[v]^-1
        self.bfs(u)[v].map(|r| G::inv(&r))
    }

    fn component_size(&self, u: usize) -> usize {
        self.bfs(u).iter().flatten().count()
    }
}

fn potential_case<G: Group>(rng: &mut Rng, rand_elem: impl Fn(&mut Rng) -> G::S, st: &mut Stress, case: usize)
where
    G::S: PartialEq + Copy + std::fmt::Debug,
{
    let n = rng.range(1, 9) as usize;
    let hidden: Vec<G::S> = (0..n).map(|_| rand_elem(rng)).collect(); // the "true" a_i
    let mut uf = PotentialUnionFind::<G>::new(n);
    let mut model = Naive::<G> { n, edges: vec![] };
    st.case();
    for step in 0..40 {
        let (u, v) = (rng.below(n), rng.below(n));
        let ctx = format!("case {case} step {step} n={n} u={u} v={v}");
        if rng.chance(1, 2) {
            // Half the time give the true relation (a_u = a_v ⊗ x  =>  x = a_v^-1 ⊗ a_u), else a random one.
            let x = if rng.chance(2, 3) { G::op(&G::inv(&hidden[v]), &hidden[u]) } else { rand_elem(rng) };
            let expect = match model.relation(u, v) {
                Some(r) => r == x,
                None => true,
            };
            assert_eq!(uf.unite(u, v, x), expect, "unite {ctx}");
            if expect {
                model.edges.push((u, v, x));
            }
        } else {
            let want = model.relation(u, v);
            assert_eq!(uf.same(u, v), want.is_some(), "same {ctx}");
            if let Some(r) = want {
                assert_eq!(uf.relation(u, v), r, "relation {ctx}");
            }
            assert_eq!(uf.size(u), model.component_size(u), "size {ctx}");
        }
        st.ops(1);
    }
}

#[test]
fn potential_matches_bfs_model() {
    let mut st = Stress::new("union_find/potential");
    for case in 0..600 * scale() {
        let mut rng = Rng::new(case as u64);
        if case % 2 == 0 {
            potential_case::<AddMod7>(&mut rng, |r| r.below(7) as u8, &mut st, case);
        } else {
            potential_case::<Perm4>(&mut rng, random_perm, &mut st, case);
        }
    }
    st.done();
}

#[test]
fn potential_relation_of_unrelated_vertices_panics() {
    let mut uf = PotentialUnionFind::<AddMod7>::new(3);
    assert!(uf.unite(0, 1, 3));
    assert!(std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| uf.relation(0, 2))).is_err());
}

#[test]
fn potential_long_chains_stay_fast_and_correct() {
    // A path of 200k vertices with relation +1 per edge, built in different orders.
    let mut st = Stress::new("union_find/potential_chain");
    let n = 200_000;
    for order in 0..3 {
        let mut uf = PotentialUnionFind::<AddMod7>::new(n);
        st.case();
        for i in 0..n - 1 {
            let (u, v) = match order {
                0 => (i + 1, i),         // a_{i+1} = a_i + 1
                1 => (n - 1 - i, n - 2 - i), // built from the far end
                _ => ((i * 7919) % (n - 1) + 1, (i * 7919) % (n - 1)), // scattered order
            };
            assert!(uf.unite(u, v, 1));
        }
        assert_eq!(uf.size(0), n);
        for i in [0, 1, 17, n / 2, n - 1] {
            assert_eq!(uf.relation(i, 0), (i % 7) as u8, "order {order} vertex {i}");
        }
        st.ops(n);
    }
    st.done();
}

// ---- rollback -------------------------------------------------------------------------------
#[test]
fn rollback_matches_recomputation() {
    let mut st = Stress::new("union_find/rollback");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(30_000 + case as u64);
        let n = rng.range(1, 15) as usize;
        let mut uf = RollbackUnionFind::new(n);
        let mut edges: Vec<(usize, usize)> = vec![]; // successful unions, in order
        let mut snaps: Vec<usize> = vec![uf.snapshot()];
        st.case();
        let label = |edges: &[(usize, usize)]| -> Vec<usize> {
            let mut lab: Vec<usize> = (0..n).collect();
            for &(x, y) in edges {
                let (a, b) = (lab[x], lab[y]);
                lab.iter_mut().for_each(|l| if *l == b { *l = a });
            }
            lab
        };
        for step in 0..60 {
            let ctx = format!("case {case} step {step} n={n}");
            match rng.below(5) {
                0..=2 => {
                    let (x, y) = (rng.below(n), rng.below(n));
                    let lab = label(&edges);
                    let merged = uf.union(x, y);
                    assert_eq!(merged, lab[x] != lab[y], "union({x},{y}) {ctx}");
                    if merged {
                        edges.push((x, y));
                    }
                }
                3 => snaps.push(uf.snapshot()),
                _ => {
                    let s = snaps[rng.below(snaps.len())];
                    if s <= edges.len() {
                        uf.rollback(s);
                        edges.truncate(s);
                        snaps.retain(|&t| t <= s);
                    }
                }
            }
            let lab = label(&edges);
            assert_eq!(uf.components(), (0..n).filter(|&i| lab[i] == i).count(), "components {ctx}");
            for x in 0..n {
                assert_eq!(uf.size(x), lab.iter().filter(|&&l| l == lab[x]).count(), "size({x}) {ctx}");
                for y in 0..n {
                    assert_eq!(uf.same(x, y), lab[x] == lab[y], "same({x},{y}) {ctx}");
                }
            }
            st.ops(n * n);
        }
    }
    st.done();
}

#[test]
fn rollback_to_empty_restores_singletons() {
    let mut uf = RollbackUnionFind::new(1000);
    for i in 0..999 {
        assert!(uf.union(i, i + 1));
    }
    assert_eq!((uf.components(), uf.size(0)), (1, 1000));
    uf.rollback(0);
    assert_eq!((uf.components(), uf.size(0), uf.same(0, 1)), (1000, 1, false));
    assert!(std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| uf.undo())).is_err());
}

// ---- range parallel -------------------------------------------------------------------------
const P: u64 = 998_244_353;

#[test]
fn range_parallel_matches_per_edge_unions() {
    let mut st = Stress::new("union_find/range_parallel");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(50_000 + case as u64);
        let n = rng.range(1, 40) as usize;
        let x: Vec<u64> = (0..n).map(|_| rng.below(P as usize) as u64).collect();
        let mut uf = RangeParallelUnionFind::new(n);
        let mut lab: Vec<usize> = (0..n).collect(); // naive component labels
        // Aggregates maintained through the merge callback: sum of x per root and X = sum of pairs.
        let mut sum: Vec<u64> = x.clone();
        let mut pairs = 0u64;
        let mut total_merges = 0usize;
        st.case();
        for step in 0..30 {
            let k = rng.below(n + 1);
            let (a, b) = (rng.below(n - k + 1), rng.below(n - k + 1));
            uf.unite_range(a, b, k, |absorbed, root| {
                pairs = (pairs + sum[absorbed] * sum[root]) % P;
                sum[root] = (sum[root] + sum[absorbed]) % P;
                total_merges += 1;
            });
            for i in 0..k {
                let (p, q) = (lab[a + i], lab[b + i]);
                lab.iter_mut().for_each(|l| if *l == q { *l = p });
            }
            let mut want = 0u64;
            for i in 0..n {
                for j in i + 1..n {
                    if lab[i] == lab[j] {
                        want = (want + x[i] * x[j]) % P;
                    }
                }
            }
            let ctx = format!("case {case} step {step} n={n} a={a} b={b} k={k}");
            assert_eq!(pairs, want, "pair sum X {ctx}");
            let comps = (0..n).filter(|&i| lab[i] == i).count();
            assert_eq!(total_merges, n - comps, "one callback per vertex-level merge {ctx}");
            for i in 0..n {
                assert_eq!(uf.size(i), lab.iter().filter(|&&l| l == lab[i]).count(), "size({i}) {ctx}");
            }
            for i in 0..n {
                for j in 0..n {
                    assert_eq!(uf.same(i, j), lab[i] == lab[j], "same({i},{j}) {ctx}");
                }
            }
            st.ops(n * n);
        }
    }
    st.done();
}

#[test]
fn range_parallel_large_matches_plain_dsu() {
    let mut st = Stress::new("union_find/range_parallel_large");
    let n = 100_000;
    let mut rng = Rng::new(77);
    let mut uf = RangeParallelUnionFind::new(n);
    let mut plain = union_find::plain(n);
    st.case();
    for _ in 0..2_000 * scale() {
        let k = rng.below(1_500);
        let (a, b) = (rng.below(n - k + 1), rng.below(n - k + 1));
        uf.unite_range(a, b, k, |_, _| {});
        for i in 0..k {
            plain.merge(a + i, b + i);
        }
        st.ops(k);
    }
    for _ in 0..50_000 {
        let (i, j) = (rng.below(n), rng.below(n));
        assert_eq!(uf.same(i, j), plain.are_connected(i, j));
        assert_eq!(uf.size(i), plain.get_size(i));
    }
    st.done();
}
