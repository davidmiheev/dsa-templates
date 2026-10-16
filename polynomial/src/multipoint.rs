//! Subproduct trees: multipoint evaluation, interpolation and the conversion from the monomial
//! basis to the Newton basis, all `O(n log^2 n)`.

use crate::fps::{derivative, div_rem, inv, mul, mul_trunc, resized};
use modint::StaticModInt;

type Mint<const M: u32> = StaticModInt<M>;

/// The products `prod (x - p_i)` over the nodes of a balanced binary tree on the points (node `v`
/// covers a range split at its midpoint; children `2v` and `2v + 1`).
pub struct SubproductTree<const M: u32> {
    n: usize,
    prods: Vec<Vec<Mint<M>>>,
}

impl<const M: u32> SubproductTree<M> {
    /// Build the tree for at least one point.
    pub fn new(points: &[Mint<M>]) -> Self {
        assert!(!points.is_empty(), "need at least one point");
        let n = points.len();
        let mut t = Self { n, prods: vec![Vec::new(); 4 * n] };
        t.build(1, 0, n, points);
        t
    }

    fn build(&mut self, v: usize, l: usize, r: usize, points: &[Mint<M>]) {
        if r - l == 1 {
            self.prods[v] = vec![-points[l], Mint::one()];
            return;
        }
        let mid = (l + r) / 2;
        self.build(2 * v, l, mid, points);
        self.build(2 * v + 1, mid, r, points);
        self.prods[v] = mul(&self.prods[2 * v], &self.prods[2 * v + 1]);
    }

    /// `prod_i (x - p_i)`.
    pub fn product(&self) -> &[Mint<M>] {
        &self.prods[1]
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// `f(p_i)` for every point, by the transposed ("middle product") algorithm.
    pub fn evaluate(&self, f: &[Mint<M>]) -> Vec<Mint<M>> {
        let (n_coef, m) = (f.len(), self.n);
        if n_coef == 0 {
            return vec![Mint::zero(); m];
        }
        let rev = |p: &[Mint<M>]| -> Vec<Mint<M>> { p.iter().rev().copied().collect() };
        // T = f^R / prod(1 - p_i x) mod x^N, window of its last m coefficients (zero below index 0)
        let q_root = rev(&self.prods[1]);
        let fr = rev(f);
        let t = mul_trunc(&fr, &inv(&q_root, n_coef), n_coef);
        let u: Vec<Mint<M>> = (0..m).map(|k| if n_coef + k >= m { t[n_coef + k - m] } else { Mint::zero() }).collect();
        let mut result = vec![Mint::<M>::zero(); m];
        self.descend(1, 0, m, u, &mut result);
        result
    }

    fn descend(&self, v: usize, l: usize, r: usize, u: Vec<Mint<M>>, out: &mut [Mint<M>]) {
        if r - l == 1 {
            out[l] = u[0];
            return;
        }
        let mid = (l + r) / 2;
        let (ml, mr) = (mid - l, r - mid);
        let rev = |p: &[Mint<M>]| -> Vec<Mint<M>> { p.iter().rev().copied().collect() };
        let q_left = rev(&self.prods[2 * v]);
        let q_right = rev(&self.prods[2 * v + 1]);
        let cl = mul(&u, &q_right);
        let ul: Vec<Mint<M>> = (0..ml).map(|k| cl[mr + k]).collect();
        let cr = mul(&u, &q_left);
        let ur: Vec<Mint<M>> = (0..mr).map(|k| cr[ml + k]).collect();
        drop(cl);
        drop(cr);
        drop(u);
        self.descend(2 * v, l, mid, ul, out);
        self.descend(2 * v + 1, mid, r, ur, out);
    }

    /// The polynomial `sum_i w_i prod_{j != i} (x - p_j)`.
    pub fn weighted_sum(&self, w: &[Mint<M>]) -> Vec<Mint<M>> {
        assert_eq!(w.len(), self.n, "need one weight per point");
        self.combine(1, 0, self.n, w)
    }

    fn combine(&self, v: usize, l: usize, r: usize, w: &[Mint<M>]) -> Vec<Mint<M>> {
        if r - l == 1 {
            return vec![w[l]];
        }
        let mid = (l + r) / 2;
        let a = self.combine(2 * v, l, mid, w);
        let b = self.combine(2 * v + 1, mid, r, w);
        let x = mul(&a, &self.prods[2 * v + 1]);
        let y = mul(&b, &self.prods[2 * v]);
        let mut s = vec![Mint::<M>::zero(); x.len().max(y.len())];
        for (i, &c) in x.iter().enumerate() {
            s[i] += c;
        }
        for (i, &c) in y.iter().enumerate() {
            s[i] += c;
        }
        s.truncate(r - l);
        s
    }

    /// Newton-basis coefficients `b` with `f = sum_i b_i prod_{j < i} (x - p_j)`, for `f` with at
    /// most `len()` coefficients.
    pub fn to_newton_basis(&self, f: &[Mint<M>]) -> Vec<Mint<M>> {
        assert!(f.len() <= self.n, "too many coefficients for the number of points");
        let mut out = vec![Mint::<M>::zero(); self.n];
        self.newton(1, 0, self.n, resized(f, self.n), &mut out);
        out
    }

    fn newton(&self, v: usize, l: usize, r: usize, f: Vec<Mint<M>>, out: &mut [Mint<M>]) {
        if r - l == 1 {
            out[l] = f[0];
            return;
        }
        let mid = (l + r) / 2;
        let (q, mut rem) = div_rem(&f, &self.prods[2 * v]);
        rem.resize(mid - l, Mint::zero());
        let mut q = q;
        q.resize(r - mid, Mint::zero());
        self.newton(2 * v, l, mid, rem, out);
        self.newton(2 * v + 1, mid, r, q, out);
    }
}

/// `f(p_i)` for all `i` (`O((n + m) log^2)`).
pub fn multipoint_evaluate<const M: u32>(f: &[Mint<M>], points: &[Mint<M>]) -> Vec<Mint<M>> {
    if points.is_empty() {
        return Vec::new();
    }
    SubproductTree::new(points).evaluate(f)
}

/// The unique polynomial of degree `< n` through the points `(x_i, y_i)` (distinct `x_i`).
///
/// # Panics
/// Panics if two `x_i` coincide.
pub fn interpolate<const M: u32>(xs: &[Mint<M>], ys: &[Mint<M>]) -> Vec<Mint<M>> {
    assert_eq!(xs.len(), ys.len(), "need one value per point");
    if xs.is_empty() {
        return Vec::new();
    }
    let tree = SubproductTree::new(xs);
    let dp = derivative(tree.product());
    let slopes = tree.evaluate(&dp);
    let w: Vec<Mint<M>> = ys
        .iter()
        .zip(&slopes)
        .map(|(&y, &s)| {
            assert!(s.val() != 0, "the interpolation points must be distinct");
            y / s
        })
        .collect();
    let mut f = tree.weighted_sum(&w);
    f.resize(xs.len(), Mint::zero());
    f
}

/// Coefficients `b` of `f = sum_i b_i prod_{j < i} (x - p_j)` (the Newton basis on `p_0, p_1, ...`).
pub fn monomial_to_newton<const M: u32>(f: &[Mint<M>], points: &[Mint<M>]) -> Vec<Mint<M>> {
    if points.is_empty() {
        return Vec::new();
    }
    SubproductTree::new(points).to_newton_basis(f)
}

