//! `c[k] = min_{i + j = k} a[i] + b[j]` when the structure of `a` (and `b`) allows more than the
//! `O(nm)` double loop.
//!
//! * [`min_plus_convex_convex`]: both convex — merge the two walks (Minkowski sum), `O(n + m)`.
//! * [`min_plus_convex_arbitrary`]: `a` convex, `b` arbitrary. The matrix `M[k][j] = a[k-j] + b[j]`
//!   (infinite outside the band `0 <= k-j < n`) is Monge, so the column of the row minimum is
//!   monotone in `k` and divide-and-conquer finds all minima in `O((n + m) log(n + m))`.
//! * [`min_plus_concave_arbitrary`]: `a` concave. Inside a triangle or rectangle that is entirely
//!   valid the matrix is *inverse* Monge, but the band boundary breaks it, so `b` is cut into blocks
//!   of length `n` and every block is handled by recursive halving of its triangles: a fully valid
//!   rectangle (monotone minima) plus two smaller triangles. `O((n + m) log n)`.
//!
//! Weights are any totally ordered additive type (`i64`, ...).
//!
//! # Example
//!
//! ```
//! use min_plus::{min_plus_concave_arbitrary, min_plus_convex_arbitrary, min_plus_convex_convex};
//!
//! assert_eq!(min_plus_convex_convex(&[0, 1, 3], &[0, 2, 5]), vec![0, 1, 3, 5, 8]);
//! assert_eq!(min_plus_convex_arbitrary(&[0, 1, 4], &[5, 0, 7]), vec![5, 0, 1, 4, 11]);
//! assert_eq!(min_plus_concave_arbitrary(&[0, 3, 4], &[5, 0, 7]), vec![5, 0, 3, 4, 11]);
//! ```

use std::ops::Add;

/// Min-plus convolution of two convex sequences (`a[i+1] - a[i]` non-decreasing, likewise `b`).
/// Both must be non-empty.
pub fn min_plus_convex_convex<W: Copy + Ord + Add<Output = W>>(a: &[W], b: &[W]) -> Vec<W> {
    assert!(!a.is_empty() && !b.is_empty(), "sequences must be non-empty");
    let (n, m) = (a.len(), b.len());
    let mut c = Vec::with_capacity(n + m - 1);
    let (mut i, mut j) = (0, 0);
    c.push(a[0] + b[0]);
    while i + j < n + m - 2 {
        if i + 1 < n && (j + 1 >= m || a[i + 1] + b[j] <= a[i] + b[j + 1]) {
            i += 1;
        } else {
            j += 1;
        }
        c.push(a[i] + b[j]);
    }
    c
}

/// Min-plus convolution of a convex `a` with an arbitrary `b` (both non-empty).
pub fn min_plus_convex_arbitrary<W: Copy + Ord + Add<Output = W>>(a: &[W], b: &[W]) -> Vec<W> {
    assert!(!a.is_empty() && !b.is_empty(), "sequences must be non-empty");
    let (n, m) = (a.len(), b.len());
    let cell = |k: usize, j: usize| -> Option<W> {
        if j <= k && k - j < n {
            Some(a[k - j] + b[j])
        } else {
            None
        }
    };
    let mut c: Vec<Option<W>> = vec![None; n + m - 1];
    // Leftmost row minima, columns non-decreasing in the row (explicit stack: no recursion).
    let mut stack = vec![(0usize, n + m - 1, 0usize, m - 1)];
    while let Some((lo, hi, clo, chi)) = stack.pop() {
        if lo >= hi {
            continue;
        }
        let mid = (lo + hi) / 2;
        let (mut best, mut arg) = (cell(mid, clo), clo);
        for j in clo + 1..=chi {
            let v = cell(mid, j);
            let better = match (v, best) {
                (Some(x), Some(y)) => x < y,
                (Some(_), None) => true,
                _ => false,
            };
            if better {
                best = v;
                arg = j;
            }
        }
        c[mid] = best;
        stack.push((lo, mid, clo, arg));
        stack.push((mid + 1, hi, arg, chi));
    }
    c.into_iter().map(|x| x.expect("every row has a valid cell")).collect()
}

/// Row minima of the fully valid rectangle `rows r0..r1`, `cols c0..c1` whose columns of the
/// minimum are non-increasing in the row (inverse Monge): `upd(row, value)` receives each minimum.
fn rect_min_nonincreasing<W: Copy + Ord>(
    r0: usize,
    r1: usize,
    c0: usize,
    c1: usize,
    f: &dyn Fn(usize, usize) -> W,
    upd: &mut dyn FnMut(usize, W),
) {
    let mut stack = vec![(r0, r1, c0, c1 - 1)];
    while let Some((lo, hi, clo, chi)) = stack.pop() {
        if lo >= hi {
            continue;
        }
        let mid = (lo + hi) / 2;
        let (mut best, mut arg) = (f(mid, clo), clo);
        for j in clo + 1..=chi {
            let v = f(mid, j);
            if v <= best {
                best = v;
                arg = j; // rightmost minimum
            }
        }
        upd(mid, best);
        stack.push((lo, mid, arg, chi)); // smaller rows: columns to the right
        stack.push((mid + 1, hi, clo, arg));
    }
}

/// Lower triangle `lo <= j <= k < hi` of an inverse-Monge matrix, by halving.
fn lower_triangle<W: Copy + Ord>(lo: usize, hi: usize, f: &dyn Fn(usize, usize) -> W, upd: &mut dyn FnMut(usize, W)) {
    if hi - lo <= 16 {
        for k in lo..hi {
            for j in lo..=k {
                upd(k, f(k, j));
            }
        }
        return;
    }
    let mid = (lo + hi) / 2;
    rect_min_nonincreasing(mid, hi, lo, mid, f, upd);
    lower_triangle(lo, mid, f, upd);
    lower_triangle(mid, hi, f, upd);
}

/// Min-plus convolution of a concave `a` with an arbitrary `b` (both non-empty).
pub fn min_plus_concave_arbitrary<W: Copy + Ord + Add<Output = W>>(a: &[W], b: &[W]) -> Vec<W> {
    assert!(!a.is_empty() && !b.is_empty(), "sequences must be non-empty");
    let (n, m) = (a.len(), b.len());
    let c = std::cell::RefCell::new(vec![None::<W>; n + m - 1]);
    let relax = |k: usize, v: W| {
        let mut c = c.borrow_mut();
        if c[k].map_or(true, |x| v < x) {
            c[k] = Some(v);
        }
    };
    for j0 in (0..m).step_by(n) {
        let j1 = (j0 + n).min(m);
        let t = j1 - j0;
        // rows j0..j1: only columns j <= k are valid, and k - j < t <= n
        lower_triangle(j0, j1, &|k, j| a[k - j] + b[j], &mut |k, v| relax(k, v));
        if n >= 2 {
            // rows k = j1 + n - 2 - kt (kt in 0..n-1) with column j = j1 - 1 - jt: valid iff jt <= kt,
            // value a[n - 1 - (kt - jt)] + b[j1 - 1 - jt]
            let f = |kt: usize, jt: usize| a[n - 1 - (kt - jt)] + b[j1 - 1 - jt];
            let mut upd = |kt: usize, v: W| relax(j1 + n - 2 - kt, v);
            lower_triangle(0, t, &f, &mut upd);
            if n - 1 > t {
                rect_min_nonincreasing(t, n - 1, 0, t, &f, &mut upd);
            }
        }
    }
    c.into_inner().into_iter().map(|x| x.expect("every output is reachable")).collect()
}
