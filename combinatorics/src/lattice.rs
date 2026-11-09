//! Nondecreasing sequences between two bounding sequences, by divide and conquer over the grid.

use convolution::convolution_mod;

const BASE_AREA: usize = 1024;

struct Ctx<const M: u32> {
    /// lower bounds of the elements (nondecreasing), length n
    l: Vec<usize>,
    /// upper bounds of the elements (nondecreasing), length n + 1 (`u[n] = u[n-1]`)
    u: Vec<usize>,
    fact: Vec<u32>,
    inv_fact: Vec<u32>,
    base_area: usize,
}

impl<const M: u32> Ctx<M> {
    fn mulm(a: u32, b: u32) -> u32 {
        (a as u64 * b as u64 % M as u64) as u32
    }

    fn binom(&self, n: usize, k: usize) -> u32 {
        Self::mulm(Self::mulm(self.fact[n], self.inv_fact[k]), self.inv_fact[n - k])
    }

    fn add(a: u32, b: u32) -> u32 {
        let s = a + b;
        if s >= M {
            s - M
        } else {
            s
        }
    }

    /// Path counts through the block of grid points `[c0, c1] x [h0, h1]`.
    ///
    /// `left[j]` / `bottom[a]` are the numbers of paths entering the point `(c0, h0 + j)` from
    /// the left, resp. `(c0 + a, h0)` from below. A point `(c, h)` is usable iff `h <= u[c]`; a
    /// step `(c - 1, h) -> (c, h)` is allowed iff `h >= l[c - 1]`. Returns the path counts at the
    /// right column `(c1, h)` and at the top row `(c, h1)`.
    fn solve(&self, c0: usize, c1: usize, h0: usize, h1: usize, left: &[u32], bottom: &[u32]) -> (Vec<u32>, Vec<u32>) {
        let (w, hh) = (c1 - c0 + 1, h1 - h0 + 1);
        let zeros = || (vec![0u32; hh], vec![0u32; w]);
        if self.u[c1] < h0 || (left.iter().all(|&x| x == 0) && bottom.iter().all(|&x| x == 0)) {
            return zeros();
        }
        if w == 1 || hh == 1 || w * hh <= self.base_area {
            return self.direct(c0, c1, h0, h1, left, bottom);
        }
        let valid_all = self.u[c0] >= h1;
        if valid_all && h0 >= self.l[c1 - 1] {
            return self.transfer(w - 1, hh - 1, left, bottom);
        }
        if h1 < self.l[c0] {
            return self.columns(c0, c1, h0, h1, left, bottom);
        }
        if w >= hh {
            let cm = (c0 + c1) / 2;
            let (right_l, top_l) = self.solve(c0, cm, h0, h1, left, &bottom[..cm - c0 + 1]);
            let left_r: Vec<u32> = (0..hh).map(|j| if h0 + j >= self.l[cm] { right_l[j] } else { 0 }).collect();
            let (right, top_r) = self.solve(cm + 1, c1, h0, h1, &left_r, &bottom[cm - c0 + 1..]);
            let mut top = top_l;
            top.extend(top_r);
            (right, top)
        } else {
            let hm = (h0 + h1) / 2;
            let (right_lo, top_lo) = self.solve(c0, c1, h0, hm, &left[..hm - h0 + 1], bottom);
            let bottom_up: Vec<u32> = (0..w).map(|a| if hm < self.u[c0 + a] { top_lo[a] } else { 0 }).collect();
            let (right_up, top) = self.solve(c0, c1, hm + 1, h1, &left[hm - h0 + 1..], &bottom_up);
            let mut right = right_lo;
            right.extend(right_up);
            (right, top)
        }
    }

    /// Cell by cell dynamic programming, `O(w h)`.
    fn direct(&self, c0: usize, c1: usize, h0: usize, h1: usize, left: &[u32], bottom: &[u32]) -> (Vec<u32>, Vec<u32>) {
        let (w, hh) = (c1 - c0 + 1, h1 - h0 + 1);
        let mut top = vec![0u32; w];
        let mut prev = vec![0u32; hh];
        let mut cur = vec![0u32; hh];
        for c in c0..=c1 {
            for j in 0..hh {
                let h = h0 + j;
                let mut v = if c == c0 {
                    left[j]
                } else if h >= self.l[c - 1] {
                    prev[j]
                } else {
                    0
                };
                if j == 0 {
                    v = Self::add(v, bottom[c - c0]);
                } else {
                    v = Self::add(v, cur[j - 1]);
                }
                cur[j] = if h > self.u[c] { 0 } else { v };
            }
            top[c - c0] = cur[hh - 1];
            std::mem::swap(&mut prev, &mut cur);
        }
        (prev, top)
    }

    /// A block where no step to the right is allowed: the columns are independent.
    fn columns(&self, c0: usize, c1: usize, h0: usize, h1: usize, left: &[u32], bottom: &[u32]) -> (Vec<u32>, Vec<u32>) {
        let (w, hh) = (c1 - c0 + 1, h1 - h0 + 1);
        let mut top = vec![0u32; w];
        let mut right = vec![0u32; hh];
        // column c0 collects the entries from the left and from below
        let mut acc = 0u32;
        let mut first = vec![0u32; hh];
        for j in 0..hh {
            acc = Self::add(acc, left[j]);
            if j == 0 {
                acc = Self::add(acc, bottom[0]);
            }
            first[j] = if h0 + j > self.u[c0] { 0 } else { acc };
        }
        top[0] = first[hh - 1];
        for c in c0 + 1..=c1 {
            if self.u[c] >= h0 {
                top[c - c0] = if self.u[c] >= h1 { bottom[c - c0] } else { 0 };
            }
        }
        for (j, r) in right.iter_mut().enumerate() {
            *r = if h0 + j <= self.u[c1] { bottom[c1 - c0] } else { 0 };
        }
        (right, top)
    }

    /// A block of usable points without blocked steps: paths counted by binomials, as
    /// convolutions. `dc = c1 - c0 >= 1`, `dh = h1 - h0 >= 1`.
    fn transfer(&self, dc: usize, dh: usize, left: &[u32], bottom: &[u32]) -> (Vec<u32>, Vec<u32>) {
        let f = &self.fact;
        let ifa = &self.inv_fact;
        // right[j] = sum_{i <= j} left[i] C(dc + j - i, dc) + sum_a bottom[a] C(dc - a + j, j)
        let k1: Vec<u32> = (0..=dh).map(|d| self.binom(dc + d, dc)).collect();
        let r1 = convolution_mod::<M>(left, &k1);
        let u: Vec<u32> = (0..=dc).map(|t| Self::mulm(bottom[dc - t], ifa[t])).collect(); // t = dc - a
        let u_rev: Vec<u32> = u.iter().rev().copied().collect();
        let fk: Vec<u32> = f[..=dc + dh].to_vec();
        let s2 = convolution_mod::<M>(&u_rev, &fk);
        let right: Vec<u32> = (0..=dh).map(|j| Self::add(r1[j], Self::mulm(ifa[j], s2[dc + j]))).collect();
        // top[s] = sum_i left[i] C(s + dh - i, s) + sum_{a <= s} bottom[a] C(s - a + dh, dh)
        let v: Vec<u32> = (0..=dh).map(|t| Self::mulm(left[dh - t], ifa[t])).collect(); // t = dh - i
        let v_rev: Vec<u32> = v.iter().rev().copied().collect();
        let s3 = convolution_mod::<M>(&v_rev, &fk);
        let k2: Vec<u32> = (0..=dc).map(|d| self.binom(d + dh, dh)).collect();
        let t4 = convolution_mod::<M>(bottom, &k2);
        let top: Vec<u32> = (0..=dc).map(|s| Self::add(Self::mulm(ifa[s], s3[dh + s]), t4[s])).collect();
        (right, top)
    }
}

/// The number of integer sequences `x_0 <= x_1 <= ... <= x_{n-1} < m` with `a_i <= x_i < b_i`,
/// modulo the NTT prime `M` (`0 <= a_i < b_i <= m`), in `O((n + m) log^2 (n + m))`.
///
/// The sequence is a monotone lattice path; the bounds, made monotone by prefix maxima of `a` and
/// suffix minima of `b`, are two staircases, and the region between them is cut recursively along
/// its longer side. Blocks lying fully inside the region are crossed with binomial convolutions.
pub fn count_increasing_sequences<const M: u32>(m: usize, a: &[usize], b: &[usize]) -> u32 {
    count_with_base::<M>(m, a, b, BASE_AREA)
}

#[doc(hidden)]
pub fn count_with_base<const M: u32>(m: usize, a: &[usize], b: &[usize], base_area: usize) -> u32 {
    let n = a.len();
    assert_eq!(n, b.len(), "length mismatch");
    assert!(a.iter().zip(b).all(|(&x, &y)| x < y && y <= m), "need a_i < b_i <= m");
    if n == 0 {
        return 1 % M;
    }
    let mut l = a.to_vec();
    for i in 1..n {
        l[i] = l[i].max(l[i - 1]);
    }
    let mut u: Vec<usize> = b.iter().map(|&y| y - 1).collect();
    for i in (0..n - 1).rev() {
        u[i] = u[i].min(u[i + 1]);
    }
    if (0..n).any(|i| l[i] > u[i]) {
        return 0;
    }
    u.push(u[n - 1]);
    let top = u[n];
    let size = n + top + 2;
    let mut fact = vec![1u32; size + 1];
    for i in 1..=size {
        fact[i] = (fact[i - 1] as u64 * i as u64 % M as u64) as u32;
    }
    let mut inv_fact = vec![1u32; size + 1];
    let mut e = M as u64 - 2;
    let (mut base, mut r) = (fact[size] as u64, 1u64);
    while e > 0 {
        if e & 1 == 1 {
            r = r * base % M as u64;
        }
        base = base * base % M as u64;
        e >>= 1;
    }
    inv_fact[size] = r as u32;
    for i in (1..=size).rev() {
        inv_fact[i - 1] = (inv_fact[i] as u64 * i as u64 % M as u64) as u32;
    }
    let ctx = Ctx::<M> { l, u, fact, inv_fact, base_area };
    let mut left = vec![0u32; top + 1];
    left[0] = 1;
    let bottom = vec![0u32; n + 1];
    let (right, _) = ctx.solve(0, n, 0, top, &left, &bottom);
    right[top]
}
