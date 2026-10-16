//! Classical polynomial transforms that are one or two convolutions: Taylor shift, shifting the
//! sampling points of a polynomial, evaluation and interpolation on a geometric sequence
//! (chirp-z), prefix sums of a polynomial (Bernoulli numbers) and products of many polynomials.

use crate::fps::{self, mul};
use modint::StaticModInt;
use std::collections::BinaryHeap;
use std::cmp::Reverse;

type Mint<const M: u32> = StaticModInt<M>;

/// `(k!, 1/k!)` for `k` in `0..=n` (requires `n < M`).
pub fn factorials<const M: u32>(n: usize) -> (Vec<Mint<M>>, Vec<Mint<M>>) {
    assert!((n as u64) < M as u64, "factorials beyond the modulus are not invertible");
    let mut fact = vec![Mint::<M>::one(); n + 1];
    for i in 1..=n {
        fact[i] = fact[i - 1] * Mint::from(i);
    }
    let mut inv = vec![Mint::<M>::one(); n + 1];
    inv[n] = fact[n].inv();
    for i in (1..=n).rev() {
        inv[i - 1] = inv[i] * Mint::from(i);
    }
    (fact, inv)
}

/// The coefficients of `f(x + c)`.
pub fn taylor_shift<const M: u32>(f: &[Mint<M>], c: Mint<M>) -> Vec<Mint<M>> {
    let n = f.len();
    if n == 0 {
        return Vec::new();
    }
    let (fact, inv_fact) = factorials::<M>(n);
    // g_k k! = sum_j (f_j j!) (c^(j-k) / (j-k)!)
    let a: Vec<Mint<M>> = (0..n).rev().map(|j| f[j] * fact[j]).collect();
    let mut b = Vec::with_capacity(n);
    let mut pw = Mint::<M>::one();
    for &fi in inv_fact.iter().take(n) {
        b.push(pw * fi);
        pw *= c;
    }
    let conv = mul(&a, &b);
    (0..n).map(|k| conv[n - 1 - k] * inv_fact[k]).collect()
}

/// Given `f(0), f(1), ..., f(n-1)` of a polynomial of degree `< n`, the values `f(c), f(c+1), ...,
/// f(c + m - 1)` (arguments taken modulo `M`).
pub fn shift_sampling_points<const M: u32>(values: &[Mint<M>], c: Mint<M>, m: usize) -> Vec<Mint<M>> {
    let n = values.len();
    assert!(n >= 1, "need at least one sample");
    let p = M as u64;
    let (fact, inv_fact) = factorials::<M>(n);
    // Lagrange weights w_j = y_j (-1)^(n-1-j) / (j! (n-1-j)!)
    let w: Vec<Mint<M>> = (0..n)
        .map(|j| {
            let t = values[j] * inv_fact[j] * inv_fact[n - 1 - j];
            if (n - 1 - j) % 2 == 1 { -t } else { t }
        })
        .collect();
    let _ = &fact;
    let mut result = vec![Mint::<M>::zero(); m];
    let x_of = |i: usize| -> u64 { (c.val() as u64 + i as u64) % p };
    let mut i = 0;
    while i < m {
        let x = x_of(i);
        if x < n as u64 {
            result[i] = values[x as usize];
            i += 1;
            continue;
        }
        // maximal run of consecutive i with x_of(i) in [n, p-1] and no wrap-around
        let mut j = i;
        while j < m && x_of(j) >= n as u64 && (j == i || x_of(j) == x_of(j - 1) + 1) {
            j += 1;
        }
        let len = j - i;
        let s = x; // first x of the run
        // 1 / (s - n + 1 + t) for t in 0..len + n - 1, by batch inversion
        let count = len + n - 1;
        let base = s + 1 - n as u64; // >= 1
        let mut prefix = vec![Mint::<M>::one(); count + 1];
        for t in 0..count {
            prefix[t + 1] = prefix[t] * Mint::new(base + t as u64);
        }
        let mut inv_all = prefix[count].inv();
        let mut iv = vec![Mint::<M>::zero(); count];
        for t in (0..count).rev() {
            iv[t] = inv_all * prefix[t];
            inv_all *= Mint::new(base + t as u64);
        }
        let conv = mul(&w, &iv);
        // P(x) = prod_{t<n} (x - t) at the first point, then P(x+1) = P(x) (x+1) / (x+1-n)
        let mut pval = Mint::<M>::one();
        for t in 0..n as u64 {
            pval *= Mint::new(s - t);
        }
        for k in 0..len {
            result[i + k] = pval * conv[n - 1 + k];
            pval = pval * Mint::new(s + k as u64 + 1) * iv[k];
        }
        i = j;
    }
    result
}

/// `f(a), f(a r), f(a r^2), ..., f(a r^(m-1))` by Bluestein's trick: `r^(ij)` is
/// `r^C(i+j,2) / (r^C(i,2) r^C(j,2))`, so the evaluation is one correlation. `O((n + m) log)`.
pub fn eval_geometric<const M: u32>(f: &[Mint<M>], a: Mint<M>, r: Mint<M>, m: usize) -> Vec<Mint<M>> {
    let n = f.len();
    let mut result = vec![Mint::<M>::zero(); m];
    if n == 0 || m == 0 {
        return result;
    }
    if r.val() == 0 {
        // points a, 0, 0, ...
        result[0] = f.iter().rev().fold(Mint::zero(), |acc, &c| acc * a + c);
        for x in result.iter_mut().skip(1) {
            *x = f[0];
        }
        return result;
    }
    let rinv = r.inv();
    let tri_pow = |base: Mint<M>, len: usize| -> Vec<Mint<M>> {
        // base^C(t, 2) for t in 0..len
        let mut v = Vec::with_capacity(len);
        let (mut cur, mut step) = (Mint::<M>::one(), Mint::<M>::one());
        for _ in 0..len {
            v.push(cur);
            cur *= step; // C(t+1, 2) = C(t, 2) + t
            step *= base;
        }
        v
    };
    let u = tri_pow(r, n + m - 1);
    let v = tri_pow(rinv, n.max(m));
    let mut apow = Mint::<M>::one();
    let mut left = Vec::with_capacity(n);
    for j in 0..n {
        left.push(f[j] * apow * v[j]);
        apow *= a;
    }
    left.reverse();
    let conv = mul(&left, &u);
    for i in 0..m {
        result[i] = v[i] * conv[n - 1 + i];
    }
    result
}

/// The coefficients of the polynomial of degree `< n` taking the value `y[i]` at `a q^i` for
/// `i < n` (the points must be distinct), via the Newton form on a geometric progression and the
/// q-binomial theorem: two convolutions.
pub fn interpolate_geometric<const M: u32>(y: &[Mint<M>], a: Mint<M>, q: Mint<M>) -> Vec<Mint<M>> {
    let n = y.len();
    if n <= 1 {
        return y.to_vec();
    }
    assert!(a.val() != 0, "n >= 2 distinct points need a non-zero start");
    // F_s = prod_{u=1}^{s} (q^u - 1), q^{C(s,2)}
    let mut f_prod = vec![Mint::<M>::one(); n];
    let mut qpow = Mint::<M>::one();
    for s in 1..n {
        qpow *= q;
        f_prod[s] = f_prod[s - 1] * (qpow - Mint::one());
        assert!(f_prod[s].val() != 0, "the points are not distinct");
    }
    let mut inv_f = vec![Mint::<M>::zero(); n];
    inv_f[n - 1] = f_prod[n - 1].inv();
    for s in (1..n).rev() {
        inv_f[s - 1] = inv_f[s] * (q.pow(s as u64) - Mint::one());
    }
    let tri = |base: Mint<M>| -> Vec<Mint<M>> {
        let mut v = Vec::with_capacity(n);
        let (mut cur, mut step) = (Mint::<M>::one(), Mint::<M>::one());
        for _ in 0..n {
            v.push(cur);
            cur *= step;
            step *= base;
        }
        v
    };
    let qp = tri(q);
    let qp_inv = if q.val() == 0 { vec![Mint::<M>::one(); n] } else { tri(q.inv()) };
    // d_k = a^{-k} q^{-C(k,2)} sum_i (y_i / F_i) (-1)^{k-i} q^{C(k-i,2)} / F_{k-i}
    let e: Vec<Mint<M>> = (0..n).map(|i| y[i] * inv_f[i]).collect();
    let g: Vec<Mint<M>> = (0..n)
        .map(|s| {
            let t = qp[s] * inv_f[s];
            if s % 2 == 1 { -t } else { t }
        })
        .collect();
    let conv = mul(&e, &g);
    let a_inv = a.inv();
    let mut apow = Mint::<M>::one();
    let mut h = Vec::with_capacity(n); // H_k = d_k F_k
    for k in 0..n {
        let d = apow * qp_inv[k] * conv[k];
        h.push(d * f_prod[k]);
        apow *= a_inv;
    }
    // c_j = (1 / F_j) sum_{k >= j} H_k L_{k-j},  L_s = (-a)^s q^{C(s,2)} / F_s
    let mut l = Vec::with_capacity(n);
    let mut pw = Mint::<M>::one();
    for s in 0..n {
        l.push(pw * qp[s] * inv_f[s]);
        pw *= -a;
    }
    h.reverse();
    let corr = mul(&h, &l);
    (0..n).map(|j| inv_f[j] * corr[n - 1 - j]).collect()
}

/// The Bernoulli numbers `B_0, ..., B_n` (convention `B_1 = -1/2`): `B_j / j!` are the
/// coefficients of `x / (e^x - 1)`, one power series inverse. `O(n log n)`.
pub fn bernoulli_numbers<const M: u32>(n: usize) -> Vec<Mint<M>> {
    let (fact, inv_fact) = factorials::<M>(n + 2);
    let denom: Vec<Mint<M>> = (0..n + 2).map(|t| inv_fact[t + 1]).collect();
    let b = fps::inv(&denom, n + 1);
    (0..=n).map(|j| b[j] * fact[j]).collect()
}

/// `g` of degree `<= n` with `g(m) = sum_{i < m} f(i)` for every integer `m >= 0`, where
/// `f = a_0 + a_1 x + ...` has `n` coefficients. From Faulhaber's formula with Bernoulli numbers:
/// `g_m = (1/m!) sum_j (B_j / j!) a_{m+j-1} (m+j-1)!`.
pub fn prefix_sum_polynomial<const M: u32>(f: &[Mint<M>]) -> Vec<Mint<M>> {
    let n = f.len();
    let (fact, inv_fact) = factorials::<M>(n + 2);
    // B_j / j! = [x^j] x / (e^x - 1) = 1 / sum_t x^t / (t+1)!
    let denom: Vec<Mint<M>> = (0..n + 2).map(|t| inv_fact[t + 1]).collect();
    let bern = fps::inv(&denom, n + 1); // B_j / j! for j in 0..=n
    // g_m = (1/m!) sum_{j=0}^{n-m} bern_j * (a_{m+j-1} (m+j-1)!)
    let t: Vec<Mint<M>> = (0..n).map(|i| f[i] * fact[i]).collect(); // t_i = a_i i!
    let mut tr = t.clone();
    tr.reverse();
    let conv = mul(&tr, &bern);
    let mut g = vec![Mint::<M>::zero(); n + 1];
    for m in 1..=n {
        // sum_j bern_j t_{m+j-1}: reversed index of t_{m-1+j} is n-1-(m-1+j) = n-m-j
        // conv[(n-m)] = sum_j tr[n-m-j] bern[j] = sum_j t[m-1+j] bern[j]
        g[m] = conv[n - m] * inv_fact[m];
    }
    g
}

/// The product of all the polynomials, merging the two shortest first (`O(D log D log k)` for total
/// degree `D`). The empty product is `[1]`.
pub fn product_all<const M: u32>(polys: Vec<Vec<Mint<M>>>) -> Vec<Mint<M>> {
    if polys.is_empty() {
        return vec![Mint::one()];
    }
    let mut store: Vec<Vec<Mint<M>>> = polys;
    let mut heap: BinaryHeap<Reverse<(usize, usize)>> = store.iter().enumerate().map(|(i, p)| Reverse((p.len(), i))).collect();
    while heap.len() > 1 {
        let Reverse((_, i)) = heap.pop().unwrap();
        let Reverse((_, j)) = heap.pop().unwrap();
        let (a, b) = (std::mem::take(&mut store[i]), std::mem::take(&mut store[j]));
        let c = mul(&a, &b);
        let len = c.len();
        store[i] = c;
        heap.push(Reverse((len, i)));
    }
    let Reverse((_, i)) = heap.pop().unwrap();
    std::mem::take(&mut store[i])
}
