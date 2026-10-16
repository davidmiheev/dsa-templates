//! A table indexed by the values `floor(n / i)`.

/// The isqrt `r` of `n` (`r * r <= n < (r + 1)^2`).
pub fn isqrt(n: u64) -> u64 {
    let mut r = (n as f64).sqrt() as u64;
    while r * r > n {
        r -= 1;
    }
    while (r + 1) * (r + 1) <= n {
        r += 1;
    }
    r
}

/// Values attached to every element of `Q_n = { floor(n / i) : 1 <= i <= n }`.
///
/// With `r = isqrt(n)`, `Q_n` is `{1, ..., r}` together with `floor(n / i)` for `i = 1..=r`; the
/// two halves are stored as `small[v]` (`v` in `0..=r`) and `large[i]` (`i` in `1..=r`, the value
/// at `floor(n / i)`). When `floor(n / r) = r` the element `r` appears in both and both copies
/// are kept equal by the algorithms of this crate.
#[derive(Clone, Debug)]
pub struct Quotients<T> {
    n: u64,
    r: u64,
    /// `small[v]` is the entry of `v`, `v <= r`.
    pub small: Vec<T>,
    /// `large[i]` is the entry of `floor(n / i)`, `1 <= i <= r` (`large[0]` is unused).
    pub large: Vec<T>,
}

impl<T: Copy> Quotients<T> {
    /// A table for `n >= 1` with every entry `fill`.
    pub fn new(n: u64, fill: T) -> Self {
        assert!(n >= 1, "n must be positive");
        let r = isqrt(n);
        Quotients { n, r, small: vec![fill; r as usize + 1], large: vec![fill; r as usize + 1] }
    }

    /// The table of `f(v)` over all `v` in `Q_n`.
    pub fn from_fn(n: u64, mut f: impl FnMut(u64) -> T) -> Self {
        let mut q = Self::new(n, f(1));
        for v in 1..=q.r {
            q.small[v as usize] = f(v);
            q.large[v as usize] = f(n / v);
        }
        q
    }

    /// The table from the entries in increasing order of their values (see [`Quotients::ascending`]).
    pub fn from_ascending(n: u64, values: &[T]) -> Self {
        let mut q = Self::new(n, values[0]);
        let r = q.r as usize;
        let top = if n / q.r == q.r { r - 1 } else { r };
        assert_eq!(values.len(), r + top, "wrong number of values");
        q.small[1..=r].copy_from_slice(&values[..r]);
        for i in (1..=top).rev() {
            q.large[i] = values[r + (top - i)];
        }
        if top < r {
            q.large[r] = q.small[r];
        }
        q
    }

    pub fn n(&self) -> u64 {
        self.n
    }

    /// `isqrt(n)`.
    pub fn root(&self) -> u64 {
        self.r
    }

    /// The entry of `v`, which has to be an element of `Q_n` (or `<= isqrt(n)`).
    #[inline(always)]
    pub fn get(&self, v: u64) -> T {
        if v <= self.r {
            self.small[v as usize]
        } else {
            debug_assert!(self.n / (self.n / v) == v, "{v} is not a floor quotient of {}", self.n);
            self.large[(self.n / v) as usize]
        }
    }

    /// The entries in increasing order of their values (each value once): the `a_0, ..., a_{k-1}`
    /// of the Library Checker problems.
    pub fn ascending(&self) -> Vec<T> {
        let mut out: Vec<T> = self.small[1..].to_vec();
        let top = if self.n / self.r == self.r { self.r - 1 } else { self.r };
        for i in (1..=top).rev() {
            out.push(self.large[i as usize]);
        }
        out
    }
}

/// The elements of `Q_n` in increasing order.
pub fn quotient_values(n: u64) -> Vec<u64> {
    let r = isqrt(n);
    let mut out: Vec<u64> = (1..=r).collect();
    let top = if n / r == r { r - 1 } else { r };
    for i in (1..=top).rev() {
        out.push(n / i);
    }
    out
}
