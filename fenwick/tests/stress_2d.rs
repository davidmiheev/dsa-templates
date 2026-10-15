//! Stress test: `Fenwick2D` against a plain list of weighted points.

use fenwick::Fenwick2D;
use std::panic::{catch_unwind, AssertUnwindSafe};

#[path = "../../stress_common.rs"]
mod common;
use common::{scale, Rng, Stress};

#[test]
fn matches_naive_rectangle_sums() {
    let mut st = Stress::new("fenwick/2d-naive");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64);
        let c = rng.range(1, 12);
        let cells: Vec<(i64, i64)> = (0..rng.range(0, 30)).map(|_| (rng.range(-c, c), rng.range(-c, c))).collect();
        let mut f = Fenwick2D::<i64, i64, i64>::new(&cells);
        let mut added: Vec<(i64, i64, i64)> = Vec::new();
        st.case();
        for step in 0..80 {
            let ctx = format!("case {case} step {step}");
            if !cells.is_empty() && rng.chance(1, 2) {
                let (x, y) = cells[rng.below(cells.len())];
                let w = rng.range(-50, 50);
                f.add(x, y, w);
                added.push((x, y, w));
            } else {
                let (xl, xh, yl, yh) = (rng.range(-c - 2, c + 2), rng.range(-c - 2, c + 2), rng.range(-c - 2, c + 2), rng.range(-c - 2, c + 2));
                let want: i64 = added.iter().filter(|p| xl <= p.0 && p.0 < xh && yl <= p.1 && p.1 < yh).map(|p| p.2).sum();
                assert_eq!(f.rect(xl, xh, yl, yh), want, "rect({xl}, {xh}, {yl}, {yh}) {ctx}");
                let want: i64 = added.iter().filter(|p| p.0 < xh && p.1 < yh).map(|p| p.2).sum();
                assert_eq!(f.prefix(xh, yh), want, "prefix({xh}, {yh}) {ctx}");
            }
            st.ops(1);
        }
    }
    st.done();
}

#[test]
fn large_scale() {
    let mut st = Stress::new("fenwick/2d-large");
    let n = 100_000;
    let mut rng = Rng::new(11);
    let cells: Vec<(i64, i64)> = (0..n).map(|_| (rng.range(0, 1_000_000_000), rng.range(0, 1_000_000_000))).collect();
    let mut f = Fenwick2D::<i64, i64, i64>::new(&cells);
    let mut sum_all = 0i64;
    st.case();
    for _ in 0..n * scale() {
        let (x, y) = cells[rng.below(n)];
        let w = rng.range(0, 1000);
        f.add(x, y, w);
        sum_all += w;
        st.ops(1);
    }
    assert_eq!(f.rect(0, 1_000_000_001, 0, 1_000_000_001), sum_all);
    st.done();
}

#[test]
fn panics_on_undeclared_cell() {
    let mut f = Fenwick2D::<u8, u8, i32>::new(&[(1, 1), (2, 2)]);
    assert!(catch_unwind(AssertUnwindSafe(|| f.add(1, 2, 5))).is_err(), "declared x, undeclared y");
    assert!(catch_unwind(AssertUnwindSafe(|| f.add(9, 1, 5))).is_err(), "undeclared x");
    let e = Fenwick2D::<u8, u8, i32>::new(&[]);
    assert_eq!((e.prefix(5, 5), e.rect(0, 9, 0, 9)), (0, 0));
}

#[test]
fn dominance_sums_match_naive() {
    use fenwick::{dominance_sums, Event};
    let mut st = Stress::new("fenwick/dominance");
    for case in 0..500 * scale() {
        let mut rng = Rng::new(case as u64 + 1000);
        let c = rng.range(1, 10);
        let events: Vec<Event<i64, i64, i64>> = (0..rng.range(0, 60))
            .map(|_| {
                let (x, y) = (rng.range(-c, c), rng.range(-c, c));
                if rng.chance(1, 2) { Event::Add { x, y, w: rng.range(-9, 9) } } else { Event::Query { x, y } }
            })
            .collect();
        let mut want = Vec::new();
        for (t, e) in events.iter().enumerate() {
            if let Event::Query { x, y } = *e {
                let s: i64 = events[..t]
                    .iter()
                    .filter_map(|p| match *p {
                        Event::Add { x: px, y: py, w } if px < x && py < y => Some(w),
                        _ => None,
                    })
                    .sum();
                want.push(s);
            }
        }
        st.case();
        assert_eq!(dominance_sums(&events), want, "case {case}");
        st.ops(events.len());
    }
    st.done();
}

#[test]
fn dominance_sums_large() {
    use fenwick::{dominance_sums, Event};
    let mut st = Stress::new("fenwick/dominance-large");
    let mut rng = Rng::new(5);
    let n = 400_000 * scale();
    let events: Vec<Event<i64, i64, i64>> = (0..n)
        .map(|_| {
            let (x, y) = (rng.range(0, 1_000_000_000), rng.range(0, 1_000_000_000));
            if rng.chance(1, 2) { Event::Add { x, y, w: 1 } } else { Event::Query { x, y } }
        })
        .collect();
    st.case();
    let ans = dominance_sums(&events);
    // Spot-check a few answers against a direct scan.
    let queries: Vec<usize> = (0..n).filter(|&t| matches!(events[t], Event::Query { .. })).collect();
    for k in (0..queries.len()).step_by(queries.len() / 20 + 1) {
        let t = queries[k];
        let Event::Query { x, y } = events[t] else { unreachable!() };
        let want = events[..t].iter().filter(|e| matches!(**e, Event::Add { x: px, y: py, .. } if px < x && py < y)).count() as i64;
        assert_eq!(ans[k], want, "query #{k}");
    }
    st.ops(n);
    st.done();
}
