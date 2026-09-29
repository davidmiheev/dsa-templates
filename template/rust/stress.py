#!/usr/bin/env python3
"""Stress test for the Rust competition template (template/rust/main.rs).

What is checked (all black-box, through real processes):

1. The example solution in main.rs against a brute-force shortest-path solver,
   in both FastScan modes (INTERACTIVE = false / true), with input tokens
   spread over lines with random whitespace (spaces, tabs, CRLF, form feeds, ...).
2. FastScan/Out extracted from main.rs against known token streams: arbitrary
   non-ASCII tokens, edge-case integers, trailing junk, and malformed input
   (which must make the program fail, not silently continue).
3. Interactive mode against a real interactor (binary search over a hidden
   number). Negative controls check that every mismatched mode combination
   really deadlocks, i.e. that the INTERACTIVE flag is needed.
4. Timings for n = q = 2*10^5 in both modes.

Usage: python3 template/rust/stress.py [--quick]
Needs cargo and network access for the ac-library-rs crate the template imports.
"""
import argparse
import heapq
import random
import re
import statistics
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEMPLATE = (HERE / "main.rs").read_text()
FLAG = "const INTERACTIVE: bool = false;"

m = re.search(r"// ---- fast I/O: begin[^\n]*\n(.*?)// ---- fast I/O: end ----", TEMPLATE, re.S)
assert m, "fast I/O markers not found in main.rs"
IO = m.group(1)
HDR = "use std::io::{self, Read, BufRead, Write, BufWriter, StdinLock, StdoutLock};\nuse std::fmt::Display;\n"
assert FLAG in TEMPLATE, "expected `%s` in main.rs" % FLAG

ECHO = HDR + IO + """
fn main() {
    let interactive = std::env::args().nth(1).as_deref() == Some("i");
    let mut cin = FastScan::new(interactive);
    let mut out = Out::new(interactive);
    let m: usize = cin.next();
    let v: Vec<String> = cin.read_vec(m);
    out.print_vec(&v, "|");
    let k: usize = cin.next();
    let w: Vec<i64> = cin.read_vec(k);
    out.print_vec(&w, " ");
}
"""

# Interactive toy problem: hidden secret in 1..=n. `? x` is answered by "<1 if secret <= x else 0> <query count>".
GUESS = HDR + IO + """
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (scan_i, out_i) = (a[1] == "1", a[2] == "1");
    let mut cin = FastScan::new(scan_i);
    let mut out = Out::new(out_i);
    let n: i64 = cin.next();
    let (mut lo, mut hi) = (1i64, n);
    while lo < hi {
        let mid = (lo + hi) / 2;
        out.println(format_args!("? {}", mid));
        let flag: i64 = cin.next();
        let _count: i64 = cin.next();
        if flag == 1 { hi = mid } else { lo = mid + 1 }
    }
    out.println(format_args!("! {}", lo));
}
"""

CARGO = """[package]
name = "template_stress"
version = "0.1.0"
edition = "2021"
[dependencies]
ac-library-rs = "0.1"
[profile.release]
opt-level = 3
"""

failures = 0


def report(name, cases, bad, extra=""):
    global failures
    failures += bad
    print(f"[stress] {name:<44} cases={cases:>6} failures={bad:>3} {extra}")


def run(exe, data, args=(), timeout=60):
    p = subprocess.run([str(exe), *args], input=data, capture_output=True, timeout=timeout)
    return p.returncode, p.stdout


SEPS = [" ", " ", "  ", "\t", "\n", "\r\n", " \n", "\n\n", "\x0b", "\x0c", " \t \r\n"]


def fmt(tokens, r, mode):
    """Join tokens into an input file with a chosen or random layout."""
    if mode == "canon":
        s = " ".join(tokens) + "\n"
    elif mode == "oneline":
        s = " ".join(tokens)
    elif mode == "perline":
        s = "\n".join(tokens) + "\n"
    else:
        s = "".join(t + r.choice(SEPS) for t in tokens[:-1]) + tokens[-1]
        if r.random() < 0.5:
            s = r.choice(SEPS) + s
        if r.random() < 0.5:
            s += r.choice(SEPS)
    return s.encode()


MODES = ["canon", "oneline", "perline", "mixed", "mixed"]


def graph(n, a, b):
    """Cycle 0..n-1 with edge i -- (i+1)%n of weight a[i], plus hub n joined to i with weight b[i]."""
    adj = [[] for _ in range(n + 1)]
    for i in range(n):
        j = (i + 1) % n
        adj[i].append((j, a[i]))
        adj[j].append((i, a[i]))
        adj[i].append((n, b[i]))
        adj[n].append((i, b[i]))
    return adj


def dijkstra(adj, s):
    d = [float("inf")] * len(adj)
    d[s] = 0
    pq = [(0, s)]
    while pq:
        du, u = heapq.heappop(pq)
        if du > d[u]:
            continue
        for v, w in adj[u]:
            if du + w < d[v]:
                d[v] = du + w
                heapq.heappush(pq, (du + w, v))
    return d


def brute(n, a, b, queries):
    adj, cache, out = graph(n, a, b), {}, []
    for l, r in queries:
        if l - 1 not in cache:
            cache[l - 1] = dijkstra(adj, l - 1)
        out.append(cache[l - 1][r - 1])
    return out


def instance(r, n, q, big=False):
    val = (lambda: r.randint(1, 10**9)) if big else (lambda: r.choice([0, 1, 2, 5, r.randint(0, 50), r.randint(0, 10**9)]))
    a, b = [val() for _ in range(n)], [val() for _ in range(n)]
    queries = []
    for _ in range(q):
        while True:
            l, rr = r.randint(1, n + 1), r.randint(1, n + 1)
            if not (l == rr == n + 1):  # the template solution does not define hub-to-hub queries
                break
        queries.append((l, rr))
    return a, b, queries


def tokens_of(n, a, b, queries):
    return [str(n), str(len(queries))] + list(map(str, a)) + list(map(str, b)) + [str(x) for p in queries for x in p]


def test_solution(bins, quick):
    r = random.Random(1)
    cases = bad = 0
    total = 150 if quick else 500
    for it in range(total):
        n = r.choice([1, 2, 3, 4, 5, 8, 13, 40]) if it < total - 10 else r.choice([150, 300])
        q = r.randint(0, 25) if n < 100 else 200
        a, b, queries = instance(r, n, q)
        expect = brute(n, a, b, queries)
        for exe in bins:
            for mode in MODES:
                cases += 1
                rc, out = run(exe, fmt(tokens_of(n, a, b, queries), r, mode))
                got = list(map(int, out.split())) if rc == 0 else None
                if got != expect:
                    bad += 1
                    if bad <= 5:
                        print(f"  MISMATCH {exe.name} mode={mode} n={n} q={q}\n    expected={expect[:8]} got={None if got is None else got[:8]} rc={rc}")
    report("example solution vs brute-force Dijkstra", cases, bad, "(both modes, random whitespace)")


ALPHA = list("abcXYZ019-+.,;:!?#$%&*()[]{}<>=/\\_~^|'\"") + ["é", "日本", "😀", "–", " "]


def test_scanner(echo, quick):
    r = random.Random(2)
    cases = bad = 0
    total = 300 if quick else 1000
    tok = lambda: "".join(r.choice(ALPHA) for _ in range(r.randint(1, 10)))
    num = lambda: r.choice([str(r.randint(-10**18, 10**18)), "0", "-0", "+7", "007", "9223372036854775807", "-9223372036854775808"])
    for _ in range(total):
        m, k = r.randint(0, 25), r.randint(0, 25)
        words, nums = [tok() for _ in range(m)], [num() for _ in range(k)]
        toks = [str(m)] + words + [str(k)] + nums
        if r.random() < 0.2:
            toks += [tok() for _ in range(r.randint(1, 3))]  # trailing junk must be ignored
        expect = ("|".join(words) + "\n" + " ".join(str(int(x)) for x in nums) + "\n").encode()
        for arg in ("b", "i"):
            for mode in MODES:
                cases += 1
                rc, out = run(echo, fmt(toks, r, mode), [arg])
                if rc != 0 or out != expect:
                    bad += 1
                    if bad <= 5:
                        print(f"  MISMATCH mode={arg}/{mode} rc={rc} out={out[:60]!r} expect={expect[:60]!r}")
    report("FastScan/Out token round-trip", cases, bad, "(both modes, incl. non-ASCII tokens)")

    malformed = [b"", b"   \n", b"3\na b\n", b"2\nx y\n1\n", b"1\nx\n1\n12x\n", b"1\nx\n1\n9223372036854775808\n",
                 b"1\nx\n1\n--5\n", b"1\nx\n1\n1.5\n", b"-1\n", b"abc\n", b"1\n\xff\n1\n5\n"]
    bad = 0
    for d in malformed:
        for arg in ("b", "i"):
            rc, _ = run(echo, d, [arg])
            if rc == 0:
                bad += 1
                print(f"  ACCEPTED malformed input {d!r} in mode {arg}")
    report("malformed input must fail", len(malformed) * 2, bad)


def interact(exe, args, secret, n=10**9, limit=3.0):
    p = subprocess.Popen([str(exe), *args], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
    timer = threading.Timer(limit, p.kill)
    timer.start()
    q = 0
    try:
        p.stdin.write(f"{n}\n")
        p.stdin.flush()
        while True:
            line = p.stdout.readline()
            if not line:
                return "HANG"
            f = line.split()
            if f[0] == "?":
                q += 1
                p.stdin.write(f"{1 if secret <= int(f[1]) else 0} {q}\n")
                p.stdin.flush()
            else:
                return "OK" if int(f[1]) == secret else "WRONG"
    except (BrokenPipeError, ValueError):
        return "HANG"
    finally:
        timer.cancel()
        p.kill()


def test_interactive(guess, quick):
    r = random.Random(3)
    secrets = [1, 2, 10**9, 10**9 - 1, 777777777] + [r.randint(1, 10**9) for _ in range(5 if quick else 30)]
    bad = sum(interact(guess, ["1", "1"], s) != "OK" for s in secrets)
    report("interactive: scan=true out=true completes", len(secrets), bad, "(real interactor, ~30 queries each)")
    if quick:
        return
    # Negative controls: every other combination must deadlock, otherwise the flag would be pointless.
    bad = 0
    for combo in (["0", "1"], ["1", "0"], ["0", "0"]):
        res = interact(guess, combo, 123456789, limit=1.5)
        if res != "HANG":
            bad += 1
            print(f"  combo scan={combo[0]} out={combo[1]} unexpectedly finished: {res}")
    report("interactive: mismatched modes deadlock", 3, bad, "(negative controls)")


def test_timing(bins, reps):
    r = random.Random(4)
    n = q = 200_000
    a, b, queries = instance(r, n, q, big=True)
    data = " ".join(tokens_of(n, a, b, queries)).encode()
    outs, rows = set(), []
    for exe in bins:
        ts = []
        for _ in range(reps):
            t = time.perf_counter()
            rc, out = run(exe, data)
            ts.append((time.perf_counter() - t) * 1000)
            assert rc == 0
        outs.add(out)
        rows.append((exe.name, statistics.median(ts), min(ts)))
    bad = 0 if len(outs) == 1 else 1
    adj = graph(n, a, b)
    got = list(map(int, next(iter(outs)).split()))
    for i in r.sample(range(q), 3):
        l, rr = queries[i]
        if dijkstra(adj, l - 1)[rr - 1] != got[i]:
            bad += 1
            print(f"  wrong answer on query {i}")
    for name, med, mn in rows:
        print(f"[timing] {name:<14} n=q=2e5: median {med:6.1f} ms  min {mn:6.1f} ms  ({reps} runs, {len(data)/1e6:.1f} MB input)")
    report("n=q=2e5: both modes agree, 3 sampled queries verified", 1 + 3, bad)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="fewer cases, skip negative controls")
    args = ap.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "src" / "bin").mkdir(parents=True)
        (tmp / "Cargo.toml").write_text(CARGO)
        (tmp / "src" / "bin" / "solve_batch.rs").write_text(TEMPLATE)
        (tmp / "src" / "bin" / "solve_inter.rs").write_text(TEMPLATE.replace(FLAG, FLAG.replace("false", "true")))
        (tmp / "src" / "bin" / "echo.rs").write_text(ECHO)
        (tmp / "src" / "bin" / "guess.rs").write_text(GUESS)
        print("building template variants (cargo build --release)...", flush=True)
        b = subprocess.run(["cargo", "build", "--release", "-q", "--manifest-path", str(tmp / "Cargo.toml")], capture_output=True, text=True)
        if b.returncode:
            sys.exit(b.stderr)
        rel = tmp / "target" / "release"
        bins = [rel / "solve_batch", rel / "solve_inter"]
        test_solution(bins, args.quick)
        test_scanner(rel / "echo", args.quick)
        test_interactive(rel / "guess", args.quick)
        test_timing(bins, 5 if args.quick else 15)
    print("\nFAILED" if failures else "\nall template stress tests passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
