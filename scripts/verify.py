#!/usr/bin/env python3
"""Independent verifier for the configurations published in this repository.

Deliberately written from scratch in Python (no shared code with the Rust
search tools) so that results do not depend on the correctness of the
searches that produced them.

Usage:
    scripts/verify.py distinct N x1,y1 x2,y2 ...    distinct pairwise distances
    scripts/verify.py no3 CODE ...                   no-three-in-line (Flammenkamp code)
    scripts/verify.py results                        everything under results/
"""
import itertools
import math
import os
import re
import sys

ALPHABET = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz#$%&@?!()[]<>{}=*+|-/~^_:;,."


def check_distinct(n, pts):
    """All points in [0, n)^2, pairwise distinct squared distances."""
    if len(set(pts)) != len(pts):
        return "repeated point"
    for x, y in pts:
        if not (0 <= x < n and 0 <= y < n):
            return f"({x},{y}) outside the {n}x{n} grid"
    seen = {}
    for p, q in itertools.combinations(pts, 2):
        d = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2
        if d in seen:
            return f"squared distance {d} repeated: {seen[d]} and {(p, q)}"
        seen[d] = (p, q)
    return None


def check_no3(code):
    """Decode a Flammenkamp configuration code and check 2n points, no three
    collinear (via cross products of all triples through each pair)."""
    body = code[1:]
    n = len(body) // 2
    pts = [(i // 2, ALPHABET.index(ch)) for i, ch in enumerate(body)]
    if len(set(pts)) != 2 * n or any(c >= n for _, c in pts):
        return n, "bad point set"
    lines = set()
    for (r1, c1), (r2, c2) in itertools.combinations(pts, 2):
        a, b = r2 - r1, c1 - c2
        c = -(a * c1 + b * r1)
        g = math.gcd(math.gcd(a, b), c)
        a, b, c = a // g, b // g, c // g
        if a < 0 or (a == 0 and b < 0):
            a, b, c = -a, -b, -c
        if (a, b, c) in lines:
            return n, "three points on a line"
        lines.add((a, b, c))
    return n, None


def parse_points(tokens):
    return [tuple(int(v) for v in t.split(",")) for t in tokens]


def verify_results(root):
    bad = 0
    pattern = re.compile(r"n=(\d+) k=(\d+): FOUND ((?:\(\d+,\d+\) ?)+)")
    for dirpath, _, files in os.walk(root):
        for name in sorted(files):
            path = os.path.join(dirpath, name)
            with open(path) as f:
                for line in f:
                    m = pattern.search(line)
                    if m:
                        n, k = int(m.group(1)), int(m.group(2))
                        pts = [tuple(map(int, p.strip("()").split(","))) for p in m.group(3).split()]
                        err = check_distinct(n, pts) or (None if len(pts) == k else "wrong size")
                        print(f"{path}: distinct n={n} k={k}: {'OK' if not err else 'FAIL ' + err}")
                        bad += bool(err)
    return bad


def main(argv):
    if len(argv) >= 3 and argv[1] == "distinct":
        err = check_distinct(int(argv[2]), parse_points(argv[3:]))
        print("OK" if not err else f"FAIL: {err}")
        return 1 if err else 0
    if len(argv) >= 3 and argv[1] == "no3":
        bad = 0
        for code in argv[2:]:
            n, err = check_no3(code)
            print(f"n={n}: {'OK' if not err else 'FAIL ' + err}")
            bad += bool(err)
        return 1 if bad else 0
    if len(argv) == 2 and argv[1] == "results":
        root = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
        return 1 if verify_results(root) else 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
