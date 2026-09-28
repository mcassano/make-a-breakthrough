#!/usr/bin/env python3
"""Measure time-to-first-solution of cdcl/no3_cdcl over several seeds.

Example: scripts/bench_cdcl.py --n 28 32 --seeds 10 -- --own-decisions
"""
import argparse
import concurrent.futures as cf
import os
import statistics
import subprocess
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "cdcl", "no3_cdcl")


def run(job):
    n, cls, seed, timeout, extra = job
    start = time.time()
    res = subprocess.run([BIN, str(n), cls, "--seed", str(seed), "--time-limit", str(timeout), *extra],
                         capture_output=True, text=True)
    elapsed = time.time() - start
    return n, seed, elapsed if res.returncode == 10 else float("inf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", required=True)
    ap.add_argument("--cls", default="rot4")
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--timeout", type=float, default=600)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("extra", nargs="*")
    args = ap.parse_args()
    jobs = [(n, args.cls, s, args.timeout, args.extra) for n in args.n for s in range(1, args.seeds + 1)]
    times = {n: [] for n in args.n}
    with cf.ThreadPoolExecutor(args.jobs) as ex:
        for n, seed, t in ex.map(run, jobs):
            times[n].append(t)
    for n in args.n:
        ts = sorted(times[n])
        print(f"n={n} median {statistics.median(ts):.1f}s runs {[round(t, 1) for t in ts]}", flush=True)


if __name__ == "__main__":
    main()
