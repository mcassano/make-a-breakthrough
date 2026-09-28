#!/usr/bin/env python3
"""Benchmark SAT solvers on no-three-in-line encodings.

Runs every (n, seed) job with a timeout, several jobs in parallel, and prints
the solve time of each run and the median per n. Satisfying assignments are
decoded and verified with `no3 extract`.

Example:
    scripts/bench_sat.py --n 36 40 44 --seeds 5 --solver kissat -- --direct 6
"""
import argparse
import concurrent.futures as cf
import os
import statistics
import subprocess
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NO3 = os.path.join(ROOT, "target", "release", "no3")


def solver_cmd(solver, seed, cnf):
    if solver == "kissat":
        return ["kissat", "-q", f"--seed={seed}", cnf]
    if solver == "cadical":
        return ["cadical", "-q", f"--seed={seed}", cnf]
    raise SystemExit(f"unknown solver {solver}")


def run(job):
    n, cls, seed, solver, timeout, cnf, workdir = job
    out = os.path.join(workdir, f"{solver}-{cls}-{n}-{seed}.out")
    start = time.time()
    try:
        with open(out, "w") as f:
            subprocess.run(solver_cmd(solver, seed, cnf), stdout=f, timeout=timeout)
    except subprocess.TimeoutExpired:
        return n, seed, None, None
    elapsed = time.time() - start
    res = subprocess.run([NO3, "extract", str(n), cls, out], capture_output=True, text=True)
    if res.returncode != 0:
        return n, seed, elapsed, "UNSAT/ERROR " + res.stderr.strip()
    return n, seed, elapsed, res.stdout.strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", required=True)
    ap.add_argument("--cls", default="rot4")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--solver", default="kissat")
    ap.add_argument("--timeout", type=float, default=600)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("cnf_args", nargs="*", help="extra arguments for `no3 cnf`")
    args = ap.parse_args()

    workdir = tempfile.mkdtemp(prefix="no3bench-")
    jobs = []
    for n in args.n:
        cnf = os.path.join(workdir, f"{args.cls}-{n}.cnf")
        subprocess.run([NO3, "cnf", str(n), args.cls, cnf, *args.cnf_args], check=True, capture_output=True)
        jobs += [(n, args.cls, s, args.solver, args.timeout, cnf, workdir) for s in range(1, args.seeds + 1)]

    times = {n: [] for n in args.n}
    with cf.ThreadPoolExecutor(args.jobs) as ex:
        for n, seed, t, code in ex.map(run, jobs):
            times[n].append(t if t is not None else float("inf"))
            shown = "timeout" if t is None else f"{t:8.1f}s"
            print(f"n={n} seed={seed} {shown} {code or ''}", flush=True)
    for n in args.n:
        print(f"n={n} median {statistics.median(times[n]):.1f}s  runs {sorted(round(t, 1) for t in times[n])}")


if __name__ == "__main__":
    main()
