# make-a-breakthrough

Exact computer searches on classical problems about points in the `n x n`
grid.

## Result: A271490(18) = 14

**No 15 points of the 18 x 18 grid have pairwise distinct distances.**

Erdős and Guy (1970) asked how many points of the `n x n` grid can have all
pairwise distances different; the answers are OEIS [A271490]. They were known
for `n <= 17` (the last value found by F. A. C. Cariboni in 2022), and it was
known that 14 points fit in the `18 x 18` grid ([A193838], a(14) = 18). An
exhaustive search shows that 15 points do not, so

* **A271490(18) = 14** (new term), and
* **A193838(15) >= 19** (15 points need at least a `19 x 19` grid; the known
  upper bound is 20).

```
$ ./target/release/distinct 18 15 --onion
n=18 k=15: NONE (1.985e12 nodes, 44294.2s, 10 threads, ...)
```

The search took 12.3 hours on a 10-core Apple M2 Pro (log:
`results/distinct/chain.log`). It is in `src/distinct.rs`, and
`docs/distinct-distances.md` explains why each reduction is exact. In short:
bitset forward checking; normalisation by translation; restriction to sets
spanning the full width, which is valid because no 15-set fits in
`17 x 17`; and an "onion" order that searches the boundary of each possible
bounding box first so its symmetries can be broken early.

[A271490]: https://oeis.org/A271490
[A193838]: https://oeis.org/A193838

### Why this can be trusted

* The same program reproduced every known value it was run on, in both
  directions:

  | n | k | known (A271490) | this program |
  |---|---|---|---|
  | 2..10 | a(n), a(n)+1 | | reproduced (unit tests) |
  | 12, 13, 14 | 11, 12, 12 | impossible | impossible |
  | 15 | 12 / 13 | possible / impossible | found / impossible |
  | 16 | 13 / 14 | possible / impossible | found / impossible |
  | 17 | 14 | impossible (Cariboni 2022) | impossible (4.3e11 nodes, 2.4 h) |

* The symmetry reductions lose nothing: on 13 x 13 and 14 x 14 the unreduced
  search and the reduced search find exactly the same classes of solutions up
  to symmetry (63 and 61,331). Unit tests check this and the exactness of
  every reduction on all small grids.
* Every configuration the program reports is re-checked by
  `scripts/verify.py`, an independent verifier sharing no code with the
  search.

## Other problems explored

### No-three-in-line (`docs/no-three-in-line.md`)

Place `2n` points in the `n x n` grid with no three on a line. The
smallest open case is `n = 75`. The repository contains a symmetric search
built on a reformulation: a quarter-turn symmetric solution is a 2-regular
multigraph on the row pairs. It includes a bitset DFS, a large-neighbourhood
search, a SAT encoding and a CaDiCaL external propagator, all verified
against the complete counts in Flammenkamp's database. The DFS needs about
17x less CPU time than a recently published CP-SAT approach at `n = 44`. The
measurements also show that `n = 75`/`78` is still about three orders of
magnitude beyond one workstation, and the doc records which ideas did and did
not shrink the search.

### No-four-on-a-circle (`src/circle.rs`)

The Erdős–Purdy problem: no four points on a common circle or line. The tools
are an exact enumeration of all circles through at least four grid points
(integer circle keys), a SAT encoding and an iterated local search. The local
search reaches the AlphaEvolve lower bounds of Ghosal et al. (2026) for
`n <= 12` but not beyond. Neither SAT nor local search settled a(11) of
[A395821]; no new result here.

[A395821]: https://oeis.org/A395821

## Building and running

Requires Rust 1.70+. The CDCL tool additionally needs CaDiCaL >= 2.0
(`brew install cadical`).

```sh
cargo build --release
cargo test --release

# Distinct distances
./target/release/distinct 14 12 --onion      # 12 points cannot fit in 14 x 14
./target/release/distinct 16 13 --onion      # finds 13 points in 16 x 16
scripts/run_distinct_chain.sh                # the runs recorded in results/

# No-three-in-line
./target/release/no3 dfs 40                  # a quarter-turn symmetric solution
./target/release/no3 count-rot4 22           # exhaustive count
./target/release/no3 verify FILE             # check configurations

# No-four-on-a-circle
./target/release/no4circle search 12         # local search
./target/release/no4circle check-file FILE   # check point sets

scripts/verify.py results                    # independent re-check
```

## Layout

| path | contents |
|---|---|
| `src/distinct.rs`, `src/bin/distinct.rs` | distinct-distance search |
| `src/geometry.rs`, `symmetry.rs`, `cnf.rs`, `dfs.rs`, `lns.rs`, `local.rs`, `format.rs`, `src/bin/no3.rs` | no-three-in-line tools |
| `cdcl/` | CaDiCaL external-propagator solver (C++) |
| `src/circle.rs`, `circle_ls.rs`, `src/bin/no4circle.rs` | no-four-on-a-circle tools |
| `results/` | logs of the exhaustive runs |
| `scripts/` | run, benchmark and verification scripts |
| `docs/` | method notes and measurements |

## License

MIT
