# No-three-in-line: search tools and measurements

The no-three-in-line problem asks for `2n` points in the `n x n` grid with no
three on a line. As of September 2026, solutions are known for every
`n <= 74` and for `n = 76`; `n = 75` is the smallest open case and `n = 78`
would be a new record ([Flammenkamp's database][db]). All large solutions
have quarter-turn symmetry (`rot4`, even `n`) or quarter-turn symmetry except
for one diagonal pair (`rct4`, odd `n`).

This directory of the repository contains the tools I built to attack that
frontier, and what I learned. **It does not contain a new record**: the
measurements below indicate that `n = 75` or `78` needs roughly three orders
of magnitude more work than one workstation can do with these methods.

[db]: https://wwwhomes.uni-bielefeld.de/achim/no3in/readme.html

## The orbit-graph formulation

Put the origin at the grid center and number the `m = n/2` pairs of mirror
rows `±y` by `1..m`. A quarter-turn orbit `{(x,y), (-y,x), (-x,-y), (y,-x)}`
with `x, y > 0` puts exactly one point in row pair `y` and one in row pair
`x`. So a `rot4` configuration with two points per row is exactly a
**2-regular multigraph on the `m` row pairs with oriented edges** (an edge
`(x, y)` is the orbit of the quadrant cell `(x, y)`), subject to collinearity
constraints. Columns take care of themselves by symmetry. The `rct4` class has
the same structure plus one "axis" orbit filling the middle row and one
diagonal pair (see `src/local.rs`).

This view gives:

* a local search whose moves (swapping two coordinates) never break the
  row/column counts (`src/local.rs`);
* a depth-first search that branches on row pairs, with the blocked orbits
  kept as one 64-bit mask per quadrant row and column (`src/dfs.rs`). Because
  the chosen set is invariant under the quarter turn, adding an orbit only
  requires walking the lines through *one* of its points.

## Correctness checks

* `no3 verify` independently checks every one of the 431,008 configurations
  in Flammenkamp's database (all valid).
* The exhaustive `rot4` count reproduces brute force for `n <= 12` and the
  complete database counts for `n = 14..22` (unit tests).
* The CaDiCaL-based solver (`cdcl/`) enumerates the same `rot4` counts and the
  database's `rct4` counts for `n = 17..23`.
* Every reported solution is re-verified by `scripts/verify.py`, which shares
  no code with the search tools.

## Measurements (Apple M2 Pro, 10 cores)

Time to first `rot4` solution:

| method | n = 32 | n = 36 | n = 40 | n = 44 |
|---|---|---|---|---|
| randomized DFS (`no3 dfs`, 10 threads) | ~0.1 s | ~1 s | 20–40 s | ~15 min (7e9 nodes) |
| large-neighbourhood search (`no3 lns`, 10 threads) | 0.4 s | ~2 s | ~40 s | ~7 min |
| kissat on the orbit CNF (`no3 cnf`), per run | – | median 179 s | > 35 min | – |
| CaDiCaL + external propagator (`cdcl/`), per run | 17–30 s | – | – | – |
| CP-SAT, mean single run ([arXiv:2602.07751][cpsat]) | – | – | – | ~40 h |

[cpsat]: https://arxiv.org/abs/2602.07751

At `n = 44` the DFS uses about 8,500 CPU-seconds against the ~145,000
CPU-seconds reported for a single CP-SAT run, roughly 17 times less.

### What changed the search tree, and what did not

* **Branching order matters a lot.** Exhaustive enumeration at `n = 28`:
  innermost row pairs first takes 11 s; outermost first takes 479 s (40x).
  Choosing the row pair with the fewest usable orbits per missing orbit is best
  (10.0M nodes vs 13.3M for plain fewest-remaining-options).
* **Solution density per node is invariant.** I sampled 15,000 random nodes at
  depth 3 for `n = 32`, counted the solutions below each exactly, and split
  them into quintiles by cheap features (usable orbits, minimum slack,
  log-product of options). Every quintile had about one solution per million
  nodes. Restricting the search (forbidding the central cells, or cells near
  the diagonals) removes nodes and solutions in the same proportion. So
  guiding the search by such features, or restricting it heuristically, does
  not help; only stronger inference does.
* **Singleton consistency** (tentatively choosing every usable orbit and
  discarding those that make some row pair unsatisfiable) shrinks the tree by
  5.6x at `n = 22` and 6.5x at `n = 28`. That gain grows too slowly to pay for
  its cost.
* **Clause learning did not shrink the tree** at the sizes I could measure:
  CaDiCaL with a lazy line propagator and the same branching heuristic needs
  about as many decisions per solution (~1e6 at `n = 32`) as the DFS needs
  nodes, at ~30x the cost per step.

### Why the frontier is out of reach here

Nodes per solution grow by a factor of about 2.5–3.3 per step `n -> n + 2`
(exhaustive counts for `n <= 30`, random restarts for `n <= 44`).
Extrapolating from `n = 44` gives 1e15–1e17 nodes for one solution at
`n = 78`, versus about 1e7 nodes/s on this machine. The known record solutions
came from GPU clusters (Riley, Prellberg) and a new SAT solver run by Heule.

## Reproducing

```sh
cargo build --release
./target/release/no3 dfs 40                        # a rot4 solution for n = 40
./target/release/no3 lns 40 --k 8 --budget 100000  # the same with LNS
./target/release/no3 count-rot4 22                 # exhaustive count (16)
./target/release/no3 cnf 40 rot4 n40.cnf           # SAT encoding
make -C cdcl && cdcl/no3_cdcl 30 rot4 --own-decisions
```
