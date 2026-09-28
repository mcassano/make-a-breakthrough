# Distinct distances in the grid

How many points of the `n x n` grid can be chosen so that all pairwise
distances are different? Erdős and Guy asked this in 1970. The answers are
OEIS [A271490] (the maximum for each `n`) and [A193838] (the smallest grid
that holds `k` such points).

| source | known before this work |
|---|---|
| A271490 | exact for `n <= 17`: 1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 10, 11, 11, 12, 13, 13 |
| A193838 | exact for `k <= 14` (a(14) = 18, Cariboni 2022); upper bounds a(15..18) <= 20, 22, 24, 27 |

So 14 points fit in the `18 x 18` grid, but whether **15** do was open, and
A271490(18) was not determined.

[A271490]: https://oeis.org/A271490
[A193838]: https://oeis.org/A193838

## Method

`src/distinct.rs` is an exhaustive depth-first search with bitsets.

* **State.** The chosen points, the set of squared distances used so far, and
  the set of forbidden cells. A cell is forbidden once its distance to some
  chosen point is used, or once it is equidistant from two chosen points
  (then two of its distances coincide). Adding a point ORs precomputed
  bitsets: the "ring" of cells at distance `d` from a point for every newly
  used `d`, and the perpendicular bisector of the new point and each old one.
  Cells are considered in a fixed order and each point is taken after the
  previous one, so every set is generated once.
* **Normalisation.** Translate the set so it touches row 0 and column 0.
* **Spanning.** If no `k`-set fits in `(n-1) x (n-1)`, a `k`-set in `n x n`
  has width `n` or height `n`, and after transposing, width `n`. So it
  suffices to search sets touching both column 0 and column `n - 1`.
* **Box symmetry** (`--boxed`). For each bounding-box height `h <= n`, rows
  are searched outside-in: `0, h-1, 1, h-2, ...`. Once rows `0` and `h-1` are
  decided, the pair of their patterns must be lexicographically least among
  its images under the two reflections of the `n x h` box. Every orbit keeps
  at least one member, and most subtrees are cut right after the first two
  rows. Long distances, which a set this large must mostly use, are also
  settled early.
* **Parallelism.** Subtrees rooted at the first two points are independent
  work items, shared by the worker threads.

Solutions are re-checked by an independent pair-by-pair test in Rust and by
`scripts/verify.py`, which shares no code with the search.

## Validation

The unit tests check the search against brute force on small grids. They
check that the spanning and box reductions are exact, i.e.
`exists(n, k) = exists(n-1, k) or exists_reduced(n, k)` for all `n <= 9` and
all relevant `k`, and they reproduce A271490 for `n <= 10`. Larger known
values were reproduced with the full program (`results/distinct/`):

| n | k | expected (A271490) | result |
|---|---|---|---|
| 12 | 11 | impossible | none, 5.5e8 nodes |
| 13 | 12 | impossible | none, 1.1e9 nodes |
| 14 | 12 | impossible | none, 8.0e9 nodes |
| 15 | 12 | possible | found |
| 15 | 13 | impossible | none (see `results/distinct`) |
| 16 | 13 | possible | found |
| 16 | 14 | impossible | see `results/distinct` |
| 17 | 14 | impossible (Cariboni 2022) | see `results/distinct` |

## Reproducing

```sh
cargo build --release
./target/release/distinct 14 12 --boxed    # ~2 minutes on 10 cores: NONE
./target/release/distinct 16 13 --boxed    # finds 13 points in 16 x 16
./target/release/distinct 18 15 --boxed    # the open case (about a day)
scripts/verify.py results                  # re-check every set in results/
```
