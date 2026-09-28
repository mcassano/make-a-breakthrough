# Distinct distances in the grid

How many points of the `n x n` grid can be chosen so that all pairwise
distances are different? Erdős and Guy asked this in 1970. The answers are
OEIS [A271490] (the maximum for each `n`) and [A193838] (the smallest grid
that holds `k` such points).

| sequence | known before this work |
|---|---|
| A271490 | exact for `n <= 17`: 1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 10, 11, 11, 12, 13, 13 |
| A193838 | exact for `k <= 14` (a(14) = 18, F. A. C. Cariboni 2022); upper bounds a(15..18) <= 20, 22, 24, 27 |

So 14 points fit in the `18 x 18` grid, but whether **15** do was open, and
A271490(18) was not determined.

[A271490]: https://oeis.org/A271490
[A193838]: https://oeis.org/A193838

## Method

`src/distinct.rs` is an exhaustive depth-first search with bitsets.

**State and propagation.** The search keeps the chosen points, the set of
squared distances used so far, and the set of forbidden cells. A cell is
forbidden once its distance to some chosen point is used, or once it is
equidistant from two chosen points (then two of its distances coincide).
Adding a point ORs precomputed bitsets: the "ring" of cells at distance `d`
from a point for every newly used `d`, and the perpendicular bisector of the
new point and each old one. Cells are visited in a fixed order and each point
comes after the previous one, so every set is generated once. A counting bound
also prunes: the pairs still to be formed need that many distinct unused
distances, and those can only occur between a chosen point and a later cell or
between two later cells (precomputed per position).

**Reductions.** All of these are exact:

1. *Translation.* Every set can be moved inside the grid so that it touches
   row 0 and column 0.
2. *Spanning.* If no `k`-set fits in `(n-1) x (n-1)`, then a `k`-set in
   `n x n` spans the full width or the full height, and after transposing,
   the full width. So it suffices to search sets whose bounding box is
   `n x h` for some `h <= n`.
3. *Onion order and box symmetry.* For each height `h`, cells are searched
   ring by ring from the boundary of the `n x h` box inwards. A set with this
   exact bounding box touches all four sides, so its first point lies on the
   boundary ring, and once the ring is decided all four sides must be hit.
   The chosen boundary cells must then be lexicographically least among their
   images under the symmetries of the box: the two reflections (and their
   product), plus the transposition and the other square symmetries when
   `h = n`. Every symmetry class keeps at least one member. Searching the
   boundary first also decides the long distances, which a large set must
   mostly use, before the interior.

On the validation cases the onion order needs 5-8x fewer nodes than row order
with mirror breaking (`13 x 13`: 2.1e8 vs 1.7e9; `14 x 14`: 2.3e9 vs 1.2e10;
`15 x 15`: 1.0e10 vs 6.1e10).

**Parallelism.** Subtrees rooted at the first two points are independent
work items shared by the worker threads.

## Correctness

* Unit tests compare the search with brute force on small grids. They check
  that every reduction is exact, i.e.
  `exists(n, k) = exists(n-1, k) or exists_reduced(n, k)` for `n <= 9` and
  all relevant `k`, in all three orders. They also reproduce A271490 for
  `n <= 10`.
* Larger known values reproduced with the full program (logs in
  `results/distinct/`):

  | n | k | A271490 says | result | nodes | time (10 threads) |
  |---|---|---|---|---|---|
  | 12 | 11 | impossible | none | 5.5e8 (row order) | 19 s |
  | 13 | 12 | impossible | none | 2.1e8 | 9 s |
  | 14 | 12 | impossible | none | 2.3e9 | 79 s |
  | 15 | 12 | possible | found | | |
  | 15 | 13 | impossible | none | 1.0e10 | 200 s |
  | 16 | 13 | possible | found | | |
  | 16 | 14 | impossible | none | 4.8e10 | 1,106 s |
  | 17 | 14 | impossible (Cariboni 2022) | see `results/distinct/chain.log` | | |

* Every set the program reports is re-checked pair by pair in Rust and by
  `scripts/verify.py`, which shares no code with the search.

## Reproducing

```sh
cargo build --release
./target/release/distinct 14 12 --onion     # about a minute on 10 cores: NONE
./target/release/distinct 16 13 --onion     # finds 13 points in 16 x 16
./target/release/distinct 18 15 --onion     # the open case
scripts/verify.py results                   # re-check every set in results/
```

`--boxed` (rows `0, h-1, 1, h-2, ...`) and `--mirror` (row order with the
spanning and mirror reductions) select the older, slower orders, which serve
as cross-checks.
