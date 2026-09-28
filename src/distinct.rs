//! Exhaustive search for point sets with pairwise distinct distances in the
//! `n x n` grid (Erdős and Guy, "Distinct distances between lattice points",
//! 1970; OEIS A271490 and A193838).
//!
//! Cells are visited in a fixed order and each chosen point comes after the
//! previous one, so every set is generated once. The search keeps
//!
//! * the set of squared distances already used, and
//! * the set of cells that can no longer be chosen: a cell is forbidden if
//!   its distance to some chosen point is already used, or if it is
//!   equidistant from two chosen points (which would repeat a distance).
//!
//! Both are bitsets, and adding a point updates the forbidden cells by OR-ing
//! precomputed bitsets: the "ring" of cells at a given distance from a point
//! and the perpendicular bisector of two points.
//!
//! Reductions (see `docs/distinct-distances.md` for the arguments):
//!
//! * Every set can be translated to touch row 0 and column 0.
//! * `span`: if no `k`-set fits in the `(n-1) x (n-1)` grid, every `k`-set of
//!   the `n x n` grid spans its full width or height; after transposing, the
//!   full width. So only sets touching column `n - 1` are needed.
//! * `mirror`: one set of each mirror pair under `x -> n - 1 - x`.
//! * [`Order::BoxedRows`] and [`Order::Onion`] search each bounding-box
//!   height `h` separately, so the box `n x h` is known. Its sides are
//!   searched first (rows `0` and `h - 1`, or the whole boundary ring), and a
//!   set is kept only if that first part is lexicographically least among
//!   its images under the symmetries of the box: 4 reflections, or all 8
//!   symmetries of the square when `h = n` in onion order.

use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

const CW: usize = 8; // cell bitset words: up to 512 cells (n <= 22)
const DW: usize = 16; // distance bitset words: squared distances < 1024

type Cells = [u64; CW];
type Dists = [u64; DW];

#[inline]
fn get(b: &[u64], i: usize) -> bool {
    b[i >> 6] >> (i & 63) & 1 == 1
}

#[inline]
fn set(b: &mut [u64], i: usize) {
    b[i >> 6] |= 1 << (i & 63);
}

#[inline]
fn or_into(dst: &mut Cells, src: &Cells) {
    for w in 0..CW {
        dst[w] |= src[w];
    }
}

/// The order in which cells are searched.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Order {
    /// Row-major over the whole `n x n` grid.
    #[default]
    Rows,
    /// Per bounding-box height `h`: rows `0, h-1, 1, h-2, ...`.
    BoxedRows,
    /// Per bounding-box height `h`: the boundary ring of the `n x h` box
    /// first, then the next ring inwards, and so on.
    Onion,
}

const LEFT: u8 = 1;
const RIGHT: u8 = 2;
const BOTTOM: u8 = 4;
const TOP: u8 = 8;

/// Precomputed geometry of an `n x h` box whose cells are searched in a
/// fixed order; `xy[i]` is the cell at position `i`.
pub struct Grid {
    pub n: usize,
    pub h: usize,
    order: Order,
    xy: Vec<(usize, usize)>,
    /// Position of cell `(x, y)` at index `y * n + x`.
    pos: Vec<usize>,
    maxd: usize,
    /// `ring[p * maxd + d]`: positions at squared distance `d` from `p`.
    ring: Vec<Cells>,
    /// `bisector[p * cells + q]`: positions equidistant from `p` and `q`.
    bisector: Vec<Cells>,
    /// `reach[p * (cells + 1) + f]`: distances from `p` to positions `>= f`.
    reach: Vec<Dists>,
    /// `among[f]`: distances between two positions that are both `>= f`.
    among: Vec<Dists>,
    /// Positions on the left, right, bottom and top side of the box.
    sides: [Cells; 4],
    /// The first point of every searched set lies before `first_end`.
    first_end: usize,
    /// Positions before `sym_end` decide the symmetry reduction.
    sym_end: usize,
}

impl Grid {
    /// The `n x n` grid in row-major order.
    pub fn new(n: usize) -> Grid {
        let xy = (0..n * n).map(|i| (i % n, i / n)).collect();
        Grid::build(n, n, Order::Rows, xy, n, n)
    }

    /// The `n x h` box with rows ordered `0, h-1, 1, h-2, ...`.
    pub fn boxed_rows(n: usize, h: usize) -> Grid {
        let mut ys = Vec::with_capacity(h);
        let (mut lo, mut hi) = (0, h);
        while lo < hi {
            ys.push(lo);
            lo += 1;
            if lo < hi {
                hi -= 1;
                ys.push(hi);
            }
        }
        let xy = ys.iter().flat_map(|&y| (0..n).map(move |x| (x, y))).collect();
        let sym_end = if h >= 2 { 2 * n } else { n };
        Grid::build(n, h, Order::BoxedRows, xy, n, sym_end)
    }

    /// The `n x h` box ordered by rings from the boundary inwards (row-major
    /// within a ring).
    pub fn onion(n: usize, h: usize) -> Grid {
        let ring = |&(x, y): &(usize, usize)| x.min(n - 1 - x).min(y).min(h - 1 - y);
        let mut xy: Vec<(usize, usize)> = (0..h).flat_map(|y| (0..n).map(move |x| (x, y))).collect();
        xy.sort_by_key(|c| (ring(c), c.1, c.0));
        let boundary = xy.iter().filter(|c| ring(c) == 0).count();
        Grid::build(n, h, Order::Onion, xy, boundary, boundary)
    }

    fn build(n: usize, h: usize, order: Order, xy: Vec<(usize, usize)>, first_end: usize, sym_end: usize) -> Grid {
        let cells = xy.len();
        assert!(
            n >= 2 && h >= 1 && cells == n * h && cells <= CW * 64,
            "grid size out of range"
        );
        let maxd = (n - 1).pow(2) + (h - 1).pow(2) + 1;
        assert!(maxd <= DW * 64);
        let dist = |p: usize, q: usize| {
            let ((px, py), (qx, qy)) = (xy[p], xy[q]);
            px.abs_diff(qx).pow(2) + py.abs_diff(qy).pow(2)
        };
        let mut pos = vec![0; cells];
        for (i, &(x, y)) in xy.iter().enumerate() {
            pos[y * n + x] = i;
        }
        let mut ring = vec![[0u64; CW]; cells * maxd];
        let mut bisector = vec![[0u64; CW]; cells * cells];
        for p in 0..cells {
            for c in 0..cells {
                set(&mut ring[p * maxd + dist(p, c)], c);
            }
            for q in 0..cells {
                if p != q {
                    for c in 0..cells {
                        if dist(p, c) == dist(q, c) {
                            set(&mut bisector[p * cells + q], c);
                        }
                    }
                }
            }
        }
        let mut reach = vec![[0u64; DW]; cells * (cells + 1)];
        let mut among = vec![[0u64; DW]; cells + 1];
        for f in (0..cells).rev() {
            for p in 0..cells {
                let mut r = reach[p * (cells + 1) + f + 1];
                if p != f {
                    set(&mut r, dist(p, f));
                }
                reach[p * (cells + 1) + f] = r;
            }
            let mut a = among[f + 1];
            for w in 0..DW {
                a[w] |= reach[f * (cells + 1) + f + 1][w];
            }
            among[f] = a;
        }
        let mut sides = [[0u64; CW]; 4];
        for (i, &(x, y)) in xy.iter().enumerate() {
            for (s, on) in [x == 0, x == n - 1, y == 0, y == h - 1].into_iter().enumerate() {
                if on {
                    set(&mut sides[s], i);
                }
            }
        }
        Grid {
            n,
            h,
            order,
            xy,
            pos,
            maxd,
            ring,
            bisector,
            reach,
            among,
            sides,
            first_end,
            sym_end,
        }
    }

    fn cells(&self) -> usize {
        self.xy.len()
    }

    /// Squared distance between two positions.
    pub fn dist(&self, p: usize, q: usize) -> usize {
        let ((px, py), (qx, qy)) = (self.xy[p], self.xy[q]);
        px.abs_diff(qx).pow(2) + py.abs_diff(qy).pow(2)
    }

    /// Side flags of a position.
    fn side_bits(&self, p: usize) -> u8 {
        let (x, y) = self.xy[p];
        ((x == 0) as u8 * LEFT)
            | ((x == self.n - 1) as u8 * RIGHT)
            | ((y == 0) as u8 * BOTTOM)
            | ((y == self.h - 1) as u8 * TOP)
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Options {
    /// Only sets touching column `n - 1` (see the module documentation).
    pub span: bool,
    /// Symmetry reduction (requires `span` in row order).
    pub mirror: bool,
    /// Stop at the first solution.
    pub first_only: bool,
    pub order: Order,
}

#[derive(Clone, Copy)]
struct Node {
    forbidden: Cells,
    used: Dists,
}

pub struct Search<'a> {
    g: &'a Grid,
    k: usize,
    opts: Options,
    /// Sides of the box every set must touch.
    required: u8,
    chosen: Vec<usize>,
    pub nodes: u64,
    pub solutions: u64,
    pub found: Option<Vec<usize>>,
}

impl<'a> Search<'a> {
    pub fn new(g: &'a Grid, k: usize, opts: Options) -> Search<'a> {
        let required = match g.order {
            Order::Rows => {
                assert!(!opts.mirror || opts.span, "mirror breaking needs span");
                LEFT | BOTTOM | if opts.span { RIGHT } else { 0 }
            }
            Order::BoxedRows | Order::Onion => {
                assert!(opts.span && opts.mirror, "box orders need span and mirror");
                LEFT | RIGHT | BOTTOM | TOP
            }
        };
        Search {
            g,
            k,
            opts,
            required,
            chosen: Vec::new(),
            nodes: 0,
            solutions: 0,
            found: None,
        }
    }

    /// The node after adding position `p` to the chosen points of `node`.
    fn extend(&self, node: &Node, p: usize) -> Node {
        let g = self.g;
        let cells = g.cells();
        let mut next = *node;
        set(&mut next.forbidden, p);
        // Cells at an already used distance from p.
        for (w, &word) in node.used.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let d = w * 64 + bits.trailing_zeros() as usize;
                bits &= bits - 1;
                or_into(&mut next.forbidden, &g.ring[p * g.maxd + d]);
            }
        }
        for &q in &self.chosen {
            // The new distance |pq| may not occur again anywhere.
            let d = g.dist(p, q);
            set(&mut next.used, d);
            or_into(&mut next.forbidden, &g.ring[p * g.maxd + d]);
            for &r in &self.chosen {
                or_into(&mut next.forbidden, &g.ring[r * g.maxd + d]);
            }
            // Cells equidistant from p and q.
            or_into(&mut next.forbidden, &g.bisector[p * cells + q]);
        }
        next
    }

    /// Whether the chosen points before `sym_end` pass the symmetry
    /// reduction (and, in box orders, touch every side of the box).
    fn canonical(&self) -> bool {
        let g = self.g;
        let (n, h) = (g.n, g.h);
        let prefix: Vec<usize> = self.chosen.iter().copied().filter(|&p| p < g.sym_end).collect();
        match g.order {
            Order::Rows => {
                let r0: Vec<usize> = prefix.iter().map(|&p| g.xy[p].0).collect();
                let mut m: Vec<usize> = r0.iter().map(|&x| n - 1 - x).collect();
                m.sort_unstable();
                r0 <= m
            }
            Order::BoxedRows | Order::Onion => {
                // The prefix holds whole sides of the box: rows 0 and h-1 in
                // boxed-row order, the entire boundary in onion order.
                let sides = prefix.iter().fold(0, |acc, &p| acc | g.side_bits(p));
                let prefix_sides = if g.order == Order::Onion {
                    LEFT | RIGHT | BOTTOM | TOP
                } else {
                    BOTTOM | TOP
                };
                if sides & prefix_sides != prefix_sides {
                    return false;
                }
                // Images of the prefix under the symmetries of the box, as
                // sorted position lists. Rows 0 and h-1 (boxed rows) and the
                // boundary ring (onion) are unions of whole symmetry orbits.
                let transposes = g.order == Order::Onion && h == n;
                let maps: &[(bool, bool, bool)] = if transposes {
                    &[
                        (false, true, false),
                        (true, false, false),
                        (true, true, false),
                        (false, false, true),
                        (false, true, true),
                        (true, false, true),
                        (true, true, true),
                    ]
                } else {
                    &[(false, true, false), (true, false, false), (true, true, false)]
                };
                maps.iter().all(|&(fx, fy, t)| {
                    let mut img: Vec<usize> = prefix
                        .iter()
                        .map(|&p| {
                            let (mut x, mut y) = g.xy[p];
                            if t {
                                (x, y) = (y, x);
                            }
                            if fx {
                                x = n - 1 - x;
                            }
                            if fy {
                                y = h - 1 - y;
                            }
                            g.pos[y * n + x]
                        })
                        .collect();
                    img.sort_unstable();
                    prefix <= img
                })
            }
        }
    }

    /// Candidate positions of `node` at index `>= from`.
    fn candidates(&self, node: &Node, from: usize) -> (Cells, usize) {
        let cells = self.g.cells();
        let mut avail = [0u64; CW];
        let mut count = 0;
        for w in 0..CW {
            let lo = w * 64;
            let mut m = !node.forbidden[w];
            if lo + 64 <= from || lo >= cells {
                m = 0;
            } else {
                if lo < from {
                    m &= !0u64 << (from - lo);
                }
                if lo + 64 > cells {
                    m &= (1u64 << (cells - lo)) - 1;
                }
            }
            avail[w] = m;
            count += m.count_ones() as usize;
        }
        (avail, count)
    }

    /// Depth-first search over positions `>= from`; `touched` holds the
    /// sides reached so far. Returns true to stop.
    fn dfs(&mut self, node: &Node, from: usize, touched: u8) -> bool {
        self.nodes += 1;
        let g = self.g;
        let need = self.k - self.chosen.len();
        if need == 0 {
            let ok =
                touched & self.required == self.required && (!self.opts.mirror || from > g.sym_end || self.canonical());
            if ok {
                self.solutions += 1;
                if self.found.is_none() {
                    self.found = Some(self.chosen.clone());
                }
                return self.opts.first_only;
            }
            return false;
        }
        let (avail, mut count) = self.candidates(node, from);
        if count < need {
            return false;
        }
        // Every required side not yet touched must still be reachable.
        let missing = self.required & !touched;
        for s in 0..4 {
            if missing >> s & 1 == 1 && !(0..CW).any(|w| avail[w] & g.sides[s][w] != 0) {
                return false;
            }
        }
        // The remaining pairs need distinct unused distances, realisable only
        // between a chosen point and a later cell or between two later cells.
        let pairs_needed = self.k * (self.k - 1) / 2 - self.chosen.len() * (self.chosen.len() - 1) / 2;
        let cells = g.cells();
        let f = from.min(cells);
        let mut possible = g.among[f];
        for &q in &self.chosen {
            let r = &g.reach[q * (cells + 1) + f];
            for w in 0..DW {
                possible[w] |= r[w];
            }
        }
        let achievable: u32 = (0..DW).map(|w| (possible[w] & !node.used[w]).count_ones()).sum();
        if (achievable as usize) < pairs_needed {
            return false;
        }
        for w in 0..CW {
            let mut bits = avail[w];
            while bits != 0 {
                let p = w * 64 + bits.trailing_zeros() as usize;
                bits &= bits - 1;
                if count < need {
                    return false;
                }
                count -= 1;
                // Leaving the prefix fixes what decides the symmetry
                // reduction; all later siblings share it.
                if self.opts.mirror && p >= g.sym_end && from <= g.sym_end && !self.canonical() {
                    return false;
                }
                let child = self.extend(node, p);
                self.chosen.push(p);
                let stop = self.dfs(&child, p + 1, touched | g.side_bits(p));
                self.chosen.pop();
                if stop {
                    return true;
                }
            }
        }
        false
    }

    /// Searches the subtree whose first point is `a` (which must precede
    /// `first_end`) and whose second point is `b > a` (if `k >= 2`).
    /// Returns true to stop.
    pub fn run_prefix(&mut self, a: usize, b: Option<usize>) -> bool {
        let g = self.g;
        assert!(a < g.first_end);
        let root = Node {
            forbidden: [0; CW],
            used: [0; DW],
        };
        self.chosen.clear();
        let n1 = self.extend(&root, a);
        self.chosen.push(a);
        let touched = g.side_bits(a);
        let stop = match b {
            None => self.dfs(&n1, a + 1, touched),
            Some(b) if !get(&n1.forbidden, b) => {
                if self.opts.mirror && b >= g.sym_end && !self.canonical() {
                    false
                } else {
                    let n2 = self.extend(&n1, b);
                    self.chosen.push(b);
                    let stop = self.dfs(&n2, b + 1, touched | g.side_bits(b));
                    self.chosen.pop();
                    stop
                }
            }
            Some(_) => false,
        };
        self.chosen.clear();
        stop
    }

    /// Knuth's unbiased estimator of the size of the search tree (row order,
    /// ignoring the mirror reduction) from one random root-to-leaf probe.
    pub fn knuth_probe(&mut self, rng: &mut crate::local::Rng) -> f64 {
        let g = self.g;
        let mut node = Node {
            forbidden: [0; CW],
            used: [0; DW],
        };
        self.chosen.clear();
        let (mut estimate, mut weight) = (1.0, 1.0);
        let (mut from, mut touched) = (0, 0u8);
        while self.chosen.len() < self.k {
            let need = self.k - self.chosen.len();
            let (avail, _) = self.candidates(&node, from);
            let mut kids = Vec::new();
            for w in 0..CW {
                let mut bits = avail[w];
                while bits != 0 {
                    let c = w * 64 + bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    if !(self.chosen.is_empty() && c >= g.first_end) {
                        kids.push(c);
                    }
                }
            }
            let missing = self.required & !touched;
            let reachable = (0..4).all(|s| missing >> s & 1 == 0 || kids.iter().any(|&c| get(&g.sides[s], c)));
            if kids.len() < need || !reachable {
                break;
            }
            // dfs() only expands a child while enough candidates remain.
            let usable = kids.len() - (need - 1);
            weight *= usable as f64;
            estimate += weight;
            let p = kids[rng.below(usable)];
            node = self.extend(&node, p);
            self.chosen.push(p);
            from = p + 1;
            touched |= g.side_bits(p);
        }
        self.chosen.clear();
        estimate
    }
}

/// Result of a full search.
pub struct Outcome {
    pub solutions: u64,
    pub nodes: u64,
    /// One solution, as `(x, y)` points.
    pub example: Option<Vec<(usize, usize)>>,
}

/// Searches all (reduced) `k`-point sets of the `n x n` grid with `threads`
/// workers, calling `progress(h, done, total, nodes, solutions)` after each
/// work item (`h` is the box height, or `n` in row order).
pub fn search(
    n: usize,
    k: usize,
    opts: Options,
    threads: usize,
    progress: impl Fn(usize, usize, usize, u64, u64) + Sync,
) -> Outcome {
    let grids: Box<dyn Iterator<Item = Grid>> = match opts.order {
        Order::Rows => Box::new(std::iter::once(Grid::new(n))),
        Order::BoxedRows => Box::new((1..=n).map(move |h| Grid::boxed_rows(n, h))),
        Order::Onion => Box::new((1..=n).map(move |h| Grid::onion(n, h))),
    };
    let mut total = Outcome {
        solutions: 0,
        nodes: 0,
        example: None,
    };
    for g in grids {
        let out = search_grid(&g, k, opts, threads, &|d, t, nodes, sols| {
            progress(g.h, d, t, nodes, sols)
        });
        total.nodes += out.nodes;
        total.solutions += out.solutions;
        if total.example.is_none() {
            total.example = out.example;
        }
        if opts.first_only && total.example.is_some() {
            break;
        }
    }
    total
}

fn search_grid(
    g: &Grid,
    k: usize,
    opts: Options,
    threads: usize,
    progress: &(dyn Fn(usize, usize, u64, u64) + Sync),
) -> Outcome {
    let cells = g.cells();
    let mut items: Vec<(usize, Option<usize>)> = Vec::new();
    for a in 0..g.first_end {
        if k == 1 {
            items.push((a, None));
        } else {
            items.extend((a + 1..cells).map(|b| (a, Some(b))));
        }
    }
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let nodes = AtomicU64::new(0);
    let solutions = AtomicU64::new(0);
    let stop = AtomicBool::new(false);
    let example: Mutex<Option<Vec<usize>>> = Mutex::new(None);
    std::thread::scope(|scope| {
        for _ in 0..threads.max(1) {
            scope.spawn(|| {
                let mut s = Search::new(g, k, opts);
                while !stop.load(Ordering::Relaxed) {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    let Some(&(a, b)) = items.get(i) else { break };
                    let (n0, s0) = (s.nodes, s.solutions);
                    let halt = s.run_prefix(a, b);
                    nodes.fetch_add(s.nodes - n0, Ordering::Relaxed);
                    solutions.fetch_add(s.solutions - s0, Ordering::Relaxed);
                    if let Some(f) = s.found.take() {
                        example.lock().unwrap().get_or_insert(f);
                    }
                    if halt {
                        stop.store(true, Ordering::Relaxed);
                    }
                    let d = done.fetch_add(1, Ordering::Relaxed) + 1;
                    progress(
                        d,
                        items.len(),
                        nodes.load(Ordering::Relaxed),
                        solutions.load(Ordering::Relaxed),
                    );
                }
            });
        }
    });
    let example = example
        .into_inner()
        .unwrap()
        .map(|v| v.iter().map(|&p| g.xy[p]).collect());
    Outcome {
        solutions: solutions.into_inner(),
        nodes: nodes.into_inner(),
        example,
    }
}

/// Independent check that the points are distinct cells of the `n x n` grid
/// with pairwise distinct distances.
pub fn check(n: usize, pts: &[(usize, usize)]) -> Result<(), String> {
    let mut seen = std::collections::HashMap::new();
    for (i, &(x1, y1)) in pts.iter().enumerate() {
        if x1 >= n || y1 >= n {
            return Err(format!("({x1}, {y1}) is outside the {n}x{n} grid"));
        }
        for &(x2, y2) in &pts[i + 1..] {
            let d = (x1 as i64 - x2 as i64).pow(2) + (y1 as i64 - y2 as i64).pow(2);
            if d == 0 {
                return Err(format!("({x1}, {y1}) is repeated"));
            }
            if let Some(prev) = seen.insert(d, ((x1, y1), (x2, y2))) {
                return Err(format!(
                    "squared distance {d} occurs for {prev:?} and {:?}",
                    ((x1, y1), (x2, y2))
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OEIS A271490(n) for n = 1..=17.
    const A271490: [usize; 17] = [1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 10, 11, 11, 12, 13, 13];

    const ROWS: Options = Options {
        span: false,
        mirror: false,
        first_only: true,
        order: Order::Rows,
    };
    const SPAN: Options = Options {
        span: true,
        mirror: true,
        first_only: true,
        order: Order::Rows,
    };
    const BOXED: Options = Options {
        order: Order::BoxedRows,
        ..SPAN
    };
    const ONION: Options = Options {
        order: Order::Onion,
        ..SPAN
    };

    fn exists(n: usize, k: usize, opts: Options) -> bool {
        let out = search(
            n,
            k,
            Options {
                first_only: true,
                ..opts
            },
            4,
            |_, _, _, _, _| {},
        );
        if let Some(ex) = &out.example {
            check(n, ex).unwrap();
            assert_eq!(ex.len(), k);
        }
        out.example.is_some()
    }

    fn count(n: usize, k: usize, opts: Options) -> u64 {
        search(
            n,
            k,
            Options {
                first_only: false,
                ..opts
            },
            2,
            |_, _, _, _, _| {},
        )
        .solutions
    }

    #[test]
    fn reproduces_a271490_up_to_n_10() {
        for n in 2..=10 {
            let a = A271490[n - 1];
            assert!(exists(n, a, ROWS), "n={n}: {a} points must fit");
            assert!(!exists(n, a + 1, ROWS), "n={n}: {} points must not fit", a + 1);
        }
    }

    #[test]
    fn reductions_are_exact() {
        // A k-set fits in n x n iff it fits in (n-1) x (n-1) or a reduced
        // k-set exists in n x n.
        for n in 3..=9 {
            for k in 2..=A271490[n - 1] + 1 {
                let full = exists(n, k, ROWS);
                let smaller = exists(n - 1, k, ROWS);
                for opts in [SPAN, BOXED, ONION] {
                    assert_eq!(full, smaller || exists(n, k, opts), "{:?} n={n} k={k}", opts.order);
                }
            }
        }
    }

    /// All k-sets of the n x n grid with distinct distances whose bounding
    /// box has width n, translated to touch row 0.
    fn spanning_sets(n: usize, k: usize) -> Vec<Vec<(usize, usize)>> {
        fn rec(n: usize, k: usize, i: usize, cur: &mut Vec<(usize, usize)>, out: &mut Vec<Vec<(usize, usize)>>) {
            if cur.len() == k {
                let xs = cur.iter().map(|p| p.0);
                if xs.clone().min() == Some(0) && xs.max() == Some(n - 1) && cur.iter().any(|p| p.1 == 0) {
                    out.push(cur.clone());
                }
                return;
            }
            for j in i..n * n {
                cur.push((j % n, j / n));
                if check(n, cur).is_ok() {
                    rec(n, k, j + 1, cur, out);
                }
                cur.pop();
            }
        }
        let mut out = Vec::new();
        rec(n, k, 0, &mut Vec::new(), &mut out);
        out
    }

    #[test]
    fn counts_agree_with_brute_force() {
        for (n, k) in [(4, 4), (5, 4), (6, 5)] {
            let sets = spanning_sets(n, k);
            // Row order without reductions: all sets touching row 0 and
            // column 0 (the spanning ones are a subset).
            let normalised = |s: &Vec<(usize, usize)>| s.iter().any(|p| p.0 == 0) && s.iter().any(|p| p.1 == 0);
            assert_eq!(
                count(
                    n,
                    k,
                    Options {
                        span: true,
                        mirror: false,
                        ..ROWS
                    }
                ),
                sets.iter().filter(|s| normalised(s)).count() as u64
            );
            // Box orders keep at least one set per symmetry orbit of the
            // bounding box and never more than all sets.
            let mut orbits = std::collections::HashSet::new();
            for s in &sets {
                let h = s.iter().map(|p| p.1).max().unwrap() + 1;
                let mut best: Option<Vec<(usize, usize)>> = None;
                for t in 0..8 {
                    if t >= 4 && h != n {
                        continue;
                    }
                    let mut img: Vec<(usize, usize)> = s
                        .iter()
                        .map(|&(x, y)| {
                            let (mut x, mut y) = if t >= 4 { (y, x) } else { (x, y) };
                            if t & 1 == 1 {
                                x = n - 1 - x;
                            }
                            if t & 2 == 2 {
                                y = h - 1 - y;
                            }
                            (x, y)
                        })
                        .collect();
                    img.sort_unstable();
                    if best.as_ref().map_or(true, |b| img < *b) {
                        best = Some(img);
                    }
                }
                orbits.insert(best.unwrap());
            }
            for opts in [BOXED, ONION] {
                let c = count(n, k, opts);
                assert!(
                    c >= orbits.len() as u64 && c <= sets.len() as u64,
                    "{:?} n={n} k={k}: {c} vs {} orbits",
                    opts.order,
                    orbits.len()
                );
            }
        }
    }
}
