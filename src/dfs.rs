//! Randomized depth-first search for quarter-turn symmetric (`rot4`)
//! configurations on even grids.
//!
//! As explained in [`crate::local`], with the origin at the grid center a
//! rot4 configuration is a set of `m = n/2` quadrant cells `(x, y)`,
//! `1 <= x, y <= m`, each standing for its orbit under the quarter turn, such
//! that every value `v` is used exactly twice among all coordinates (the
//! "degree" of the row pair `v`). Viewing `(x, y)` as an edge between `x` and
//! `y`, the configuration is a 2-regular multigraph with oriented edges.
//!
//! The search keeps, for every quadrant cell, whether its orbit is *blocked*:
//! it lies on a line through two chosen points, or choosing it would put
//! three points on a line with an already chosen orbit. Because the chosen
//! set is invariant under the quarter turn, the blocked set is too, so when
//! an orbit with representative `P` is added it suffices to walk the lines
//! through `P` and each other chosen point: the lines through the other three
//! points of the orbit are their rotations and block the same orbits.
//!
//! Blocked flags are stored twice, as one 64-bit mask per quadrant row and
//! per quadrant column, so the number of usable edges at a vertex is two
//! popcounts. The search branches on the vertex with the fewest usable edges
//! and restarts with a fresh random ordering after a node budget.

use crate::geometry::{gcd, Cell};
use crate::local::Rng;

/// Line lookup restricted to lines containing at least three grid cells.
struct LongLines {
    n: usize,
    /// For each difference vector: `(a, b, base)` of its primitive
    /// direction, or `base = u32::MAX` if lines in that direction have fewer
    /// than three cells.
    dir: Vec<(i32, i32, u32)>,
    /// Offsets into `orbits` for each line id (`len = num_lines + 1`).
    start: Vec<u32>,
    /// Orbit ids of the cells of every line, concatenated.
    orbits: Vec<u16>,
}

impl LongLines {
    fn new(n: usize, orbit_of: impl Fn(Cell) -> u16) -> LongLines {
        let n1 = n as i64 - 1;
        let w = 2 * n - 1;
        let max_step = n1 / 2;
        let mut dir = vec![(0, 0, u32::MAX); w * w];
        let mut bases = std::collections::HashMap::new();
        let mut num_lines: u64 = 0;
        for dr in -n1..=n1 {
            for dc in -n1..=n1 {
                if dr == 0 && dc == 0 {
                    continue;
                }
                let g = gcd(dr, dc);
                let (mut a, mut b) = (dr / g, dc / g);
                if a < 0 || (a == 0 && b < 0) {
                    a = -a;
                    b = -b;
                }
                if a.max(b.abs()) > max_step {
                    continue;
                }
                let base = *bases.entry((a, b)).or_insert_with(|| {
                    let base = num_lines;
                    num_lines += ((a + b.abs()) * n1 + 1) as u64;
                    base
                });
                dir[(dr + n1) as usize * w + (dc + n1) as usize] = (a as i32, b as i32, base as u32);
            }
        }
        let mut lists: Vec<Vec<u16>> = vec![Vec::new(); num_lines as usize];
        let mut ll = LongLines {
            n,
            dir,
            start: Vec::new(),
            orbits: Vec::new(),
        };
        for line in crate::geometry::lines(n, 3) {
            let id = ll.line(line[0], line[1]).expect("long line");
            lists[id] = line.iter().map(|&p| orbit_of(p)).collect();
        }
        ll.start.reserve(lists.len() + 1);
        for list in &lists {
            ll.start.push(ll.orbits.len() as u32);
            ll.orbits.extend_from_slice(list);
        }
        ll.start.push(ll.orbits.len() as u32);
        ll
    }

    /// Id of the line through two distinct cells if it has >= 3 cells.
    #[inline]
    fn line(&self, p: Cell, q: Cell) -> Option<usize> {
        let n1 = self.n as i32 - 1;
        let (pr, pc) = (p.0 as i32, p.1 as i32);
        let (dr, dc) = (q.0 as i32 - pr, q.1 as i32 - pc);
        let (a, b, base) = self.dir[((dr + n1) * (2 * n1 + 1) + dc + n1) as usize];
        if base == u32::MAX {
            return None;
        }
        let k = b * pr - a * pc + a * n1 + if b < 0 { -b * n1 } else { 0 };
        Some(base as usize + k as usize)
    }

    #[inline]
    fn orbits_on(&self, line: usize) -> &[u16] {
        &self.orbits[self.start[line] as usize..self.start[line + 1] as usize]
    }
}

/// Search state at one node of the tree; copied on descent.
#[derive(Clone)]
pub(crate) struct State {
    /// `row[y-1]` bit `x-1`: orbit `(x, y)` is blocked.
    pub(crate) row: Vec<u64>,
    /// `col[x-1]` bit `y-1`: orbit `(x, y)` is blocked.
    pub(crate) col: Vec<u64>,
    /// Remaining degree of each vertex.
    pub(crate) need: Vec<u8>,
    /// Bit `v-1` set iff `need[v-1] > 0` and the vertex is still open.
    pub(crate) active: u64,
    /// Total missing degree of vertices given up on (large neighbourhood
    /// search only; always 0 in the exact search).
    pub(crate) lost: u32,
}

pub struct Rot4Dfs {
    pub n: usize,
    pub(crate) m: usize,
    lines: LongLines,
    /// The four cells of each orbit; `cells[o][0]` is the quadrant cell.
    pub(crate) cells: Vec<[Cell; 4]>,
    /// Orbits that cannot be chosen together with each orbit because one of
    /// its lines through two of its own points contains them.
    conflicts: Vec<Vec<u16>>,
    pub(crate) chosen: Vec<u16>,
    pub nodes: u64,
}

impl Rot4Dfs {
    pub fn new(n: usize) -> Result<Rot4Dfs, String> {
        if n % 2 != 0 || n < 4 {
            return Err("rot4 search needs an even n >= 4".into());
        }
        let m = n / 2;
        if m > 64 {
            return Err("n > 128 is not supported".into());
        }
        // Quadrant coordinate v in 1..=m is the half-integer offset v - 1/2
        // from the center; grid index m - 1 + v (positive) or m - v (negative).
        let coord = |u: i32| -> usize {
            if u > 0 {
                (m as i32 - 1 + u) as usize
            } else {
                (m as i32 + u) as usize
            }
        };
        let cell = |x: i32, y: i32| -> Cell { (coord(y), coord(x)) };
        let mut cells = Vec::with_capacity(m * m);
        for y in 1..=m as i32 {
            for x in 1..=m as i32 {
                cells.push([cell(x, y), cell(-y, x), cell(-x, -y), cell(y, -x)]);
            }
        }
        let orbit_id = |(r, c): Cell| -> u16 {
            // Rotate the cell into the quadrant x, y > 0.
            let (mut x, mut y) = (
                if c >= m {
                    (c + 1 - m) as i32
                } else {
                    c as i32 - m as i32
                },
                if r >= m {
                    (r + 1 - m) as i32
                } else {
                    r as i32 - m as i32
                },
            );
            while !(x > 0 && y > 0) {
                (x, y) = (y, -x);
            }
            ((y - 1) as usize * m + (x - 1) as usize) as u16
        };
        let lines = LongLines::new(n, orbit_id);
        let mut conflicts = vec![Vec::new(); m * m];
        for (o, pts) in cells.iter().enumerate() {
            for j in 1..3 {
                if let Some(l) = lines.line(pts[0], pts[j]) {
                    for &q in lines.orbits_on(l) {
                        if q as usize != o {
                            conflicts[o].push(q);
                            conflicts[q as usize].push(o as u16);
                        }
                    }
                }
            }
        }
        for c in &mut conflicts {
            c.sort_unstable();
            c.dedup();
        }
        Ok(Rot4Dfs {
            n,
            m,
            lines,
            cells,
            conflicts,
            chosen: Vec::new(),
            nodes: 0,
        })
    }

    pub(crate) fn block(&self, s: &mut State, o: u16) {
        let (x, y) = (o as usize % self.m, o as usize / self.m);
        s.row[y] |= 1 << x;
        s.col[x] |= 1 << y;
    }

    /// Adds orbit `o` to the chosen set and blocks everything it rules out.
    pub(crate) fn choose(&mut self, s: &mut State, o: u16) {
        let (x, y) = (o as usize % self.m, o as usize / self.m);
        s.need[x] -= 1;
        s.need[y] -= 1;
        for v in [x, y] {
            if s.need[v] == 0 {
                s.active &= !(1 << v);
            }
        }
        let p = self.cells[o as usize][0];
        for &q in &self.conflicts[o as usize] {
            self.block(s, q);
        }
        for &c in &self.chosen {
            for &q in &self.cells[c as usize] {
                if let Some(l) = self.lines.line(p, q) {
                    for &b in self.lines.orbits_on(l) {
                        self.block(s, b);
                    }
                }
            }
        }
        self.block(s, o);
        self.chosen.push(o);
    }

    /// Usable edges at vertex `v` (loops count twice, as they fill both
    /// slots), restricted to still-active partner vertices.
    pub(crate) fn usable(&self, s: &State, v: usize) -> u32 {
        let mut free_row = !s.row[v] & s.active;
        let free_col = !s.col[v] & s.active & !(1 << v);
        if s.need[v] < 2 {
            free_row &= !(1 << v);
        }
        free_row.count_ones() + free_col.count_ones() + (free_row >> v & 1) as u32
    }

    /// Candidate orbits at vertex `v`.
    pub(crate) fn candidates(&self, s: &State, v: usize, out: &mut Vec<u16>) {
        out.clear();
        let m = self.m;
        let mut free_row = !s.row[v] & s.active;
        if s.need[v] < 2 {
            free_row &= !(1 << v);
        }
        let mut free_col = !s.col[v] & s.active & !(1 << v);
        while free_row != 0 {
            let x = free_row.trailing_zeros() as usize;
            free_row &= free_row - 1;
            out.push((v * m + x) as u16);
        }
        while free_col != 0 {
            let y = free_col.trailing_zeros() as usize;
            free_col &= free_col - 1;
            out.push((y * m + v) as u16);
        }
    }

    /// The open vertex with the fewest usable edges per missing edge (ties
    /// broken at random), or `None` if some vertex can no longer be
    /// completed. Of the orderings we measured this gives the smallest trees;
    /// branching on the outermost rows first is about 40x worse.
    fn pick_vertex(&self, s: &State, rng: Option<&mut Rng>) -> Option<usize> {
        let mut best = (u32::MAX, 0usize);
        let mut ties = 0;
        let mut rng = rng;
        let mut act = s.active;
        while act != 0 {
            let v = act.trailing_zeros() as usize;
            act &= act - 1;
            let u = self.usable(s, v);
            let need = s.need[v] as u32;
            if u < need {
                return None;
            }
            let score = u * 2 / need;
            if score < best.0 {
                best = (score, v);
                ties = 1;
            } else if score == best.0 {
                ties += 1;
                if let Some(r) = rng.as_deref_mut() {
                    if r.below(ties) == 0 {
                        best = (score, v);
                    }
                }
            }
        }
        Some(best.1)
    }

    /// Depth-first search below `s`. Returns `Some(true)` when a solution
    /// is completed, `Some(false)` when the subtree is exhausted, `None`
    /// when the node budget ran out.
    fn search(&mut self, s: &State, rng: &mut Rng, budget: u64) -> Option<bool> {
        self.nodes += 1;
        if self.nodes > budget {
            return None;
        }
        if s.active == 0 {
            return Some(true);
        }
        let Some(v) = self.pick_vertex(s, Some(rng)) else {
            return Some(false);
        };
        let mut cands = Vec::new();
        self.candidates(s, v, &mut cands);
        for i in (1..cands.len()).rev() {
            cands.swap(i, rng.below(i + 1));
        }
        let mut base = s.clone();
        for &o in &cands {
            let mut child = base.clone();
            self.choose(&mut child, o);
            match self.search(&child, rng, budget) {
                Some(false) => {}
                other => return other,
            }
            self.chosen.pop();
            // Later siblings exclude this orbit, so each set is visited once.
            self.block(&mut base, o);
            if self.usable(&base, v) < base.need[v] as u32 {
                break;
            }
        }
        Some(false)
    }

    pub(crate) fn root(&self) -> State {
        let m = self.m;
        State {
            row: vec![0; m],
            col: vec![0; m],
            need: vec![2; m],
            active: if m == 64 { !0 } else { (1u64 << m) - 1 },
            lost: 0,
        }
    }

    /// One randomized attempt with a node budget. Returns the solution's
    /// points, or `None` if the budget ran out or the (sub)tree is empty.
    pub fn attempt(&mut self, rng: &mut Rng, budget: u64) -> Option<Vec<Cell>> {
        self.chosen.clear();
        self.nodes = 0;
        let root = self.root();
        if self.search(&root, rng, budget) == Some(true) {
            let mut pts: Vec<Cell> = self.chosen.iter().flat_map(|&o| self.cells[o as usize]).collect();
            pts.sort_unstable();
            return Some(pts);
        }
        None
    }

    /// Counts all rot4 solutions (as labelled configurations: a solution
    /// and its mirror image count twice unless they coincide).
    pub fn count_all(&mut self) -> u64 {
        fn rec(d: &mut Rot4Dfs, s: &State, count: &mut u64) {
            d.nodes += 1;
            if s.active == 0 {
                *count += 1;
                return;
            }
            let Some(v) = d.pick_vertex(s, None) else { return };
            let mut cands = Vec::new();
            d.candidates(s, v, &mut cands);
            let mut base = s.clone();
            for &o in &cands {
                let mut child = base.clone();
                d.choose(&mut child, o);
                rec(d, &child, count);
                d.chosen.pop();
                d.block(&mut base, o);
            }
        }
        self.chosen.clear();
        self.nodes = 0;
        let root = self.root();
        let mut count = 0;
        rec(self, &root, &mut count);
        count
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{check_no_three_in_line, check_solution};
    use crate::symmetry::{Class, Model};

    /// Counts rot4 solutions by brute force over the SAT model's variables.
    fn brute_force_count(n: usize) -> u64 {
        fn rec(model: &Model, i: usize, m: usize, chosen: &mut Vec<usize>, count: &mut u64) {
            if chosen.len() == m {
                if check_solution(model.n, &model.points(chosen)).is_ok() {
                    *count += 1;
                }
                return;
            }
            if model.vars.len() - i < m - chosen.len() {
                return;
            }
            chosen.push(i);
            if check_no_three_in_line(model.n, &model.points(chosen)).is_ok() {
                rec(model, i + 1, m, chosen, count);
            }
            chosen.pop();
            rec(model, i + 1, m, chosen, count);
        }
        let model = Model::new(n, Class::Rot4).unwrap();
        let mut count = 0;
        rec(&model, 0, n / 2, &mut Vec::new(), &mut count);
        count
    }

    #[test]
    fn exhaustive_count_matches_brute_force() {
        for n in [4, 6, 8, 10, 12] {
            let mut d = Rot4Dfs::new(n).unwrap();
            assert_eq!(d.count_all(), brute_force_count(n), "n={n}");
        }
    }

    #[test]
    fn exhaustive_count_matches_flammenkamp_database() {
        // 2 * (#rot4 classes) + (#fully symmetric classes) from the database
        // of known configurations, which is complete for these n.
        for (n, expected) in [(14, 26), (16, 26), (18, 14), (20, 32), (22, 16)] {
            assert_eq!(Rot4Dfs::new(n).unwrap().count_all(), expected, "n={n}");
        }
    }

    #[test]
    fn randomized_search_finds_valid_solutions() {
        for n in [10, 20, 30] {
            let mut d = Rot4Dfs::new(n).unwrap();
            let mut rng = Rng::new(5);
            let pts = (0..1000).find_map(|_| d.attempt(&mut rng, 100_000)).expect("solution");
            check_solution(n, &pts).unwrap();
        }
    }
}
