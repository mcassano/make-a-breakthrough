//! Local search over the "slot" representation of symmetric configurations.
//!
//! Put the origin at the grid center. Under quarter-turn symmetry every
//! row `y` and its mirror row `-y` hold the same number of points, and a
//! quarter-turn orbit `{(x,y), (-y,x), (-x,-y), (y,-x)}` with `x, y > 0`
//! places exactly one point in each of the row pairs `|y|` and `|x|`.
//! Numbering the `m = n/2` row pairs `1..=m`, a rot4 configuration with two
//! points in every row is therefore the same thing as a sequence of the
//! multiset `{1, 1, 2, 2, ..., m, m}` read as `m` ordered pairs `(x, y)`.
//!
//! Exchanging two entries of the sequence keeps every row (and, by symmetry,
//! every column) at exactly two points, so the search only has to eliminate
//! collinear triples. For odd `n` and the `rct4` class the middle row is
//! filled by one axis orbit `{(s,0), (0,s), (-s,0), (0,-s)}` and the
//! remaining symmetry defect is one diagonal pair `{(d,d), (-d,-d)}`; both
//! consume one entry each, and the other `m - 1` ordered pairs must be
//! off-diagonal. See the README for the derivation.

use crate::geometry::{gcd, Cell};
use crate::symmetry::Class;

/// Small, fast, deterministic PRNG (xoshiro256**).
pub struct Rng([u64; 4]);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ (x >> 31)
        };
        Rng([next(), next(), next(), next()])
    }

    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.0;
        let result = s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// Uniform integer in `0..k`.
    pub fn below(&mut self, k: usize) -> usize {
        (((self.next_u64() >> 32) * k as u64) >> 32) as usize
    }

    /// Uniform float in `[0, 1)`.
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

/// Dense index of every line of the grid, addressed through the difference
/// vector of two of its cells.
pub struct LineIndex {
    n: usize,
    /// For a difference `(dr, dc)` (offset by `n - 1`): the primitive
    /// direction `(a, b)` and the base index of that direction's lines.
    dir: Vec<(i32, i32, u32)>,
    pub num_lines: usize,
}

impl LineIndex {
    pub fn new(n: usize) -> LineIndex {
        let w = 2 * n - 1;
        let n1 = n as i64 - 1;
        let mut dir = vec![(0, 0, 0); w * w];
        let mut base_of = std::collections::HashMap::new();
        let mut next_base: u64 = 0;
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
                // Lines of direction (a, b) are labelled by k = b*r - a*c,
                // which ranges over an interval of length (a + |b|)(n - 1).
                let base = *base_of.entry((a, b)).or_insert_with(|| {
                    let base = next_base;
                    next_base += ((a + b.abs()) * n1 + 1) as u64;
                    base
                });
                let idx = (dr + n1) as usize * w + (dc + n1) as usize;
                dir[idx] = (a as i32, b as i32, u32::try_from(base).expect("too many lines"));
            }
        }
        LineIndex {
            n,
            dir,
            num_lines: next_base as usize,
        }
    }

    /// Index of the line through two distinct cells.
    #[inline]
    pub fn line(&self, p: Cell, q: Cell) -> usize {
        let n1 = self.n as i32 - 1;
        let (pr, pc) = (p.0 as i32, p.1 as i32);
        let (dr, dc) = (q.0 as i32 - pr, q.1 as i32 - pc);
        let (a, b, base) = self.dir[((dr + n1) * (2 * n1 + 1) + dc + n1) as usize];
        // Shift k so that it starts at 0 for this direction.
        let k = b * pr - a * pc + a * n1 + if b < 0 { -b * n1 } else { 0 };
        base as usize + k as usize
    }
}

/// How the slot sequence maps to points.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Layout {
    /// Even `n`, quarter turn: `m` ordered pairs.
    Rot4,
    /// Odd `n`, nearly quarter turn: `m - 1` ordered pairs, an axis slot and
    /// a diagonal slot.
    Rct4,
}

/// The configuration state and incremental collinearity bookkeeping.
pub struct Search {
    pub n: usize,
    m: usize,
    layout: Layout,
    /// The slot values (row-pair numbers `1..=m`).
    slots: Vec<u16>,
    /// Number of slot groups (each yields a point orbit).
    groups: usize,
    index: LineIndex,
    /// Number of pairs of current points on each line.
    pairs: Vec<u16>,
    /// `points_of_pairs[p] = k` where `p = k(k-1)/2`, else 0.
    count_of_pairs: Vec<u16>,
    /// Current points as cells, grouped by slot group (fixed size per group).
    points: Vec<Vec<Cell>>,
    /// Twice the number of collinear triples.
    pub cost: i64,
}

impl Search {
    pub fn new(n: usize, class: Class, rng: &mut Rng) -> Result<Search, String> {
        let layout = match class {
            Class::Rot4 if n % 2 == 0 => Layout::Rot4,
            Class::Rct4 if n % 2 == 1 => Layout::Rct4,
            _ => {
                return Err(format!(
                    "local search supports rot4 (even n) and rct4 (odd n), not {} with n={n}",
                    class.name()
                ))
            }
        };
        let m = n / 2;
        let slots: Vec<u16> = (1..=m as u16).flat_map(|v| [v, v]).collect();
        let groups = match layout {
            Layout::Rot4 => m,
            Layout::Rct4 => m + 1,
        };
        let index = LineIndex::new(n);
        let max_pairs = n * (n - 1) / 2;
        let mut count_of_pairs = vec![0u16; max_pairs + 1];
        for k in 2..=n {
            count_of_pairs[k * (k - 1) / 2] = k as u16;
        }
        let mut s = Search {
            n,
            m,
            layout,
            slots,
            groups,
            pairs: vec![0; index.num_lines],
            index,
            count_of_pairs,
            points: vec![Vec::new(); groups],
            cost: 0,
        };
        // Random legal starting state (Fisher-Yates shuffles until legal).
        loop {
            for i in (1..s.slots.len()).rev() {
                s.slots.swap(i, rng.below(i + 1));
            }
            if s.is_legal() {
                break;
            }
        }
        for g in 0..groups {
            s.place_group(g);
        }
        Ok(s)
    }

    /// Grid cell of the centered coordinates `(x, y)` (x: column, y: row).
    /// For even `n` a value `v` stands for the half-integer offset `v - 1/2`.
    fn cell(&self, x: i32, y: i32) -> Cell {
        let conv = |u: i32| -> usize {
            if self.n % 2 == 0 {
                (if u > 0 {
                    self.m as i32 - 1 + u
                } else {
                    self.m as i32 + u
                }) as usize
            } else {
                (self.m as i32 + u) as usize
            }
        };
        (conv(y), conv(x))
    }

    fn rot4_orbit(&self, x: i32, y: i32) -> Vec<Cell> {
        vec![self.cell(x, y), self.cell(-y, x), self.cell(-x, -y), self.cell(y, -x)]
    }

    /// The points of slot group `g` under the current slot values, and
    /// whether the group is structurally invalid.
    fn group_points(&self, g: usize) -> (Vec<Cell>, bool) {
        let v = |i: usize| self.slots[i] as i32;
        match self.layout {
            Layout::Rot4 => (self.rot4_orbit(v(2 * g), v(2 * g + 1)), false),
            Layout::Rct4 => {
                let pairs = self.m - 1;
                if g < pairs {
                    let (x, y) = (v(2 * g), v(2 * g + 1));
                    (self.rot4_orbit(x, y), x == y)
                } else if g == pairs {
                    (self.rot4_orbit(v(2 * pairs), 0), false)
                } else {
                    let d = v(2 * pairs + 1);
                    (vec![self.cell(d, d), self.cell(-d, -d)], false)
                }
            }
        }
    }

    /// The slot group containing slot `i`.
    fn group_of(&self, i: usize) -> usize {
        match self.layout {
            Layout::Rot4 => i / 2,
            Layout::Rct4 => (i / 2).min(self.m - 1) + usize::from(i == 2 * self.m - 1),
        }
    }

    /// Adds a point, updating the cost by the number of new collinear triples.
    ///
    /// If the line through `p` and `q` already holds `k` points, adding `p`
    /// creates `C(k, 2)` triples; that line is met once for each of its `k`
    /// points, so summing `k - 1` over all `q` gives twice the new triples.
    fn add_point(&mut self, p: Cell) {
        let mut delta = 0i64;
        for &q in self.points.iter().flatten() {
            let l = self.index.line(p, q);
            delta += self.count_of_pairs[self.pairs[l] as usize].max(1) as i64 - 1;
        }
        for &q in self.points.iter().flatten() {
            self.pairs[self.index.line(p, q)] += 1;
        }
        self.cost += delta;
    }

    /// Removes a point that is not in `self.points`; the mirror of `add_point`.
    fn remove_point(&mut self, p: Cell) {
        for &q in self.points.iter().flatten() {
            self.pairs[self.index.line(p, q)] -= 1;
        }
        let mut delta = 0i64;
        for &q in self.points.iter().flatten() {
            let l = self.index.line(p, q);
            delta += self.count_of_pairs[self.pairs[l] as usize].max(1) as i64 - 1;
        }
        self.cost -= delta;
    }

    fn place_group(&mut self, g: usize) {
        let (pts, _) = self.group_points(g);
        for p in pts {
            self.add_point(p);
            self.points[g].push(p);
        }
    }

    fn unplace_group(&mut self, g: usize) {
        while let Some(p) = self.points[g].pop() {
            self.remove_point(p);
        }
    }

    /// The orbit label of group `g`: its slot values, used to detect groups
    /// that would produce the same points.
    fn group_key(&self, g: usize) -> (u16, u16, u8) {
        let v = |i: usize| self.slots[i];
        match self.layout {
            Layout::Rot4 => (v(2 * g), v(2 * g + 1), 0),
            Layout::Rct4 if g < self.m - 1 => (v(2 * g), v(2 * g + 1), 0),
            Layout::Rct4 if g == self.m - 1 => (v(2 * g), 0, 1),
            Layout::Rct4 => (v(2 * g - 1), 0, 2),
        }
    }

    /// Whether group `g` is a legal orbit that no other group duplicates.
    fn group_ok(&self, g: usize) -> bool {
        let key = self.group_key(g);
        if self.group_points(g).1 {
            return false;
        }
        (0..self.groups).all(|h| h == g || self.group_key(h) != key)
    }

    /// Whether every group is legal and distinct.
    pub fn is_legal(&self) -> bool {
        (0..self.groups).all(|g| self.group_ok(g))
    }

    /// Whether exchanging slots `i` and `j` keeps the state legal.
    pub fn swap_is_legal(&mut self, i: usize, j: usize) -> bool {
        self.slots.swap(i, j);
        let ok = self.group_ok(self.group_of(i)) && self.group_ok(self.group_of(j));
        self.slots.swap(i, j);
        ok
    }

    /// Exchanges the values of slots `i` and `j`, updating the cost.
    pub fn swap_slots(&mut self, i: usize, j: usize) {
        let (gi, gj) = (self.group_of(i), self.group_of(j));
        self.unplace_group(gi);
        if gj != gi {
            self.unplace_group(gj);
        }
        self.slots.swap(i, j);
        self.place_group(gi);
        if gj != gi {
            self.place_group(gj);
        }
    }

    pub fn num_slots(&self) -> usize {
        self.slots.len()
    }

    pub fn slot(&self, i: usize) -> u16 {
        self.slots[i]
    }

    /// All current points.
    pub fn points(&self) -> Vec<Cell> {
        let mut pts: Vec<Cell> = self.points.iter().flatten().copied().collect();
        pts.sort_unstable();
        pts
    }

    /// Recomputes the cost from scratch (for testing the incremental updates).
    pub fn recompute_cost(&self) -> i64 {
        let pts = self.points();
        let mut per_line = std::collections::HashMap::new();
        for (i, &p) in pts.iter().enumerate() {
            for &q in &pts[i + 1..] {
                *per_line.entry(self.index.line(p, q)).or_insert(0usize) += 1;
            }
        }
        // A line with k points has k(k-1)/2 pairs and contributes
        // 2 * C(k, 3) = k(k-1)(k-2)/3 to the cost.
        per_line
            .values()
            .map(|&pairs| {
                let k = self.count_of_pairs[pairs] as i64;
                k * (k - 1) * (k - 2) / 3
            })
            .sum()
    }
}

/// Simulated annealing over slot exchanges. Returns the points of a solution
/// if one is found within `max_steps` steps.
pub fn anneal(search: &mut Search, rng: &mut Rng, max_steps: u64, t0: f64, t1: f64) -> Option<Vec<Cell>> {
    let k = search.num_slots();
    let mut t = t0;
    let cool = (t1 / t0).powf(1.0 / max_steps as f64);
    for _ in 0..max_steps {
        if search.cost == 0 {
            return Some(search.points());
        }
        let i = rng.below(k);
        let j = rng.below(k);
        if search.slot(i) == search.slot(j) || !search.swap_is_legal(i, j) {
            continue;
        }
        let before = search.cost;
        search.swap_slots(i, j);
        let delta = search.cost - before;
        if delta > 0 && rng.unit() >= (-(delta as f64) / t).exp() {
            search.swap_slots(i, j);
        }
        t *= cool;
    }
    (search.cost == 0).then(|| search.points())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{check_solution, line_key};

    #[test]
    fn line_index_agrees_with_line_key() {
        let n = 9;
        let idx = LineIndex::new(n);
        let cells: Vec<Cell> = (0..n).flat_map(|r| (0..n).map(move |c| (r, c))).collect();
        let mut by_key = std::collections::HashMap::new();
        for &p in &cells {
            for &q in &cells {
                if p != q {
                    let l = idx.line(p, q);
                    assert!(l < idx.num_lines);
                    assert_eq!(*by_key.entry(line_key(p, q)).or_insert(l), l);
                }
            }
        }
        let distinct: std::collections::HashSet<usize> = by_key.values().copied().collect();
        assert_eq!(distinct.len(), by_key.len());
    }

    #[test]
    fn incremental_cost_matches_recomputation() {
        for (n, class) in [(12, Class::Rot4), (13, Class::Rct4)] {
            let mut rng = Rng::new(7);
            let mut s = Search::new(n, class, &mut rng).unwrap();
            assert_eq!(s.cost, s.recompute_cost());
            for _ in 0..500 {
                let (i, j) = (rng.below(s.num_slots()), rng.below(s.num_slots()));
                if !s.swap_is_legal(i, j) {
                    continue;
                }
                s.swap_slots(i, j);
                assert!(s.is_legal());
                assert_eq!(s.cost, s.recompute_cost());
            }
        }
    }

    #[test]
    fn every_state_has_two_points_per_row_and_column() {
        for (n, class) in [(10, Class::Rot4), (11, Class::Rct4)] {
            let mut rng = Rng::new(3);
            let s = Search::new(n, class, &mut rng).unwrap();
            let pts = s.points();
            assert_eq!(pts.len(), 2 * n);
            for i in 0..n {
                assert_eq!(pts.iter().filter(|p| p.0 == i).count(), 2);
                assert_eq!(pts.iter().filter(|p| p.1 == i).count(), 2);
            }
        }
    }

    #[test]
    fn annealing_solves_small_grids() {
        // rct4 has no solutions for n = 11, 13, 15.
        for (n, class) in [
            (10, Class::Rot4),
            (16, Class::Rot4),
            (9, Class::Rct4),
            (17, Class::Rct4),
        ] {
            let mut rng = Rng::new(1);
            let found = (0..50).find_map(|_| {
                let mut s = Search::new(n, class, &mut rng).unwrap();
                anneal(&mut s, &mut rng, 200_000, 3.0, 0.2)
            });
            let pts = found.unwrap_or_else(|| panic!("no solution for n={n}"));
            check_solution(n, &pts).unwrap();
        }
    }
}
