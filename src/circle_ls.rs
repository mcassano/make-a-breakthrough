//! Iterated local search for large no-four-on-a-circle sets.
//!
//! The state is always a valid set `S` (at most three points on every
//! circle or line). A curve with three points of `S` is *full*; a cell
//! outside `S` is *free* if no full curve passes through it, and free cells
//! can be added. When none is free the search makes plateau moves: it picks
//! a cell blocked by exactly one full curve, removes one of that curve's
//! three points and adds the cell. Removed points are tabu for a while so
//! the walk does not undo itself. Occasionally it forces in a cell blocked
//! by several curves (removing one point from each), which shrinks `S`
//! temporarily and lets the search escape plateaus.
//!
//! This is the scheme of the iterated local search of Andrade, Resende and
//! Werneck for maximum independent set, applied to the 4-uniform
//! hypergraph of concyclic/collinear quadruples.

use crate::circle::{curves, Point};
use crate::local::Rng;

pub struct Instance {
    pub n: usize,
    /// Cells of each curve (indices `x * n + y`).
    curve_cells: Vec<u32>,
    curve_start: Vec<u32>,
    /// Curves through each cell.
    cell_curves: Vec<u32>,
    cell_start: Vec<u32>,
}

impl Instance {
    pub fn new(n: usize) -> Instance {
        let cs = curves(n);
        let mut curve_cells = Vec::new();
        let mut curve_start = vec![0u32];
        let mut deg = vec![0u32; n * n];
        for c in &cs {
            for &v in c {
                deg[v as usize] += 1;
            }
            curve_cells.extend_from_slice(c);
            curve_start.push(curve_cells.len() as u32);
        }
        let mut cell_start = vec![0u32; n * n + 1];
        for v in 0..n * n {
            cell_start[v + 1] = cell_start[v] + deg[v];
        }
        let mut fill = cell_start.clone();
        let mut cell_curves = vec![0u32; cell_start[n * n] as usize];
        for (ci, c) in cs.iter().enumerate() {
            for &v in c {
                cell_curves[fill[v as usize] as usize] = ci as u32;
                fill[v as usize] += 1;
            }
        }
        Instance {
            n,
            curve_cells,
            curve_start,
            cell_curves,
            cell_start,
        }
    }

    fn curves_of(&self, v: usize) -> &[u32] {
        &self.cell_curves[self.cell_start[v] as usize..self.cell_start[v + 1] as usize]
    }

    fn cells_of(&self, c: usize) -> &[u32] {
        &self.curve_cells[self.curve_start[c] as usize..self.curve_start[c + 1] as usize]
    }

    pub fn num_curves(&self) -> usize {
        self.curve_start.len() - 1
    }
}

pub struct Search<'a> {
    inst: &'a Instance,
    in_set: Vec<bool>,
    /// Points of `S` on each curve.
    count: Vec<u8>,
    /// For cells outside `S`: number of full curves through the cell.
    blocked: Vec<u32>,
    /// Iteration until which a cell may not be added (after removal).
    tabu_until: Vec<u64>,
    set: Vec<u32>,
    free: Vec<u32>,
    free_pos: Vec<u32>,
    pub iter: u64,
}

const NONE: u32 = u32::MAX;

impl<'a> Search<'a> {
    pub fn new(inst: &'a Instance) -> Search<'a> {
        let cells = inst.n * inst.n;
        let mut s = Search {
            inst,
            in_set: vec![false; cells],
            count: vec![0; inst.num_curves()],
            blocked: vec![0; cells],
            tabu_until: vec![0; cells],
            set: Vec::new(),
            free: Vec::new(),
            free_pos: vec![NONE; cells],
            iter: 0,
        };
        for v in 0..cells {
            s.mark_free(v);
        }
        s
    }

    pub fn size(&self) -> usize {
        self.set.len()
    }

    pub fn points(&self) -> Vec<Point> {
        let n = self.inst.n as u32;
        let mut pts: Vec<Point> = self.set.iter().map(|&v| ((v / n) as i64, (v % n) as i64)).collect();
        pts.sort_unstable();
        pts
    }

    fn mark_free(&mut self, v: usize) {
        if self.free_pos[v] == NONE {
            self.free_pos[v] = self.free.len() as u32;
            self.free.push(v as u32);
        }
    }

    fn unmark_free(&mut self, v: usize) {
        let p = self.free_pos[v];
        if p != NONE {
            let last = self.free.pop().unwrap();
            if last as usize != v {
                self.free[p as usize] = last;
                self.free_pos[last as usize] = p;
            }
            self.free_pos[v] = NONE;
        }
    }

    fn add(&mut self, v: usize) {
        debug_assert!(!self.in_set[v] && self.blocked[v] == 0);
        self.in_set[v] = true;
        self.unmark_free(v);
        self.set.push(v as u32);
        let inst = self.inst;
        for &c in inst.curves_of(v) {
            let c = c as usize;
            self.count[c] += 1;
            if self.count[c] == 3 {
                for &w in inst.cells_of(c) {
                    let w = w as usize;
                    if !self.in_set[w] {
                        self.blocked[w] += 1;
                        if self.blocked[w] == 1 {
                            self.unmark_free(w);
                        }
                    }
                }
            }
        }
    }

    fn remove(&mut self, u: usize) {
        debug_assert!(self.in_set[u]);
        self.in_set[u] = false;
        let pos = self.set.iter().position(|&x| x as usize == u).unwrap();
        self.set.swap_remove(pos);
        let inst = self.inst;
        let mut blocked_u = 0;
        for &c in inst.curves_of(u) {
            let c = c as usize;
            if self.count[c] == 3 {
                for &w in inst.cells_of(c) {
                    let w = w as usize;
                    if !self.in_set[w] && w != u {
                        self.blocked[w] -= 1;
                        if self.blocked[w] == 0 {
                            self.mark_free(w);
                        }
                    }
                }
            }
            self.count[c] -= 1;
            if self.count[c] == 3 {
                blocked_u += 1;
            }
        }
        self.blocked[u] = blocked_u;
        if blocked_u == 0 {
            self.mark_free(u);
        }
    }

    /// Adds free, non-tabu cells in random order until none is left.
    fn greedy_fill(&mut self, rng: &mut Rng) {
        loop {
            // Collect non-tabu free cells.
            let mut picked = None;
            let len = self.free.len();
            if len == 0 {
                return;
            }
            let start = rng.below(len);
            for i in 0..len {
                let v = self.free[(start + i) % len] as usize;
                if self.tabu_until[v] <= self.iter {
                    picked = Some(v);
                    break;
                }
            }
            match picked {
                Some(v) => self.add(v),
                None => return,
            }
        }
    }

    /// Force `v` into the set by removing one point of every full curve
    /// through it. Returns the removed points.
    fn force_add(&mut self, v: usize, rng: &mut Rng) -> Vec<usize> {
        let inst = self.inst;
        let mut removed = Vec::new();
        loop {
            // Find a full curve through v (if any) and drop one of its points.
            let mut full = None;
            for &c in inst.curves_of(v) {
                if self.count[c as usize] >= 3 {
                    full = Some(c as usize);
                    break;
                }
            }
            let Some(c) = full else { break };
            let members: Vec<usize> = inst
                .cells_of(c)
                .iter()
                .map(|&w| w as usize)
                .filter(|&w| self.in_set[w])
                .collect();
            let u = members[rng.below(members.len())];
            self.remove(u);
            removed.push(u);
        }
        self.add(v);
        removed
    }

    /// Runs the iterated local search for `iters` iterations; calls
    /// `report` whenever the set reaches a new maximum size.
    pub fn run(&mut self, rng: &mut Rng, iters: u64, tabu: u64, mut report: impl FnMut(&Self)) -> Vec<Point> {
        let cells = self.inst.n * self.inst.n;
        self.greedy_fill(rng);
        let mut best = self.size();
        let mut best_pts = self.points();
        report(self);
        let mut since_best = 0u64;
        while self.iter < iters {
            self.iter += 1;
            since_best += 1;
            // Candidate: a random non-member cell, preferring those blocked
            // by a single full curve (plateau swaps).
            let mut v = NONE as usize;
            for _ in 0..64 {
                let w = rng.below(cells);
                if !self.in_set[w] && self.tabu_until[w] <= self.iter && self.blocked[w] <= 1 {
                    v = w;
                    break;
                }
            }
            if v == NONE as usize {
                loop {
                    let w = rng.below(cells);
                    if !self.in_set[w] {
                        v = w;
                        break;
                    }
                }
            }
            let removed = self.force_add(v, rng);
            for &u in &removed {
                self.tabu_until[u] = self.iter + tabu + rng.below(tabu as usize + 1) as u64;
            }
            self.greedy_fill(rng);
            if self.size() > best {
                best = self.size();
                best_pts = self.points();
                since_best = 0;
                report(self);
            } else if since_best > 50_000 && self.size() + 3 < best {
                // Drifted too far below the best: shake more gently by
                // continuing from the current state (keep the walk simple).
                since_best = 0;
            }
        }
        best_pts
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circle::check;

    #[test]
    fn search_keeps_sets_valid_and_finds_known_optima() {
        // Known optimum values from OEIS A395821.
        for (n, opt) in [(5, 9), (6, 11), (7, 14)] {
            let inst = Instance::new(n);
            let mut rng = Rng::new(3);
            let mut s = Search::new(&inst);
            let best = s.run(&mut rng, 20_000, 5, |st| {
                check(n, &st.points()).unwrap();
            });
            check(n, &best).unwrap();
            assert_eq!(best.len(), opt, "n={n}");
        }
    }
}
