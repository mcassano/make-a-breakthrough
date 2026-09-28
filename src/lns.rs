//! Large neighbourhood search for rot4 configurations.
//!
//! The search keeps a valid partial configuration (no three points in line,
//! at most two points per row) and tries to complete it. Each step removes
//! every orbit touching a few randomly chosen row pairs and re-solves the
//! freed part with the depth-first machinery of [`crate::dfs`], now
//! minimising the total missing degree ("deficit") instead of requiring an
//! exact completion. A step is accepted if the deficit does not increase;
//! a configuration with deficit 0 is a solution.

use crate::dfs::{Rot4Dfs, State};
use crate::geometry::Cell;
use crate::local::Rng;

pub struct Lns {
    pub dfs: Rot4Dfs,
    /// Current orbits.
    pub orbits: Vec<u16>,
    pub deficit: u32,
    /// Nodes spent in repairs.
    pub nodes: u64,
    best: (u32, Vec<u16>),
    budget: u64,
}

impl Lns {
    pub fn new(n: usize, rng: &mut Rng, budget: u64) -> Result<Lns, String> {
        let dfs = Rot4Dfs::new(n)?;
        let mut lns = Lns {
            dfs,
            orbits: Vec::new(),
            deficit: 0,
            nodes: 0,
            best: (0, Vec::new()),
            budget,
        };
        let root = lns.dfs.root();
        lns.dfs.chosen.clear();
        lns.best = (u32::MAX, Vec::new());
        lns.repair(&root, rng);
        lns.orbits = std::mem::take(&mut lns.best.1);
        lns.deficit = lns.best.0;
        Ok(lns)
    }

    fn deficit_bound(&self, s: &State) -> u32 {
        let mut lb = s.lost;
        let mut act = s.active;
        while act != 0 {
            let v = act.trailing_zeros() as usize;
            act &= act - 1;
            lb += (s.need[v] as u32).saturating_sub(self.dfs.usable(s, v));
        }
        lb
    }

    /// Branch and bound: extend `s` to minimise the final deficit, recording
    /// improvements over `self.best` (which must be pre-set to the bound to
    /// beat). Stops after `self.budget` nodes.
    fn repair(&mut self, s: &State, rng: &mut Rng) {
        self.nodes += 1;
        self.dfs.nodes += 1;
        if self.dfs.nodes > self.budget {
            return;
        }
        let lb = self.deficit_bound(s);
        if lb >= self.best.0 {
            return;
        }
        if s.active == 0 {
            self.best = (s.lost, self.dfs.chosen.clone());
            return;
        }
        // Most constrained open vertex (smallest slack, possibly negative).
        let mut best_v = (i64::MAX, 0usize);
        let mut ties = 0;
        let mut act = s.active;
        while act != 0 {
            let v = act.trailing_zeros() as usize;
            act &= act - 1;
            let slack = self.dfs.usable(s, v) as i64 - s.need[v] as i64;
            if slack < best_v.0 {
                best_v = (slack, v);
                ties = 1;
            } else if slack == best_v.0 {
                ties += 1;
                if rng.below(ties) == 0 {
                    best_v = (slack, v);
                }
            }
        }
        let v = best_v.1;
        let mut cands = Vec::new();
        self.dfs.candidates(s, v, &mut cands);
        for i in (1..cands.len()).rev() {
            cands.swap(i, rng.below(i + 1));
        }
        let mut base = s.clone();
        for &o in &cands {
            let mut child = base.clone();
            self.dfs.choose(&mut child, o);
            self.repair(&child, rng);
            self.dfs.chosen.pop();
            if self.best.0 <= s.lost {
                return;
            }
            self.dfs.block(&mut base, o);
        }
        // Give up on v: its missing degree becomes part of the deficit.
        let mut child = base;
        child.lost += child.need[v] as u32;
        child.active &= !(1 << v);
        self.repair(&child, rng);
    }

    /// One destroy-and-repair step over the deficient row pairs plus random
    /// ones, `k` in total. Returns true if the step was accepted.
    pub fn step(&mut self, rng: &mut Rng, k: usize) -> bool {
        let m = self.dfs.m;
        let mut deg = vec![0u8; m];
        for &o in &self.orbits {
            deg[o as usize % m] += 1;
            deg[o as usize / m] += 1;
        }
        let mut freed = 0u64;
        for v in 0..m {
            if deg[v] < 2 {
                freed |= 1 << v;
            }
        }
        while (freed.count_ones() as usize) < k.min(m) {
            freed |= 1 << rng.below(m);
        }
        let kept: Vec<u16> = self
            .orbits
            .iter()
            .copied()
            .filter(|&o| freed >> (o as usize % m) & 1 == 0 && freed >> (o as usize / m) & 1 == 0)
            .collect();
        let mut s = self.dfs.root();
        self.dfs.chosen.clear();
        for &o in &kept {
            self.dfs.choose(&mut s, o);
        }
        self.best = (self.deficit + 1, Vec::new());
        self.dfs.nodes = 0;
        self.repair(&s, rng);
        if self.best.0 <= self.deficit {
            self.orbits = std::mem::take(&mut self.best.1);
            self.deficit = self.best.0;
            true
        } else {
            false
        }
    }

    pub fn points(&self) -> Vec<Cell> {
        let mut pts: Vec<Cell> = self.orbits.iter().flat_map(|&o| self.dfs.cells[o as usize]).collect();
        pts.sort_unstable();
        pts
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{check_no_three_in_line, check_solution};

    #[test]
    fn lns_keeps_valid_states_and_solves_small_grids() {
        for n in [16, 24] {
            let mut rng = Rng::new(11);
            let mut lns = Lns::new(n, &mut rng, 10_000).unwrap();
            for _ in 0..20_000 {
                if lns.deficit == 0 {
                    break;
                }
                lns.step(&mut rng, 4);
                check_no_three_in_line(n, &lns.points()).unwrap();
            }
            assert_eq!(lns.deficit, 0, "n={n}");
            check_solution(n, &lns.points()).unwrap();
        }
    }
}
