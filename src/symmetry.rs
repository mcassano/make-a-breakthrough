//! Symmetry classes and the reduced 0/1 model they induce.
//!
//! Imposing a symmetry group `G` on a configuration turns the problem on the
//! `n^2` cells into a problem on the `G`-orbits of cells: a configuration is a
//! set of orbits. Every constraint of the original problem ("at most two
//! points on this line", "exactly two points in this row") becomes a weighted
//! constraint over orbit variables, where the weight of an orbit is the number
//! of its cells on the line. Constraints that are images of each other under
//! `G` coincide, so each is kept once.

use crate::geometry::{lines, Cell};
use std::collections::HashSet;

/// The symmetry classes of Flammenkamp's database that this crate can search.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    /// No imposed symmetry.
    Iden,
    /// Invariant under the half turn.
    Rot2,
    /// Invariant under the quarter turn (requires even `n`).
    Rot4,
    /// "Nearly quarter turn" (odd `n`): quarter-turn invariant except for one
    /// half-turn pair on the main diagonal. For odd `n` this is the only way
    /// to combine quarter-turn symmetry with `2n = 2 mod 4` points.
    Rct4,
}

impl Class {
    pub fn parse(s: &str) -> Result<Class, String> {
        match s {
            "iden" => Ok(Class::Iden),
            "rot2" => Ok(Class::Rot2),
            "rot4" => Ok(Class::Rot4),
            "rct4" => Ok(Class::Rct4),
            _ => Err(format!("unknown symmetry class {s:?} (use iden, rot2, rot4 or rct4)")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Class::Iden => "iden",
            Class::Rot2 => "rot2",
            Class::Rot4 => "rot4",
            Class::Rct4 => "rct4",
        }
    }

    /// The symmetry character used in the database encoding, assuming the
    /// configuration has no symmetry beyond the imposed one.
    pub fn code_char(self) -> char {
        match self {
            Class::Iden => '.',
            Class::Rot2 => ':',
            Class::Rot4 => 'o',
            Class::Rct4 => 'c',
        }
    }
}

pub fn rot90(n: usize, (r, c): Cell) -> Cell {
    (c, n - 1 - r)
}

pub fn rot180(n: usize, (r, c): Cell) -> Cell {
    (n - 1 - r, n - 1 - c)
}

/// A weighted "sum of variables" constraint: `sum weight * x_var`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Sum {
    /// `(variable, weight)` pairs sorted by variable, weights >= 1.
    pub terms: Vec<(usize, u32)>,
}

impl Sum {
    pub fn total_weight(&self) -> u32 {
        self.terms.iter().map(|&(_, w)| w).sum()
    }
}

/// The reduced problem for one grid size and symmetry class.
///
/// A solution is a 0/1 assignment to `vars` such that every `rows` sum equals
/// 2, every `lines` sum is at most 2, and each `exactly_one` group has
/// exactly one true variable.
pub struct Model {
    pub n: usize,
    pub class: Class,
    /// The cells of each variable (an orbit, or a special pair for `Rct4`).
    pub vars: Vec<Vec<Cell>>,
    /// Row and column constraints (each must sum to exactly 2).
    pub rows: Vec<Sum>,
    /// Line constraints (each must sum to at most 2), for lines of >= 3 cells.
    pub lines: Vec<Sum>,
    pub exactly_one: Vec<Vec<usize>>,
    /// `cell_var[r * n + c]` is the variable covering the cell, if any.
    pub cell_var: Vec<Option<usize>>,
}

impl Model {
    pub fn new(n: usize, class: Class) -> Result<Model, String> {
        if n < 2 {
            return Err("n must be at least 2".into());
        }
        match class {
            Class::Rot4 if n % 2 == 1 => return Err("rot4 requires even n".into()),
            Class::Rct4 if n % 2 == 0 => return Err("rct4 requires odd n".into()),
            _ => {}
        }
        let mut vars: Vec<Vec<Cell>> = Vec::new();
        let mut cell_var = vec![None; n * n];
        let center = (n % 2 == 1).then_some((n / 2, n / 2));
        let on_diagonal = |(r, c): Cell| r == c || r + c == n - 1;
        for r in 0..n {
            for c in 0..n {
                let p = (r, c);
                if cell_var[r * n + c].is_some() {
                    continue;
                }
                let orbit: Vec<Cell> = match class {
                    Class::Iden => vec![p],
                    Class::Rot2 => {
                        // An odd number of points would be needed on the center.
                        if Some(p) == center {
                            continue;
                        }
                        vec![p, rot180(n, p)]
                    }
                    Class::Rot4 => {
                        let q = rot90(n, p);
                        vec![p, q, rot90(n, q), rot90(n, rot90(n, q))]
                    }
                    Class::Rct4 => {
                        if on_diagonal(p) {
                            // Only the half-turn pairs of the main diagonal (below).
                            if r == c && Some(p) != center {
                                vec![p, rot180(n, p)]
                            } else {
                                continue;
                            }
                        } else {
                            let q = rot90(n, p);
                            vec![p, q, rot90(n, q), rot90(n, rot90(n, q))]
                        }
                    }
                };
                for &(rr, cc) in &orbit {
                    cell_var[rr * n + cc] = Some(vars.len());
                }
                vars.push(orbit);
            }
        }

        let mut exactly_one = Vec::new();
        if class == Class::Rct4 {
            // Exactly two points on the diagonals, forming one half-turn pair
            // on the main diagonal (the anti-diagonal case is its mirror image).
            let group: Vec<usize> = (0..n / 2).map(|a| cell_var[a * n + a].unwrap()).collect();
            exactly_one.push(group);
        }

        let mut model = Model {
            n,
            class,
            vars,
            rows: Vec::new(),
            lines: Vec::new(),
            exactly_one,
            cell_var,
        };
        let mut seen = HashSet::new();
        for i in 0..n {
            let row: Vec<Cell> = (0..n).map(|j| (i, j)).collect();
            let col: Vec<Cell> = (0..n).map(|j| (j, i)).collect();
            for cells in [row, col] {
                let sum = model.sum_over(&cells);
                if seen.insert(sum.clone()) {
                    model.rows.push(sum);
                }
            }
        }
        for line in lines(n, 3) {
            let sum = model.sum_over(&line);
            if sum.total_weight() > 2 && seen.insert(sum.clone()) {
                model.lines.push(sum);
            }
        }
        Ok(model)
    }

    /// The weighted sum over the variables covering `cells`.
    pub fn sum_over(&self, cells: &[Cell]) -> Sum {
        let mut terms: Vec<(usize, u32)> = Vec::new();
        for &(r, c) in cells {
            if let Some(v) = self.cell_var[r * self.n + c] {
                terms.push((v, 1));
            }
        }
        terms.sort_unstable();
        terms.dedup_by(|b, a| {
            if a.0 == b.0 {
                a.1 += b.1;
                true
            } else {
                false
            }
        });
        Sum { terms }
    }

    /// Expands a set of true variables into grid points.
    pub fn points(&self, true_vars: &[usize]) -> Vec<Cell> {
        let mut pts: Vec<Cell> = true_vars.iter().flat_map(|&v| self.vars[v].iter().copied()).collect();
        pts.sort_unstable();
        pts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orbits_partition_the_allowed_cells() {
        for (n, class) in [(10, Class::Rot4), (11, Class::Rct4), (9, Class::Rot2), (6, Class::Iden)] {
            let m = Model::new(n, class).unwrap();
            let covered: usize = m.vars.iter().map(Vec::len).sum();
            let distinct: HashSet<Cell> = m.vars.iter().flatten().copied().collect();
            assert_eq!(covered, distinct.len(), "{class:?}: orbits overlap");
            for (v, cells) in m.vars.iter().enumerate() {
                for &(r, c) in cells {
                    assert_eq!(m.cell_var[r * n + c], Some(v));
                }
            }
        }
    }

    #[test]
    fn rot4_orbits_have_four_cells() {
        let m = Model::new(12, Class::Rot4).unwrap();
        assert_eq!(m.vars.len(), 36);
        assert!(m.vars.iter().all(|o| o.len() == 4));
        // Quarter-turn symmetry makes every column constraint a row constraint.
        assert_eq!(m.rows.len(), 6);
    }

    #[test]
    fn rct4_excludes_center_and_anti_diagonal() {
        let n = 11;
        let m = Model::new(n, Class::Rct4).unwrap();
        assert!(m.cell_var[(n / 2) * n + n / 2].is_none());
        for a in 0..n {
            assert!(a == n - 1 - a || m.cell_var[a * n + (n - 1 - a)].is_none());
        }
        assert_eq!(
            m.exactly_one,
            vec![(0..5).map(|a| m.cell_var[a * n + a].unwrap()).collect::<Vec<_>>()]
        );
    }
}
