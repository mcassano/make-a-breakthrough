//! CNF (DIMACS) encoding of a [`Model`].
//!
//! Variable `v` of the model is SAT variable `v + 1`. Line constraints
//! (`sum <= 2`) over few variables are encoded directly by forbidding every
//! minimal overweight combination; longer ones use a sequential counter.
//! Row constraints (`sum == 2`) use the same encodings plus direct "at least
//! two" clauses.

use crate::symmetry::{Model, Sum};
use std::io::{self, Write};

pub type Lit = i32;

pub struct Cnf {
    pub num_vars: u32,
    pub clauses: Vec<Vec<Lit>>,
}

impl Cnf {
    pub fn new(num_vars: u32) -> Cnf {
        Cnf {
            num_vars,
            clauses: Vec::new(),
        }
    }

    pub fn new_var(&mut self) -> Lit {
        self.num_vars += 1;
        self.num_vars as Lit
    }

    pub fn add(&mut self, clause: Vec<Lit>) {
        self.clauses.push(clause);
    }

    /// Encodes `sum(inputs) >= k` with a sequential counter (inputs may
    /// repeat a literal to give it a larger weight).
    pub fn at_least(&mut self, inputs: &[Lit], k: usize) {
        if k == 0 {
            return;
        }
        if inputs.len() < k {
            self.add(vec![]);
            return;
        }
        // r[j] after input i: "at least j+1 of the first i+1 inputs are true",
        // with the implications pointing downwards (r[j] -> ...), which is
        // what forcing the final r[k-1] needs.
        let mut prev: Vec<Lit> = Vec::new();
        for (i, &x) in inputs.iter().enumerate() {
            let width = (i + 1).min(k);
            let cur: Vec<Lit> = (0..width).map(|_| self.new_var()).collect();
            for j in 0..width {
                // cur[j] -> prev[j] or (x and prev[j-1])
                let mut c = vec![-cur[j]];
                if j < prev.len() {
                    c.push(prev[j]);
                }
                c.push(x);
                self.add(c);
                if j > 0 {
                    let mut c = vec![-cur[j]];
                    if j < prev.len() {
                        c.push(prev[j]);
                    }
                    c.push(prev[j - 1]);
                    self.add(c);
                }
            }
            prev = cur;
        }
        self.add(vec![prev[k - 1]]);
    }

    /// Encodes `sum(inputs) <= k` with a sequential counter (inputs may
    /// repeat a literal to give it a larger weight).
    pub fn at_most(&mut self, inputs: &[Lit], k: usize) {
        if inputs.len() <= k {
            return;
        }
        if k == 0 {
            for &x in inputs {
                self.add(vec![-x]);
            }
            return;
        }
        // s[j] after input i: "at least j+1 of the first i+1 inputs are true".
        let mut prev: Vec<Lit> = Vec::new();
        for (i, &x) in inputs.iter().enumerate() {
            if !prev.is_empty() && prev.len() == k {
                self.add(vec![-x, -prev[k - 1]]);
            }
            if i + 1 == inputs.len() {
                break;
            }
            let width = (i + 1).min(k);
            let cur: Vec<Lit> = (0..width).map(|_| self.new_var()).collect();
            for j in 0..width {
                if j < prev.len() {
                    self.add(vec![-prev[j], cur[j]]);
                }
                if j == 0 {
                    self.add(vec![-x, cur[0]]);
                } else if j - 1 < prev.len() {
                    self.add(vec![-x, -prev[j - 1], cur[j]]);
                }
            }
            prev = cur;
        }
    }

    pub fn write_dimacs<W: Write>(&self, mut out: W, comments: &[String]) -> io::Result<()> {
        for c in comments {
            writeln!(out, "c {c}")?;
        }
        writeln!(out, "p cnf {} {}", self.num_vars, self.clauses.len())?;
        let mut line = String::new();
        for clause in &self.clauses {
            line.clear();
            for lit in clause {
                line.push_str(&lit.to_string());
                line.push(' ');
            }
            line.push('0');
            writeln!(out, "{line}")?;
        }
        Ok(())
    }
}

fn lit(var: usize) -> Lit {
    var as Lit + 1
}

/// Encoding choices.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Constraints with at most this many variables are encoded directly
    /// (one clause per forbidden pair or triple) instead of by a counter.
    pub direct_max_vars: usize,
    /// Add clauses that remove the mirror image of each solution.
    pub break_mirror: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            direct_max_vars: 6,
            break_mirror: true,
        }
    }
}

/// Encodes `sum <= 2`.
fn at_most_two(cnf: &mut Cnf, sum: &Sum, opts: &Options) {
    if sum.terms.len() <= opts.direct_max_vars {
        at_most_two_direct(cnf, sum);
    } else {
        at_most_two_counter(cnf, sum);
    }
}

/// Encodes `sum <= 2` by forbidding each minimal set of terms of weight >= 3.
fn at_most_two_direct(cnf: &mut Cnf, sum: &Sum) {
    let t = &sum.terms;
    for (i, &(a, wa)) in t.iter().enumerate() {
        if wa > 2 {
            cnf.add(vec![-lit(a)]);
            continue;
        }
        for (j, &(b, wb)) in t.iter().enumerate().skip(i + 1) {
            if wb > 2 {
                continue;
            }
            if wa + wb > 2 {
                cnf.add(vec![-lit(a), -lit(b)]);
                continue;
            }
            for &(c, wc) in &t[j + 1..] {
                if wc == 1 {
                    cnf.add(vec![-lit(a), -lit(b), -lit(c)]);
                }
            }
        }
    }
}

/// Encodes `sum <= 2` with a sequential counter over the terms.
///
/// `s1[i]` / `s2[i]` mean "the weighted sum of the first `i + 1` terms is at
/// least 1 / at least 2". Only the upward implications are needed.
fn at_most_two_counter(cnf: &mut Cnf, sum: &Sum) {
    // Terms of weight > 2 can never be chosen.
    for &(v, w) in &sum.terms {
        if w > 2 {
            cnf.add(vec![-lit(v)]);
        }
    }
    let terms: Vec<(usize, u32)> = sum.terms.iter().copied().filter(|&(_, w)| w <= 2).collect();
    if terms.iter().map(|&(_, w)| w).sum::<u32>() <= 2 {
        return;
    }
    let k = terms.len();
    let mut prev: Option<(Lit, Lit)> = None;
    for (i, &(v, w)) in terms.iter().enumerate() {
        let x = lit(v);
        let last = i + 1 == k;
        // Forbid overflow past 2.
        if let Some((p1, p2)) = prev {
            let blocker = if w == 2 { p1 } else { p2 };
            cnf.add(vec![-x, -blocker]);
        }
        if last {
            break;
        }
        let (s1, s2) = (cnf.new_var(), cnf.new_var());
        cnf.add(vec![-x, s1]);
        if w == 2 {
            cnf.add(vec![-x, s2]);
        }
        if let Some((p1, p2)) = prev {
            cnf.add(vec![-p1, s1]);
            cnf.add(vec![-p2, s2]);
            if w == 1 {
                cnf.add(vec![-p1, -x, s2]);
            }
        }
        prev = Some((s1, s2));
    }
}

/// Encodes `sum >= 2`: some term is true, and a weight-1 term is never the
/// only true term.
fn at_least_two(cnf: &mut Cnf, sum: &Sum) {
    let all: Vec<Lit> = sum.terms.iter().map(|&(v, _)| lit(v)).collect();
    cnf.add(all.clone());
    for &(v, w) in &sum.terms {
        if w == 1 {
            let mut clause = vec![-lit(v)];
            clause.extend(all.iter().copied().filter(|&l| l != lit(v)));
            cnf.add(clause);
        }
    }
}

pub fn encode(model: &Model, opts: &Options) -> Cnf {
    let mut cnf = Cnf {
        num_vars: model.vars.len() as u32,
        clauses: Vec::new(),
    };
    for row in &model.rows {
        at_most_two(&mut cnf, row, opts);
        at_least_two(&mut cnf, row);
    }
    for line in &model.lines {
        at_most_two(&mut cnf, line, opts);
    }
    if opts.break_mirror {
        break_mirror(&mut cnf, model);
    }
    for group in &model.exactly_one {
        cnf.add(group.iter().map(|&v| lit(v)).collect());
        for (i, &a) in group.iter().enumerate() {
            for &b in &group[i + 1..] {
                cnf.add(vec![-lit(a), -lit(b)]);
            }
        }
    }
    cnf
}

/// Removes one of each pair of mirror-image solutions.
///
/// Let the top row hold points in columns `c1 < c2`. For every class handled
/// here some reflection `R` of the grid maps the class to itself and sends
/// the top row to a row with columns `{n-1-c2, n-1-c1}`: the vertical-axis
/// reflection for `iden`, `rot2` and `rot4`; for `rct4` (whose diagonal pair
/// lies on the main diagonal) the transpose, which maps the top row to the
/// left column, and the quarter turns relate the left column to the top row.
/// Of `C` and `R(C)` at least one has `c1 + c2 <= n - 1`, so requiring that
/// loses no solution up to symmetry. For `rct4` a top-row point in the corner
/// `(0, 0)` gives `c1 = 0` and the constraint holds for both.
fn break_mirror(cnf: &mut Cnf, model: &Model) {
    let n = model.n;
    for c1 in 0..n {
        for c2 in c1 + 1..n {
            if c1 + c2 <= n - 1 {
                continue;
            }
            if let (Some(a), Some(b)) = (model.cell_var[c1], model.cell_var[c2]) {
                if a != b {
                    cnf.add(vec![-lit(a), -lit(b)]);
                }
            }
        }
    }
}

/// Extracts the true model variables from a solver's `v ...` output lines.
pub fn parse_solver_output(model: &Model, output: &str) -> Result<Vec<usize>, String> {
    if !output.lines().any(|l| l.trim() == "s SATISFIABLE") {
        return Err("solver output does not report SATISFIABLE".into());
    }
    let mut true_vars = Vec::new();
    for line in output.lines().filter(|l| l.starts_with("v ")) {
        for tok in line[2..].split_whitespace() {
            let x: i64 = tok.parse().map_err(|_| format!("bad literal {tok:?}"))?;
            if x > 0 && (x as usize) <= model.vars.len() {
                true_vars.push(x as usize - 1);
            }
        }
    }
    Ok(true_vars)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symmetry::Sum;

    /// Brute-force check of the cardinality encodings on all assignments.
    fn check_encoding(weights: &[u32], exact: bool, direct: bool) {
        let k = weights.len();
        let sum = Sum {
            terms: weights.iter().enumerate().map(|(v, &w)| (v, w)).collect(),
        };
        let mut cnf = Cnf {
            num_vars: k as u32,
            clauses: Vec::new(),
        };
        at_most_two(
            &mut cnf,
            &sum,
            &Options {
                direct_max_vars: if direct { k } else { 0 },
                break_mirror: false,
            },
        );
        if exact {
            at_least_two(&mut cnf, &sum);
        }
        let aux = cnf.num_vars as usize - k;
        for mask in 0u32..(1 << k) {
            let total: u32 = (0..k).filter(|&i| mask >> i & 1 == 1).map(|i| weights[i]).sum();
            let expected = if exact { total == 2 } else { total <= 2 };
            // Satisfiable iff some assignment of the auxiliary variables works.
            let sat = (0u64..(1 << aux)).any(|aux_mask| {
                let value = |l: Lit| {
                    let v = l.unsigned_abs() as usize - 1;
                    let b = if v < k {
                        mask >> v & 1 == 1
                    } else {
                        aux_mask >> (v - k) & 1 == 1
                    };
                    if l > 0 {
                        b
                    } else {
                        !b
                    }
                };
                cnf.clauses.iter().all(|c| c.iter().any(|&l| value(l)))
            });
            assert_eq!(sat, expected, "weights {weights:?} mask {mask:b} exact {exact}");
        }
    }

    #[test]
    fn cardinality_encodings_are_exact() {
        for weights in [
            vec![1, 1, 1],
            vec![1, 1, 1, 1, 1],
            vec![2, 1, 1],
            vec![1, 2, 1, 1],
            vec![1, 1, 2],
            vec![2, 2, 1],
            vec![1, 1, 1, 2, 1],
            vec![3, 1, 1],
        ] {
            for direct in [false, true] {
                check_encoding(&weights, false, direct);
                check_encoding(&weights, true, direct);
            }
        }
    }
}

#[cfg(test)]
mod counter_tests {
    use super::*;

    fn satisfiable_with(cnf: &Cnf, inputs: usize, mask: u32) -> bool {
        let aux = cnf.num_vars as usize - inputs;
        (0u64..(1 << aux)).any(|aux_mask| {
            let value = |l: Lit| {
                let v = l.unsigned_abs() as usize - 1;
                let b = if v < inputs {
                    mask >> v & 1 == 1
                } else {
                    aux_mask >> (v - inputs) & 1 == 1
                };
                if l > 0 {
                    b
                } else {
                    !b
                }
            };
            cnf.clauses.iter().all(|c| c.iter().any(|&l| value(l)))
        })
    }

    #[test]
    fn sequential_counters_are_exact() {
        for n in 1..=4usize {
            // Literal 1 appears twice, i.e. has weight 2.
            let mut inputs: Vec<Lit> = (1..=n as Lit).collect();
            inputs.push(1);
            for k in 0..=n + 2 {
                for at_least in [false, true] {
                    let mut cnf = Cnf::new(n as u32);
                    if at_least {
                        cnf.at_least(&inputs, k)
                    } else {
                        cnf.at_most(&inputs, k)
                    }
                    for mask in 0u32..(1 << n) {
                        let total = (0..n).filter(|&i| mask >> i & 1 == 1).count() + (mask & 1) as usize;
                        let expected = if at_least { total >= k } else { total <= k };
                        assert_eq!(
                            satisfiable_with(&cnf, n, mask),
                            expected,
                            "n={n} k={k} at_least={at_least} mask={mask:b}"
                        );
                    }
                }
            }
        }
    }
}
