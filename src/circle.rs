//! The no-four-on-a-circle problem of Erdős and Purdy: choose as many points
//! of the `n x n` grid as possible so that no four lie on a common circle or
//! a common line.
//!
//! Circles are handled exactly: the circle (or line) through three points
//! `(x_i, y_i)` is `A (x^2 + y^2) + D x + E y + F = 0`, where `(A, D, E, F)`
//! are 3x3 determinants of the points. Dividing by the gcd and fixing the
//! sign gives a canonical integer key; `A = 0` exactly when the three points
//! are collinear, in which case the key describes their line. So "no four
//! points on a circle or line" is the same as "no key shared by four
//! points", and both kinds of constraints are handled uniformly.

use crate::geometry::gcd;
use std::collections::HashMap;

pub type Point = (i64, i64);

/// Canonical key of the circle or line through three distinct points.
pub fn key3(p: Point, q: Point, r: Point) -> [i64; 4] {
    let s = |(x, y): Point| x * x + y * y;
    let det = |a: [i64; 3], b: [i64; 3], c: [i64; 3]| {
        a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
    };
    let (sp, sq, sr) = (s(p), s(q), s(r));
    let a = det([p.0, p.1, 1], [q.0, q.1, 1], [r.0, r.1, 1]);
    let d = -det([sp, p.1, 1], [sq, q.1, 1], [sr, r.1, 1]);
    let e = det([sp, p.0, 1], [sq, q.0, 1], [sr, r.0, 1]);
    let f = -det([sp, p.0, p.1], [sq, q.0, q.1], [sr, r.0, r.1]);
    normalize([a, d, e, f])
}

fn normalize(mut k: [i64; 4]) -> [i64; 4] {
    let g = k.iter().fold(0, |g, &v| gcd(g, v));
    if g == 0 {
        return k;
    }
    let first = *k.iter().find(|&&v| v != 0).unwrap();
    let sign = if first < 0 { -1 } else { 1 };
    for v in &mut k {
        *v = *v / g * sign;
    }
    k
}

/// Every circle or line containing at least four points of the `n x n`
/// grid, as the sorted list of point indices (`x * n + y`) on it.
///
/// Each curve is found from its lexicographically smallest point `p`: for
/// every pair `q < r` of larger points the key of `(p, q, r)` is computed,
/// and keys shared by at least three such pairs' points are curves with at
/// least four grid points. A global set removes curves seen from a smaller
/// pivot (whose own point list would then be incomplete).
pub fn curves(n: usize) -> Vec<Vec<u32>> {
    let pts: Vec<Point> = (0..n as i64).flat_map(|x| (0..n as i64).map(move |y| (x, y))).collect();
    let mut seen: std::collections::HashSet<[i64; 4]> = std::collections::HashSet::new();
    let mut out = Vec::new();
    let mut local: HashMap<[i64; 4], Vec<u32>> = HashMap::new();
    for (i, &p) in pts.iter().enumerate() {
        local.clear();
        for j in i + 1..pts.len() {
            for k in j + 1..pts.len() {
                let key = key3(p, pts[j], pts[k]);
                let e = local.entry(key).or_default();
                e.push(j as u32);
                e.push(k as u32);
            }
        }
        for (key, mut others) in local.drain() {
            others.sort_unstable();
            others.dedup();
            if others.len() >= 3 && seen.insert(key) {
                let mut curve = vec![i as u32];
                curve.extend(others);
                out.push(curve);
            }
        }
    }
    out
}

/// Checks that no four of `points` lie on a common circle or line, and that
/// all points are distinct cells of the `n x n` grid. Independent of
/// [`curves`]: it examines every triple of the given points directly.
pub fn check(n: usize, points: &[Point]) -> Result<(), String> {
    for &(x, y) in points {
        if x < 0 || y < 0 || x >= n as i64 || y >= n as i64 {
            return Err(format!("({x}, {y}) is outside the {n}x{n} grid"));
        }
    }
    let mut sorted = points.to_vec();
    sorted.sort_unstable();
    if sorted.windows(2).any(|w| w[0] == w[1]) {
        return Err("repeated point".into());
    }
    // A curve through k of the points is hit by C(k, 3) triples; four points
    // on one curve give four triples with the same key.
    let mut count: HashMap<[i64; 4], u32> = HashMap::new();
    for i in 0..points.len() {
        for j in i + 1..points.len() {
            for k in j + 1..points.len() {
                let c = count.entry(key3(points[i], points[j], points[k])).or_insert(0);
                *c += 1;
                if *c >= 4 {
                    return Err(format!(
                        "four points on the curve through {:?}, {:?}, {:?}",
                        points[i], points[j], points[k]
                    ));
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn brute_force_curve_count(n: usize) -> usize {
        // Count 4-point subsets on a common curve and compare with sum C(k, 4).
        let pts: Vec<Point> = (0..n as i64).flat_map(|x| (0..n as i64).map(move |y| (x, y))).collect();
        let mut quads = 0;
        for a in 0..pts.len() {
            for b in a + 1..pts.len() {
                for c in b + 1..pts.len() {
                    let key = key3(pts[a], pts[b], pts[c]);
                    for d in c + 1..pts.len() {
                        if key3(pts[a], pts[b], pts[d]) == key {
                            quads += 1;
                        }
                    }
                }
            }
        }
        quads
    }

    #[test]
    fn keys_identify_circles_and_lines() {
        // Four points of a unit square are concyclic.
        assert_eq!(key3((0, 0), (1, 0), (0, 1)), key3((1, 1), (1, 0), (0, 1)));
        // Collinear points share a line key with A = 0.
        let k = key3((0, 0), (1, 1), (2, 2));
        assert_eq!(k[0], 0);
        assert_eq!(k, key3((3, 3), (1, 1), (5, 5)));
        assert_ne!(key3((0, 0), (1, 0), (0, 1)), key3((0, 0), (2, 0), (0, 1)));
    }

    #[test]
    fn curve_enumeration_matches_brute_force() {
        for n in 2..=7 {
            let quads: usize = curves(n)
                .iter()
                .map(|c| {
                    let k = c.len();
                    k * (k - 1) * (k - 2) * (k - 3) / 24
                })
                .sum();
            assert_eq!(quads, brute_force_curve_count(n), "n={n}");
        }
    }

    #[test]
    fn checker_accepts_published_example_and_rejects_squares() {
        // The OEIS A395821 example for n = 3.
        assert!(check(3, &[(0, 2), (1, 1), (1, 2), (1, 0), (2, 0)]).is_ok());
        assert!(check(3, &[(0, 0), (0, 1), (1, 0), (1, 1)]).is_err());
        assert!(check(4, &[(0, 0), (1, 1), (2, 2), (3, 3)]).is_err());
    }
}

/// Symmetry imposed on configurations. Only the half turn is useful: any
/// quarter-turn orbit is four concyclic points, and a reflection maps two
/// points to an isosceles trapezium, which is cyclic.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sym {
    None,
    Rot2,
}

/// The SAT model: orbit variables and weighted "at most 3" constraints.
pub struct Model {
    pub n: usize,
    pub sym: Sym,
    /// Grid points of each orbit (variable `i` is SAT variable `i + 1`).
    pub orbits: Vec<Vec<Point>>,
    /// Deduplicated constraints `sum weight * x <= 3` as sorted
    /// `(orbit, weight)` lists; only those with total weight >= 4.
    pub curves: Vec<Vec<(u32, u32)>>,
}

impl Model {
    pub fn new(n: usize, sym: Sym) -> Model {
        let mut orbit_of = vec![u32::MAX; n * n];
        let mut orbits: Vec<Vec<Point>> = Vec::new();
        for x in 0..n {
            for y in 0..n {
                if orbit_of[x * n + y] != u32::MAX {
                    continue;
                }
                let mut o = vec![(x as i64, y as i64)];
                if sym == Sym::Rot2 {
                    let (xx, yy) = (n - 1 - x, n - 1 - y);
                    if (xx, yy) != (x, y) {
                        o.push((xx as i64, yy as i64));
                    }
                }
                for &(a, b) in &o {
                    orbit_of[a as usize * n + b as usize] = orbits.len() as u32;
                }
                orbits.push(o);
            }
        }
        let mut seen = std::collections::HashSet::new();
        let mut out = Vec::new();
        for curve in curves(n) {
            let mut w: Vec<(u32, u32)> = Vec::new();
            for &i in &curve {
                let o = orbit_of[i as usize];
                match w.iter_mut().find(|t| t.0 == o) {
                    Some(t) => t.1 += 1,
                    None => w.push((o, 1)),
                }
            }
            w.sort_unstable();
            if w.iter().map(|t| t.1).sum::<u32>() >= 4 && seen.insert(w.clone()) {
                out.push(w);
            }
        }
        Model {
            n,
            sym,
            orbits,
            curves: out,
        }
    }

    /// CNF for "a valid configuration with at least `k` points".
    pub fn cnf(&self, k: usize) -> crate::cnf::Cnf {
        use crate::cnf::{Cnf, Lit};
        let mut cnf = Cnf::new(self.orbits.len() as u32);
        let lit = |o: u32| o as Lit + 1;
        for curve in &self.curves {
            // Orbits with weight >= 4 on one curve are forbidden outright.
            let terms: Vec<(u32, u32)> = curve
                .iter()
                .copied()
                .filter(|&(o, w)| {
                    if w >= 4 {
                        cnf.add(vec![-lit(o)]);
                    }
                    w < 4
                })
                .collect();
            if terms.iter().map(|t| t.1).sum::<u32>() < 4 {
                continue;
            }
            if terms.len() <= 8 {
                // Forbid every minimal set of total weight >= 4.
                let t = terms.len();
                for mask in 1u32..(1 << t) {
                    let chosen: Vec<usize> = (0..t).filter(|&i| mask >> i & 1 == 1).collect();
                    let weight: u32 = chosen.iter().map(|&i| terms[i].1).sum();
                    let minimal = chosen.iter().all(|&i| weight - terms[i].1 < 4);
                    if weight >= 4 && minimal {
                        cnf.add(chosen.iter().map(|&i| -lit(terms[i].0)).collect());
                    }
                }
            } else {
                let inputs: Vec<Lit> = terms
                    .iter()
                    .flat_map(|&(o, w)| std::iter::repeat(lit(o)).take(w as usize))
                    .collect();
                cnf.at_most(&inputs, 3);
            }
        }
        let inputs: Vec<Lit> = (0..self.orbits.len() as u32)
            .flat_map(|o| std::iter::repeat(lit(o)).take(self.orbits[o as usize].len()))
            .collect();
        cnf.at_least(&inputs, k);
        cnf
    }
}
