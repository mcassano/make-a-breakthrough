//! Grid geometry: cells, collinearity, and enumeration of grid lines.

/// A cell of the `n x n` grid, addressed by `(row, col)` with `0 <= row, col < n`.
pub type Cell = (usize, usize);

pub fn gcd(mut a: i64, mut b: i64) -> i64 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// Canonical key of the line through two distinct cells.
///
/// The direction is reduced to a primitive vector with a positive first
/// non-zero component, and the line is identified by that direction together
/// with the invariant `dc * r - dr * c`.
pub fn line_key(p: Cell, q: Cell) -> (i64, i64, i64) {
    assert_ne!(p, q, "a line needs two distinct cells");
    let (pr, pc) = (p.0 as i64, p.1 as i64);
    let (mut dr, mut dc) = (q.0 as i64 - pr, q.1 as i64 - pc);
    let g = gcd(dr, dc);
    dr /= g;
    dc /= g;
    if dr < 0 || (dr == 0 && dc < 0) {
        dr = -dr;
        dc = -dc;
    }
    (dr, dc, dc * pr - dr * pc)
}

/// Returns every grid line of the `n x n` grid that contains at least
/// `min_points` cells, each given as the list of its cells in order.
///
/// Every line is produced exactly once: for each primitive direction `d` we
/// start only from cells whose predecessor `p - d` lies outside the grid.
pub fn lines(n: usize, min_points: usize) -> Vec<Vec<Cell>> {
    assert!(min_points >= 2);
    let n_i = n as i64;
    // A line with k points in direction (dr, dc) spans (k-1)*max(|dr|,|dc|) <= n-1.
    let max_step = (n_i - 1) / (min_points as i64 - 1);
    let inside = |r: i64, c: i64| r >= 0 && r < n_i && c >= 0 && c < n_i;
    let mut out = Vec::new();
    for dr in 0..=max_step {
        for dc in -max_step..=max_step {
            if (dr == 0 && dc <= 0) || gcd(dr, dc) != 1 {
                continue;
            }
            for r in 0..n_i {
                for c in 0..n_i {
                    if inside(r - dr, c - dc) {
                        continue;
                    }
                    let mut cells = Vec::new();
                    let (mut rr, mut cc) = (r, c);
                    while inside(rr, cc) {
                        cells.push((rr as usize, cc as usize));
                        rr += dr;
                        cc += dc;
                    }
                    if cells.len() >= min_points {
                        out.push(cells);
                    }
                }
            }
        }
    }
    out
}

/// Checks that `points` is a valid no-three-in-line configuration on the
/// `n x n` grid: all points are distinct cells of the grid and no three are
/// collinear. Returns a description of the first violation found.
///
/// This check is deliberately independent of the line enumeration used by the
/// search code: it looks at every pair of points directly.
pub fn check_no_three_in_line(n: usize, points: &[Cell]) -> Result<(), String> {
    for &(r, c) in points {
        if r >= n || c >= n {
            return Err(format!("cell ({r}, {c}) is outside the {n}x{n} grid"));
        }
    }
    let mut seen = std::collections::HashMap::new();
    for (i, &p) in points.iter().enumerate() {
        for &q in &points[i + 1..] {
            if p == q {
                return Err(format!("cell {p:?} is used twice"));
            }
            if let Some((a, b)) = seen.insert(line_key(p, q), (p, q)) {
                return Err(format!("{a:?}, {b:?}, {p:?}, {q:?} lie on one line"));
            }
        }
    }
    Ok(())
}

/// Checks the full no-three-in-line property: exactly `2n` points, no three
/// collinear (which also forces exactly two points in every row and column).
pub fn check_solution(n: usize, points: &[Cell]) -> Result<(), String> {
    if points.len() != 2 * n {
        return Err(format!("expected {} points, found {}", 2 * n, points.len()));
    }
    check_no_three_in_line(n, points)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn brute_force_line_count(n: usize, min_points: usize) -> usize {
        let cells: Vec<Cell> = (0..n).flat_map(|r| (0..n).map(move |c| (r, c))).collect();
        let mut counts = std::collections::HashMap::new();
        for (i, &p) in cells.iter().enumerate() {
            for &q in &cells[i + 1..] {
                *counts.entry(line_key(p, q)).or_insert(0usize) += 1;
            }
        }
        // A line with k points contributes k(k-1)/2 pairs.
        let min_pairs = min_points * (min_points - 1) / 2;
        counts.values().filter(|&&pairs| pairs >= min_pairs).count()
    }

    #[test]
    fn line_enumeration_matches_brute_force() {
        for n in 2..=12 {
            for m in 2..=4 {
                assert_eq!(lines(n, m).len(), brute_force_line_count(n, m), "n={n} m={m}");
            }
        }
    }

    #[test]
    fn enumerated_lines_are_collinear_and_maximal() {
        let n = 9;
        for line in lines(n, 3) {
            let k = line_key(line[0], line[1]);
            assert!(line.windows(2).all(|w| line_key(w[0], w[1]) == k));
            let on_line = (0..n)
                .flat_map(|r| (0..n).map(move |c| (r, c)))
                .filter(|&p| p == line[0] || line_key(line[0], p) == k)
                .count();
            assert_eq!(on_line, line.len());
        }
    }

    #[test]
    fn checker_detects_collinear_triples() {
        assert!(check_no_three_in_line(5, &[(0, 0), (1, 2), (2, 4)]).is_err());
        assert!(check_no_three_in_line(5, &[(0, 0), (1, 2), (2, 3)]).is_ok());
        assert!(check_no_three_in_line(5, &[(0, 0), (0, 0)]).is_err());
        assert!(check_no_three_in_line(5, &[(0, 5)]).is_err());
    }

    #[test]
    fn two_by_two_and_classic_solutions() {
        assert!(check_solution(2, &[(0, 0), (0, 1), (1, 0), (1, 1)]).is_ok());
        // Solution for n = 4.
        let pts = [(0, 1), (0, 2), (1, 0), (1, 3), (2, 0), (2, 3), (3, 1), (3, 2)];
        assert!(check_solution(4, &pts).is_ok());
    }
}
