//! Flammenkamp's text encoding of no-three-in-line configurations.
//!
//! A configuration on the `n x n` grid is written as one symmetry-class
//! character followed by `2n` characters: for each row from top to bottom,
//! the column indices of its two points from left to right. Column indices
//! are written in a 90-character alphabet, which covers `n <= 90`.
//! See <https://wwwhomes.uni-bielefeld.de/achim/no3in/readme.html>.

use crate::geometry::Cell;

pub const ALPHABET: &str = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz#$%&@?!()[]<>{}=*+|-/~^_:;,.";

/// Symmetry-class characters in the order used by the database:
/// iden, rot2, dia1, ort1, rot4, rct4, dia2, ort2, full.
pub const CLASS_CHARS: &str = ".:/-ocx+*";

pub fn max_n() -> usize {
    ALPHABET.len()
}

/// Encodes a full configuration (exactly two points per row).
pub fn encode(n: usize, class_char: char, points: &[Cell]) -> Result<String, String> {
    if n > max_n() {
        return Err(format!("n = {n} exceeds the alphabet size {}", max_n()));
    }
    let mut rows = vec![Vec::new(); n];
    for &(r, c) in points {
        if r >= n || c >= n {
            return Err(format!("cell ({r}, {c}) is outside the grid"));
        }
        rows[r].push(c);
    }
    let alphabet: Vec<char> = ALPHABET.chars().collect();
    let mut out = String::with_capacity(2 * n + 1);
    out.push(class_char);
    for (r, cols) in rows.iter_mut().enumerate() {
        if cols.len() != 2 {
            return Err(format!("row {r} has {} points, expected 2", cols.len()));
        }
        cols.sort_unstable();
        out.extend(cols.iter().map(|&c| alphabet[c]));
    }
    Ok(out)
}

/// Decodes a configuration, returning `(class_char, n, points)`.
pub fn decode(code: &str) -> Result<(char, usize, Vec<Cell>), String> {
    let mut chars = code.trim().chars();
    let class_char = chars.next().ok_or("empty configuration")?;
    let body: Vec<char> = chars.collect();
    if body.is_empty() || body.len() % 2 != 0 {
        return Err(format!("expected an even number of positions, got {}", body.len()));
    }
    let n = body.len() / 2;
    let mut points = Vec::with_capacity(body.len());
    for (i, ch) in body.iter().enumerate() {
        let c = ALPHABET
            .find(*ch)
            .ok_or_else(|| format!("invalid position character {ch:?}"))?;
        if c >= n {
            return Err(format!("column {c} is outside the {n}x{n} grid"));
        }
        points.push((i / 2, c));
    }
    Ok((class_char, n, points))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::check_solution;

    #[test]
    fn alphabet_has_ninety_distinct_characters() {
        let mut chars: Vec<char> = ALPHABET.chars().collect();
        assert_eq!(chars.len(), 90);
        chars.sort_unstable();
        chars.dedup();
        assert_eq!(chars.len(), 90);
    }

    #[test]
    fn decodes_a_published_configuration() {
        // A rot4 solution for n = 8 quoted in the notes of the database.
        let code = "o2323676701014545";
        let (class, n, points) = decode(code).unwrap();
        assert_eq!((class, n), ('o', 8));
        assert!(check_solution(n, &points).is_ok());
        assert_eq!(encode(n, class, &points).unwrap(), code);
    }
}
