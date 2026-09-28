//! Exhaustive and heuristic search tools for three extremal problems on the
//! `n x n` grid:
//!
//! * distinct distances (Erdős–Guy): [`distinct`];
//! * no three in line: [`geometry`], [`symmetry`], [`cnf`], [`dfs`],
//!   [`lns`], [`local`], [`format`];
//! * no four on a circle (Erdős–Purdy): [`circle`], [`circle_ls`].

// Bitset code indexes several parallel word arrays in one loop.
#![allow(clippy::needless_range_loop, clippy::comparison_chain)]

pub mod circle;
pub mod circle_ls;
pub mod cnf;
pub mod dfs;
pub mod distinct;
pub mod format;
pub mod geometry;
pub mod lns;
pub mod local;
pub mod symmetry;
