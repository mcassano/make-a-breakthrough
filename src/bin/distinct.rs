//! Decide whether k points of the n x n grid can have pairwise distinct
//! distances, by exhaustive search.

use no3::distinct::{check, normal_form, search, spans, Options, Order};
use std::process::ExitCode;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

const USAGE: &str = "\
usage: distinct <n> <k> [--threads T] [--count] [--classes] [--span] [--mirror] [--boxed | --onion] [--quiet]
  Searches the k-point subsets of the n x n grid with pairwise distinct
  distances, normalised by translation to touch row 0 and column 0.
  --span    only sets that also touch column n-1 (exact when no k-set fits
            in the (n-1) x (n-1) grid; see src/distinct.rs)
  --mirror  additionally keep one set of each mirror pair (needs --span)
  --boxed   search each bounding-box height with rows ordered outside-in and
            break both reflections of the box (implies --span --mirror)
  --onion   like --boxed, but search the boundary ring of the box first and
            break all its symmetries (8 for a square box)
  --count   count all solutions instead of stopping at the first
  --classes also count the spanning solutions up to symmetry (a cross-check
            between orders: all must report the same number)";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (Some(n), Some(k)) = (
        args.first().and_then(|s| s.parse::<usize>().ok()),
        args.get(1).and_then(|s| s.parse::<usize>().ok()),
    ) else {
        eprintln!("{USAGE}");
        return ExitCode::FAILURE;
    };
    let flag = |f: &str| args.iter().any(|a| a == f);
    let threads = args
        .iter()
        .position(|a| a == "--threads")
        .and_then(|i| args.get(i + 1)?.parse().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |p| p.get()));
    let order = if flag("--onion") {
        Order::Onion
    } else if flag("--boxed") {
        Order::BoxedRows
    } else {
        Order::Rows
    };
    let reduced = order != Order::Rows || flag("--mirror");
    let classes = flag("--classes");
    let opts = Options {
        span: flag("--span") || reduced,
        mirror: reduced,
        first_only: !flag("--count") && !classes,
        order,
        collect: classes,
    };
    let quiet = flag("--quiet");
    let t = Instant::now();
    let last = AtomicU64::new(0);
    let out = search(n, k, opts, threads, |h, done, total, nodes, sols| {
        let secs = t.elapsed().as_secs();
        if !quiet && secs >= last.load(Ordering::Relaxed) + 60 {
            last.store(secs, Ordering::Relaxed);
            eprintln!(
                "[{secs}s] box height {h}: {done}/{total} work items, {nodes:.3e} nodes, {sols} solutions",
                nodes = nodes as f64
            );
        }
    });
    let secs = t.elapsed().as_secs_f64();
    if classes {
        let set: std::collections::BTreeSet<_> =
            out.all.iter().filter(|s| spans(n, s)).map(|s| normal_form(s)).collect();
        println!(
            "n={n} k={k}: {} solutions, {} classes of spanning sets under symmetry",
            out.solutions,
            set.len()
        );
    }
    let summary = format!(
        "{:.3e} nodes, {secs:.1}s, {threads} threads, {opts:?}",
        out.nodes as f64
    );
    match &out.example {
        Some(pts) => {
            if let Err(e) = check(n, pts) {
                eprintln!("INTERNAL ERROR: invalid example: {e}");
                return ExitCode::FAILURE;
            }
            let list = pts
                .iter()
                .map(|(x, y)| format!("({x},{y})"))
                .collect::<Vec<_>>()
                .join(" ");
            let count = if opts.first_only {
                String::new()
            } else {
                format!("{} solutions; ", out.solutions)
            };
            println!("n={n} k={k}: FOUND {list} ({count}{summary})");
        }
        None => println!("n={n} k={k}: NONE ({summary})"),
    }
    ExitCode::SUCCESS
}
