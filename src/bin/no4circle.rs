//! Command-line tools for the no-four-on-a-circle problem.

use no3::circle::{check, curves, Model, Point, Sym};
use no3::circle_ls::{Instance, Search};
use no3::local::Rng;
use std::io::BufWriter;
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

const USAGE: &str = "\
usage:
  no4circle curves <n>                   count the circles/lines with >= 4 grid points
  no4circle check <n> <x,y> <x,y> ...    verify a point set
  no4circle check-file <file>            verify lines of the form `n: x,y x,y ...`
  no4circle cnf <n> <k> <none|rot2> <out.cnf>   SAT encoding: >= k points, no four on a circle/line
  no4circle extract <n> <none|rot2> <solver.out>  decode, verify and print a solver's model
  no4circle search <n> [--threads T] [--iters I] [--tabu K] [--seed S] [--target M]
                                         iterated local search; prints each new best set";

fn parse_sym(s: &str) -> Result<Sym, String> {
    match s {
        "none" => Ok(Sym::None),
        "rot2" => Ok(Sym::Rot2),
        _ => Err(format!("unknown symmetry {s:?}")),
    }
}

fn parse_point(s: &str) -> Result<Point, String> {
    let (a, b) = s.split_once(',').ok_or_else(|| format!("bad point {s:?}"))?;
    let x = a.trim().parse().map_err(|_| format!("bad point {s:?}"))?;
    let y = b.trim().parse().map_err(|_| format!("bad point {s:?}"))?;
    Ok((x, y))
}

fn run(args: &[String]) -> Result<(), String> {
    match args.first().map(String::as_str) {
        Some("curves") => {
            let n: usize = args.get(1).and_then(|s| s.parse().ok()).ok_or(USAGE)?;
            let t = Instant::now();
            let cs = curves(n);
            let lines = cs.iter().filter(|c| c.len() >= 4).count();
            let max = cs.iter().map(Vec::len).max().unwrap_or(0);
            let incidences: usize = cs.iter().map(Vec::len).sum();
            println!(
                "n={n}: {lines} curves with >= 4 grid points ({incidences} incidences, max {max} points) in {:.1}s",
                t.elapsed().as_secs_f64()
            );
            Ok(())
        }
        Some("check") => {
            let n: usize = args.get(1).and_then(|s| s.parse().ok()).ok_or(USAGE)?;
            let pts = args[2..]
                .iter()
                .map(|s| parse_point(s))
                .collect::<Result<Vec<_>, _>>()?;
            check(n, &pts)?;
            println!(
                "valid: {} points on the {n}x{n} grid, no four on a circle or line",
                pts.len()
            );
            Ok(())
        }
        Some("check-file") => {
            let path = args.get(1).ok_or(USAGE)?;
            let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
            let mut bad = 0;
            for line in text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#')) {
                let (n, rest) = line.split_once(':').ok_or_else(|| format!("bad line {line:?}"))?;
                let n: usize = n.trim().parse().map_err(|_| format!("bad n in {line:?}"))?;
                let pts = rest
                    .split_whitespace()
                    .map(parse_point)
                    .collect::<Result<Vec<_>, _>>()?;
                match check(n, &pts) {
                    Ok(()) => println!("n={n:<3} {:>3} points  valid", pts.len()),
                    Err(e) => {
                        bad += 1;
                        println!("n={n:<3} {:>3} points  INVALID: {e}", pts.len());
                    }
                }
            }
            if bad > 0 {
                return Err(format!("{bad} invalid configuration(s)"));
            }
            Ok(())
        }
        Some("cnf") => {
            let [_, n, k, sym, out] = args else {
                return Err(USAGE.into());
            };
            let n: usize = n.parse().map_err(|_| USAGE)?;
            let k: usize = k.parse().map_err(|_| USAGE)?;
            let model = Model::new(n, parse_sym(sym)?);
            let cnf = model.cnf(k);
            let file = std::fs::File::create(out).map_err(|e| e.to_string())?;
            let comments = vec![format!(
                "no-four-on-a-circle n={n} k={k} sym={sym}; {} orbits, {} curves",
                model.orbits.len(),
                model.curves.len()
            )];
            cnf.write_dimacs(BufWriter::new(file), &comments)
                .map_err(|e| e.to_string())?;
            eprintln!("{}; {} vars, {} clauses", comments[0], cnf.num_vars, cnf.clauses.len());
            Ok(())
        }
        Some("extract") => {
            let [_, n, sym, path] = args else {
                return Err(USAGE.into());
            };
            let n: usize = n.parse().map_err(|_| USAGE)?;
            let model = Model::new(n, parse_sym(sym)?);
            let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
            if !text.lines().any(|l| l.trim() == "s SATISFIABLE") {
                return Err("solver output does not report SATISFIABLE".into());
            }
            let mut pts: Vec<Point> = Vec::new();
            for line in text.lines().filter(|l| l.starts_with("v ")) {
                for tok in line[2..].split_whitespace() {
                    let v: i64 = tok.parse().map_err(|_| format!("bad literal {tok}"))?;
                    if v > 0 && (v as usize) <= model.orbits.len() {
                        pts.extend(&model.orbits[v as usize - 1]);
                    }
                }
            }
            pts.sort_unstable();
            check(n, &pts)?;
            println!(
                "{n}: {}",
                pts.iter()
                    .map(|(x, y)| format!("{x},{y}"))
                    .collect::<Vec<_>>()
                    .join(" ")
            );
            eprintln!("valid configuration with {} points", pts.len());
            Ok(())
        }
        Some("search") => {
            let n: usize = args.get(1).and_then(|s| s.parse().ok()).ok_or(USAGE)?;
            let opts = &args[2..];
            let get = |name: &str, default: u64| -> Result<u64, String> {
                match opts.iter().position(|a| a == name) {
                    Some(i) => opts
                        .get(i + 1)
                        .and_then(|v| v.parse().ok())
                        .ok_or(format!("bad {name}")),
                    None => Ok(default),
                }
            };
            let threads = get(
                "--threads",
                std::thread::available_parallelism().map_or(1, |p| p.get()) as u64,
            )? as usize;
            let iters = get("--iters", 1_000_000)?;
            let tabu = get("--tabu", 7)?;
            let seed = get("--seed", 0)?;
            let target = get("--target", u64::MAX)? as usize;
            let t0 = Instant::now();
            let inst = Instance::new(n);
            eprintln!(
                "n={n}: {} curves, setup {:.1}s",
                inst.num_curves(),
                t0.elapsed().as_secs_f64()
            );
            let best = AtomicUsize::new(0);
            let out = Mutex::new(());
            std::thread::scope(|scope| {
                for t in 0..threads {
                    let (inst, best, out) = (&inst, &best, &out);
                    scope.spawn(move || {
                        let mut rng = Rng::new(seed.wrapping_mul(7919).wrapping_add(t as u64 + 1));
                        while best.load(Ordering::Relaxed) < target {
                            let mut s = Search::new(inst);
                            s.run(&mut rng, iters, tabu, |st| {
                                if st.size() > best.load(Ordering::Relaxed) {
                                    let _g = out.lock().unwrap();
                                    if st.size() > best.load(Ordering::Relaxed) {
                                        best.store(st.size(), Ordering::Relaxed);
                                        let pts = st.points();
                                        check(n, &pts).expect("local search produced an invalid set");
                                        println!(
                                            "{n}: {}",
                                            pts.iter()
                                                .map(|(x, y)| format!("{x},{y}"))
                                                .collect::<Vec<_>>()
                                                .join(" ")
                                        );
                                        eprintln!(
                                            "[{:.1}s] thread {t}: {} points (iteration {})",
                                            t0.elapsed().as_secs_f64(),
                                            pts.len(),
                                            st.iter
                                        );
                                    }
                                }
                            });
                            if target == usize::MAX {
                                break;
                            }
                        }
                    });
                }
            });
            Ok(())
        }
        _ => Err(USAGE.into()),
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}
