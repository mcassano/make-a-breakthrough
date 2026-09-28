use no3::cnf;
use no3::dfs::Rot4Dfs;
use no3::format;
use no3::geometry::check_solution;
use no3::lns::Lns;
use no3::local::{anneal, Rng, Search};
use no3::symmetry::{Class, Model};
use std::fs;
use std::io::{self, BufWriter};
use std::process::ExitCode;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

const USAGE: &str = "\
usage:
  no3 cnf <n> <class> [out.cnf] [--direct K] [--no-mirror]
                                        write the SAT encoding (class: iden, rot2, rot4, rct4)
  no3 extract <n> <class> <solver.out>  decode a SAT solver's model, verify it, print its code
  no3 verify <file|code>...             verify encoded configurations (one per line in files)
  no3 show <code>                       draw an encoded configuration
  no3 anneal <n> <class> [--threads T] [--steps S] [--t0 X] [--t1 Y] [--seed K] [--count C]
                                        simulated annealing (class: rot4 for even n, rct4 for odd n)
  no3 dfs <n> [--threads T] [--budget B] [--seed K] [--count C]
                                        randomized restarting tree search for rot4 solutions (even n)
  no3 lns <n> [--threads T] [--k K] [--budget B] [--seed S] [--count C] [--restart R]
                                        large neighbourhood search for rot4 solutions (even n)
  no3 count-rot4 <n>                    count all rot4 solutions (labelled, small n only)";

fn parse_n(s: &str) -> Result<usize, String> {
    let n: usize = s.parse().map_err(|_| format!("invalid grid size {s:?}"))?;
    if n > format::max_n() {
        return Err(format!(
            "grid size {n} exceeds the encodable maximum {}",
            format::max_n()
        ));
    }
    Ok(n)
}

fn cmd_cnf(args: &[String]) -> Result<(), String> {
    let [n, class, rest @ ..] = args else {
        return Err(USAGE.into());
    };
    let model = Model::new(parse_n(n)?, Class::parse(class)?)?;
    let opts = cnf::Options {
        direct_max_vars: flag(rest, "--direct", cnf::Options::default().direct_max_vars)?,
        break_mirror: !rest.iter().any(|a| a == "--no-mirror"),
    };
    let rest: Vec<&String> = rest
        .iter()
        .filter(|a| !a.starts_with("--") && a.parse::<usize>().is_err())
        .collect();
    let cnf = cnf::encode(&model, &opts);
    let comments = vec![
        format!("no-three-in-line n={} class={}", model.n, model.class.name()),
        format!(
            "{} orbit variables, {} row constraints, {} line constraints",
            model.vars.len(),
            model.rows.len(),
            model.lines.len()
        ),
    ];
    eprintln!("{}; {} vars, {} clauses", comments[1], cnf.num_vars, cnf.clauses.len());
    let res = match rest.first() {
        Some(path) => cnf.write_dimacs(
            BufWriter::new(fs::File::create(path).map_err(|e| e.to_string())?),
            &comments,
        ),
        None => cnf.write_dimacs(BufWriter::new(io::stdout().lock()), &comments),
    };
    res.map_err(|e| e.to_string())
}

fn cmd_extract(args: &[String]) -> Result<(), String> {
    let [n, class, path] = args else {
        return Err(USAGE.into());
    };
    let model = Model::new(parse_n(n)?, Class::parse(class)?)?;
    let output = fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let true_vars = cnf::parse_solver_output(&model, &output)?;
    let points = model.points(&true_vars);
    check_solution(model.n, &points)?;
    println!("{}", format::encode(model.n, model.class.code_char(), &points)?);
    Ok(())
}

fn verify_code(code: &str) -> Result<usize, String> {
    let (_, n, points) = format::decode(code)?;
    check_solution(n, &points)?;
    Ok(n)
}

fn cmd_verify(args: &[String]) -> Result<(), String> {
    if args.is_empty() {
        return Err(USAGE.into());
    }
    let mut failures = 0;
    for arg in args {
        let codes: Vec<String> = match fs::read_to_string(arg) {
            Ok(text) => text
                .lines()
                .map(str::trim)
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .map(String::from)
                .collect(),
            Err(_) => vec![arg.clone()],
        };
        for code in codes {
            match verify_code(&code) {
                Ok(n) => println!("valid   n={n:<3} {code}"),
                Err(e) => {
                    failures += 1;
                    println!("INVALID {e}: {code}");
                }
            }
        }
    }
    if failures > 0 {
        return Err(format!("{failures} invalid configuration(s)"));
    }
    Ok(())
}

fn cmd_show(args: &[String]) -> Result<(), String> {
    let [code] = args else { return Err(USAGE.into()) };
    let (_, n, points) = format::decode(code)?;
    let mut grid = vec![vec!['.'; n]; n];
    for (r, c) in points {
        grid[r][c] = '#';
    }
    for row in grid {
        println!("{}", row.iter().map(|ch| ch.to_string()).collect::<Vec<_>>().join(" "));
    }
    Ok(())
}

fn flag<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> Result<T, String> {
    match args.iter().position(|a| a == name) {
        Some(i) => args
            .get(i + 1)
            .and_then(|v| v.parse().ok())
            .ok_or_else(|| format!("invalid or missing value for {name}")),
        None => Ok(default),
    }
}

fn cmd_anneal(args: &[String]) -> Result<(), String> {
    let [n, class, opts @ ..] = args else {
        return Err(USAGE.into());
    };
    let (n, class) = (parse_n(n)?, Class::parse(class)?);
    let default_threads = std::thread::available_parallelism().map_or(1, |p| p.get());
    let threads: usize = flag(opts, "--threads", default_threads)?;
    let steps: u64 = flag(opts, "--steps", 2_000_000)?;
    let t0: f64 = flag(opts, "--t0", 4.0)?;
    let t1: f64 = flag(opts, "--t1", 0.5)?;
    let seed: u64 = flag(opts, "--seed", 0)?;
    let count: usize = flag(opts, "--count", 1)?;
    Search::new(n, class, &mut Rng::new(0))?; // validate parameters early

    let start = Instant::now();
    let restarts = AtomicU64::new(0);
    let found = AtomicUsize::new(0);
    let best = AtomicU64::new(u64::MAX);
    let out = Mutex::new(io::stdout());
    std::thread::scope(|scope| {
        for t in 0..threads {
            let (restarts, found, best, out) = (&restarts, &found, &best, &out);
            scope.spawn(move || {
                let mut rng = Rng::new(seed.wrapping_mul(1_000_003).wrapping_add(t as u64 + 1));
                while found.load(Ordering::Relaxed) < count {
                    let mut search = Search::new(n, class, &mut rng).expect("validated");
                    let result = anneal(&mut search, &mut rng, steps, t0, t1);
                    let r = restarts.fetch_add(1, Ordering::Relaxed) + 1;
                    let prev_best = best.fetch_min(search.cost as u64, Ordering::Relaxed);
                    let secs = start.elapsed().as_secs_f64();
                    if let Some(points) = result {
                        check_solution(n, &points).expect("search returned an invalid configuration");
                        let code = format::encode(n, class.code_char(), &points).expect("encodable");
                        let mut out = out.lock().unwrap();
                        use std::io::Write;
                        writeln!(out, "{code}").ok();
                        out.flush().ok();
                        eprintln!("[{secs:.1}s] solution after {r} restarts");
                        found.fetch_add(1, Ordering::Relaxed);
                    } else if (search.cost as u64) < prev_best || r % 1000 == 0 {
                        eprintln!(
                            "[{secs:.1}s] restarts {r}, best remaining triples {}",
                            best.load(Ordering::Relaxed) / 2
                        );
                    }
                }
            });
        }
    });
    Ok(())
}

fn cmd_dfs(args: &[String]) -> Result<(), String> {
    let [n, opts @ ..] = args else { return Err(USAGE.into()) };
    let n = parse_n(n)?;
    let default_threads = std::thread::available_parallelism().map_or(1, |p| p.get());
    let threads: usize = flag(opts, "--threads", default_threads)?;
    let budget: u64 = flag(opts, "--budget", 1_000_000)?;
    let seed: u64 = flag(opts, "--seed", 0)?;
    let count: usize = flag(opts, "--count", 1)?;
    Rot4Dfs::new(n)?;

    let start = Instant::now();
    let restarts = AtomicU64::new(0);
    let nodes = AtomicU64::new(0);
    let found = AtomicUsize::new(0);
    let out = Mutex::new(io::stdout());
    std::thread::scope(|scope| {
        for t in 0..threads {
            let (restarts, nodes, found, out) = (&restarts, &nodes, &found, &out);
            scope.spawn(move || {
                let mut rng = Rng::new(seed.wrapping_mul(1_000_003).wrapping_add(t as u64 + 1));
                let mut dfs = Rot4Dfs::new(n).expect("validated");
                while found.load(Ordering::Relaxed) < count {
                    let result = dfs.attempt(&mut rng, budget);
                    let r = restarts.fetch_add(1, Ordering::Relaxed) + 1;
                    let total = nodes.fetch_add(dfs.nodes, Ordering::Relaxed) + dfs.nodes;
                    let secs = start.elapsed().as_secs_f64();
                    if let Some(points) = result {
                        check_solution(n, &points).expect("search returned an invalid configuration");
                        let code = format::encode(n, 'o', &points).expect("encodable");
                        let mut out = out.lock().unwrap();
                        use std::io::Write;
                        writeln!(out, "{code}").ok();
                        out.flush().ok();
                        eprintln!("[{secs:.1}s] solution after {r} restarts, {total} nodes");
                        found.fetch_add(1, Ordering::Relaxed);
                    } else if r.is_power_of_two() {
                        eprintln!(
                            "[{secs:.1}s] restarts {r}, {total} nodes, {:.2e} nodes/s",
                            total as f64 / secs
                        );
                    }
                }
            });
        }
    });
    Ok(())
}

fn cmd_lns(args: &[String]) -> Result<(), String> {
    let [n, opts @ ..] = args else { return Err(USAGE.into()) };
    let n = parse_n(n)?;
    let default_threads = std::thread::available_parallelism().map_or(1, |p| p.get());
    let threads: usize = flag(opts, "--threads", default_threads)?;
    let k: usize = flag(opts, "--k", 6)?;
    let budget: u64 = flag(opts, "--budget", 20_000)?;
    let seed: u64 = flag(opts, "--seed", 0)?;
    let count: usize = flag(opts, "--count", 1)?;
    let restart_after: u64 = flag(opts, "--restart", 200_000)?;
    Rot4Dfs::new(n)?;

    let start = Instant::now();
    let found = AtomicUsize::new(0);
    let total_steps = AtomicU64::new(0);
    let total_nodes = AtomicU64::new(0);
    let out = Mutex::new(io::stdout());
    std::thread::scope(|scope| {
        for t in 0..threads {
            let (found, out, total_steps, total_nodes) = (&found, &out, &total_steps, &total_nodes);
            scope.spawn(move || {
                let mut rng = Rng::new(seed.wrapping_mul(1_000_003).wrapping_add(t as u64 + 1));
                while found.load(Ordering::Relaxed) < count {
                    let mut lns = Lns::new(n, &mut rng, budget).expect("validated");
                    let mut best = lns.deficit;
                    let mut since_improvement = 0u64;
                    let mut steps = 0u64;
                    while lns.deficit > 0 && since_improvement < restart_after && found.load(Ordering::Relaxed) < count
                    {
                        let before = lns.nodes;
                        lns.step(&mut rng, k);
                        total_nodes.fetch_add(lns.nodes - before, Ordering::Relaxed);
                        total_steps.fetch_add(1, Ordering::Relaxed);
                        steps += 1;
                        since_improvement += 1;
                        if lns.deficit < best {
                            best = lns.deficit;
                            since_improvement = 0;
                            if t == 0 {
                                eprintln!(
                                    "[{:.1}s] thread 0: deficit {best} after {steps} steps",
                                    start.elapsed().as_secs_f64()
                                );
                            }
                        }
                    }
                    if lns.deficit == 0 {
                        let points = lns.points();
                        check_solution(n, &points).expect("search returned an invalid configuration");
                        let code = format::encode(n, 'o', &points).expect("encodable");
                        let mut out = out.lock().unwrap();
                        use std::io::Write;
                        writeln!(out, "{code}").ok();
                        out.flush().ok();
                        eprintln!(
                            "[{:.1}s] solution; all threads so far: {} steps, {:.3e} repair nodes",
                            start.elapsed().as_secs_f64(),
                            total_steps.load(Ordering::Relaxed),
                            total_nodes.load(Ordering::Relaxed) as f64
                        );
                        found.fetch_add(1, Ordering::Relaxed);
                    }
                }
            });
        }
    });
    Ok(())
}

fn cmd_count_rot4(args: &[String]) -> Result<(), String> {
    let [n] = args else { return Err(USAGE.into()) };
    let mut dfs = Rot4Dfs::new(parse_n(n)?)?;
    println!("{}", dfs.count_all());
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("cnf") => cmd_cnf(&args[1..]),
        Some("extract") => cmd_extract(&args[1..]),
        Some("verify") => cmd_verify(&args[1..]),
        Some("show") => cmd_show(&args[1..]),
        Some("anneal") => cmd_anneal(&args[1..]),
        Some("dfs") => cmd_dfs(&args[1..]),
        Some("count-rot4") => cmd_count_rot4(&args[1..]),
        Some("lns") => cmd_lns(&args[1..]),
        _ => Err(USAGE.into()),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}
