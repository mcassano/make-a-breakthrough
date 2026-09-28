//! Knuth estimate of the search-tree size for `distinct n k [span]`.
use no3::distinct::{Grid, Options, Search};
use no3::local::Rng;
fn main() {
    let a: Vec<usize> = std::env::args().skip(1).map(|x| x.parse().unwrap()).collect();
    let (n, k, probes, span) = (a[0], a[1], a[2], a.get(3) == Some(&1));
    let g = Grid::new(n);
    let mut s = Search::new(
        &g,
        k,
        Options {
            span,
            ..Options::default()
        },
    );
    let mut rng = Rng::new(12345);
    let mut sum = 0.0;
    for i in 1..=probes {
        sum += s.knuth_probe(&mut rng);
        if i % (probes / 4).max(1) == 0 {
            println!("after {i} probes: estimated nodes {:.3e}", sum / i as f64);
        }
    }
}
