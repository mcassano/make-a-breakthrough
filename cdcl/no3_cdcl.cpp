// CDCL search for symmetric no-three-in-line configurations.
//
// The problem is solved by CaDiCaL with an external propagator (IPASIR-UP):
//
//  * One SAT variable per symmetry orbit of grid cells (see src/symmetry.rs
//    and src/local.rs in this repository for the derivation of the model).
//  * Row constraints ("exactly two points per row") and the conflicts
//    between an orbit and the orbits on lines through two of its own points
//    are given to the solver as ordinary clauses.
//  * The "at most two points per line" constraints are never written down.
//    Whenever two chosen orbits span a line, the propagator forces every
//    other orbit on that line to false, explaining each step with the
//    ternary clause (-a | -c | -b). This replaces millions of clauses (or
//    counter variables) with a few bitset-free line walks per assignment.
//  * Optionally the propagator makes the decisions itself: it picks the
//    row pair with the fewest remaining options and tries one of them.
//
// Supported classes: rot4 (even n) and rct4 (odd n). For rct4 the chosen
// set is invariant only under the half turn, so lines through each point
// of a chosen orbit are walked (not just through its representative).
//
// Build: see the Makefile. Usage: no3_cdcl <n> <rot4|rct4> [options]

#include <cadical.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <tuple>
#include <vector>

namespace {

using Cell = std::pair<int, int>; // (row, col)

int gcd(int a, int b) {
  a = std::abs(a);
  b = std::abs(b);
  while (b) {
    int t = a % b;
    a = b;
    b = t;
  }
  return a;
}

// ---------------------------------------------------------------------------
// The symmetry model: orbit variables and the cells they cover.

struct Model {
  int n = 0;
  bool rct4 = false;
  std::vector<std::vector<Cell>> orbits;  // cells of each variable
  std::vector<int> cell_orbit;            // n*n entries, -1 if excluded
  std::vector<std::vector<int>> rows;     // row r: orbit per point slot
  std::vector<int> diagonal_pairs;        // rct4: the exactly-one group

  Cell rot90(Cell p) const { return {p.second, n - 1 - p.first}; }
  Cell rot180(Cell p) const { return {n - 1 - p.first, n - 1 - p.second}; }

  explicit Model(int n_, bool rct4_) : n(n_), rct4(rct4_) {
    cell_orbit.assign(n * n, -1);
    Cell center{n / 2, n / 2};
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++) {
        Cell p{r, c};
        if (cell_orbit[r * n + c] >= 0) continue;
        std::vector<Cell> orbit;
        bool diag = r == c || r + c == n - 1;
        if (rct4 && diag) {
          if (r != c || (n % 2 == 1 && p == center)) continue;
          orbit = {p, rot180(p)};
        } else {
          Cell q = rot90(p);
          orbit = {p, q, rot90(q), rot90(rot90(q))};
        }
        for (Cell x : orbit) cell_orbit[x.first * n + x.second] = (int)orbits.size();
        orbits.push_back(orbit);
      }
    if (rct4)
      for (int a = 0; a < n / 2; a++) diagonal_pairs.push_back(cell_orbit[a * n + a]);
    rows.assign(n, {});
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++)
        if (cell_orbit[r * n + c] >= 0) rows[r].push_back(cell_orbit[r * n + c]);
  }
};

// ---------------------------------------------------------------------------
// Lines with at least three cells, indexed by (primitive direction, offset).

struct Lines {
  int n;
  std::vector<int> dir_base;              // per difference vector, -1 if short
  std::vector<std::pair<int, int>> dir_ab;
  std::vector<uint32_t> start;            // per line id
  std::vector<int32_t> orbit_ids;         // concatenated, deduplicated

  int line(Cell p, Cell q) const {
    int n1 = n - 1, w = 2 * n - 1;
    int dr = q.first - p.first, dc = q.second - p.second;
    int idx = (dr + n1) * w + (dc + n1);
    int base = dir_base[idx];
    if (base < 0) return -1;
    auto [a, b] = dir_ab[idx];
    int k = b * p.first - a * p.second + a * n1 + (b < 0 ? -b * n1 : 0);
    return base + k;
  }

  Lines(const Model &m) : n(m.n) {
    int n1 = n - 1, w = 2 * n - 1, max_step = n1 / 2;
    dir_base.assign(w * w, -1);
    dir_ab.assign(w * w, {0, 0});
    std::map<std::pair<int, int>, int> bases;
    int num = 0;
    for (int dr = -n1; dr <= n1; dr++)
      for (int dc = -n1; dc <= n1; dc++) {
        if (!dr && !dc) continue;
        int g = gcd(dr, dc), a = dr / g, b = dc / g;
        if (a < 0 || (a == 0 && b < 0)) a = -a, b = -b;
        if (std::max(a, std::abs(b)) > max_step) continue;
        auto it = bases.find({a, b});
        if (it == bases.end()) {
          it = bases.emplace(std::make_pair(a, b), num).first;
          num += (a + std::abs(b)) * n1 + 1;
        }
        int idx = (dr + n1) * w + (dc + n1);
        dir_base[idx] = it->second;
        dir_ab[idx] = {a, b};
      }
    std::vector<std::vector<int32_t>> lists(num);
    // Enumerate maximal lines of >= 3 cells and record their orbits.
    for (auto [ab, base] : bases) {
      auto [a, b] = ab;
      for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++) {
          int pr = r - a, pc = c - b;
          if (pr >= 0 && pr < n && pc >= 0 && pc < n) continue;
          std::vector<Cell> cells;
          for (int rr = r, cc = c; rr >= 0 && rr < n && cc >= 0 && cc < n; rr += a, cc += b)
            cells.push_back({rr, cc});
          if (cells.size() < 3) continue;
          std::vector<int32_t> ids;
          for (Cell x : cells) {
            int o = m.cell_orbit[x.first * n + x.second];
            if (o >= 0) ids.push_back(o);
          }
          std::sort(ids.begin(), ids.end());
          ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
          lists[line(cells[0], cells[1])] = ids;
        }
    }
    start.reserve(num + 1);
    for (auto &l : lists) {
      start.push_back((uint32_t)orbit_ids.size());
      orbit_ids.insert(orbit_ids.end(), l.begin(), l.end());
    }
    start.push_back((uint32_t)orbit_ids.size());
  }
};

// ---------------------------------------------------------------------------
// Independent final check: 2n points, no three collinear.

bool check_solution(int n, const std::vector<Cell> &pts, std::string &why) {
  if ((int)pts.size() != 2 * n) {
    why = "wrong number of points";
    return false;
  }
  std::map<std::tuple<long, long, long>, int> seen;
  for (size_t i = 0; i < pts.size(); i++)
    for (size_t j = i + 1; j < pts.size(); j++) {
      long dr = pts[j].first - pts[i].first, dc = pts[j].second - pts[i].second;
      if (!dr && !dc) {
        why = "duplicate point";
        return false;
      }
      long g = gcd((int)dr, (int)dc);
      dr /= g, dc /= g;
      if (dr < 0 || (dr == 0 && dc < 0)) dr = -dr, dc = -dc;
      auto key = std::make_tuple(dr, dc, dc * pts[i].first - dr * pts[i].second);
      if (seen.count(key)) {
        why = "three points on a line";
        return false;
      }
      seen[key] = 1;
    }
  return true;
}

const char *ALPHABET =
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz#$%&@?!()[]<>{}=*+|-/~^_:;,.";

std::string encode(char cls, std::vector<Cell> pts) {
  std::sort(pts.begin(), pts.end());
  std::string s(1, cls);
  for (Cell p : pts) s += ALPHABET[p.second];
  return s;
}

// ---------------------------------------------------------------------------
// The external propagator.

class Propagator : public CaDiCaL::ExternalPropagator {
public:
  Propagator(const Model &m, const Lines &l, bool decide, uint64_t seed)
      : M(m), L(l), decide_(decide), rng(seed) {
    int k = (int)M.orbits.size();
    value.assign(k, 0);
    reason.assign(k, {-1, -1});
    // Row pairs: for rot4/rct4 the constraints of row r and its mirror
    // coincide, so decisions only look at rows r < n/2 (+ the middle row).
    for (int r = 0; r < (M.n + 1) / 2; r++) decision_rows.push_back(r);
    level_start.push_back(0);
  }

  // Assignment bookkeeping --------------------------------------------------
  void notify_assignment(const std::vector<int> &lits) override {
    for (int lit : lits) {
      int o = std::abs(lit) - 1;
      if (o >= (int)value.size() || value[o]) continue;
      value[o] = lit > 0 ? 1 : -1;
      assigned.push_back(o);
      if (lit > 0) chosen.push_back(o);
    }
  }

  void notify_new_decision_level() override { level_start.push_back(assigned.size()); }

  // level_start[l] is the trail size when decision level l began.
  void notify_backtrack(size_t new_level) override {
    size_t keep = new_level + 1 < level_start.size() ? level_start[new_level + 1] : assigned.size();
    while (assigned.size() > keep) {
      int o = assigned.back();
      assigned.pop_back();
      if (value[o] > 0) chosen.pop_back();
      value[o] = 0;
    }
    level_start.resize(new_level + 1);
    processed = std::min(processed, chosen.size());
    pending.clear();
    pending_pos = 0;
  }

  // Propagation ---------------------------------------------------------------
  // For every newly chosen orbit a and every earlier chosen orbit c, walk the
  // lines through (a point of) a and each point of c. Any other orbit on such
  // a line must be false; if one is true we report the violated clause.
  void scan() {
    if (!conflict.empty()) return;
    while (processed < chosen.size()) {
      int a = chosen[processed];
      const auto &pa = M.orbits[a];
      // For rot4 the chosen set is quarter-turn invariant, so lines through
      // the representative cover all points of a. For rct4 walk all points.
      size_t na = M.rct4 ? pa.size() : 1;
      for (size_t idx = 0; idx < processed; idx++) {
        int c = chosen[idx];
        for (size_t i = 0; i < na; i++)
          for (Cell q : M.orbits[c]) {
            int line = L.line(pa[i], q);
            if (line < 0) continue;
            for (uint32_t t = L.start[line]; t < L.start[line + 1]; t++) {
              int b = L.orbit_ids[t];
              if (b == a || b == c) continue;
              if (value[b] > 0) {
                // Leave `processed` at a: after backtracking it is rescanned.
                conflict = {-(a + 1), -(c + 1), -(b + 1)};
                return;
              }
              if (value[b] == 0) pending.push_back({b, a, c});
            }
          }
      }
      processed++;
    }
  }

  int cb_propagate() override {
    scan();
    if (!conflict.empty()) return 0; // delivered as an external clause
    while (pending_pos < pending.size()) {
      auto &p = pending[pending_pos];
      if (value[p[0]] == 0) {
        reason[p[0]] = {p[1], p[2]};
        pending_pos++;
        propagations++;
        return -(p[0] + 1);
      }
      pending_pos++;
    }
    pending.clear();
    pending_pos = 0;
    return 0;
  }

  int cb_add_reason_clause_lit(int propagated_lit) override {
    if (reason_pos == 0) {
      auto [a, c] = reason[std::abs(propagated_lit) - 1];
      reason_buf = {propagated_lit, -(a + 1), -(c + 1)};
    }
    if (reason_pos < reason_buf.size()) return reason_buf[reason_pos++];
    reason_pos = 0;
    return 0;
  }

  bool cb_has_external_clause(bool &forgettable) override {
    scan();
    forgettable = true;
    return !conflict.empty();
  }

  int cb_add_external_clause_lit() override {
    if (conflict_pos < conflict.size()) return conflict[conflict_pos++];
    conflict.clear();
    conflict_pos = 0;
    return 0;
  }

  bool cb_check_found_model(const std::vector<int> &model) override {
    std::vector<Cell> pts;
    std::vector<int> owner;
    for (int lit : model)
      if (lit > 0 && lit <= (int)M.orbits.size())
        for (Cell p : M.orbits[lit - 1]) pts.push_back(p), owner.push_back(lit - 1);
    std::string why;
    if (check_solution(M.n, pts, why)) return true;
    // Every violation should have been caught during propagation; if one
    // slips through, block the offending orbits explicitly.
    std::fprintf(stderr, "c model rejected (%s); adding a blocking clause\n", why.c_str());
    for (size_t i = 0; i < pts.size() && conflict.empty(); i++)
      for (size_t j = i + 1; j < pts.size() && conflict.empty(); j++)
        for (size_t t = j + 1; t < pts.size(); t++) {
          long cross = (long)(pts[j].first - pts[i].first) * (pts[t].second - pts[i].second) -
                       (long)(pts[j].second - pts[i].second) * (pts[t].first - pts[i].first);
          if (cross == 0) {
            std::vector<int> os{owner[i], owner[j], owner[t]};
            std::sort(os.begin(), os.end());
            os.erase(std::unique(os.begin(), os.end()), os.end());
            for (int o : os) conflict.push_back(-(o + 1));
            break;
          }
        }
    return false;
  }

  // Decisions ----------------------------------------------------------------
  // Pick the row with the smallest slack (unassigned options minus missing
  // points) and choose one of its unassigned orbits at random.
  int cb_decide() override {
    if (!decide_) return 0;
    int best_row = -1, best_slack = 1 << 30, ties = 0;
    for (int r : decision_rows) {
      int have = 0, free = 0;
      for (int o : M.rows[r]) {
        if (value[o] > 0) have++;
        else if (value[o] == 0) free++;
      }
      int need = 2 - have;
      if (need <= 0 || free == 0) continue;
      int slack = free - need;
      if (slack < best_slack) {
        best_slack = slack, best_row = r, ties = 1;
      } else if (slack == best_slack && (int)(rng() % ++ties) == 0) {
        best_row = r;
      }
    }
    if (best_row < 0) return 0;
    std::vector<int> free;
    for (int o : M.rows[best_row])
      if (value[o] == 0) free.push_back(o);
    decisions++;
    return free[rng() % free.size()] + 1;
  }

  uint64_t propagations = 0, decisions = 0;

private:
  const Model &M;
  const Lines &L;
  bool decide_;
  std::mt19937_64 rng;
  std::vector<int> value;        // per orbit: 1 true, -1 false, 0 unassigned
  std::vector<int> assigned;     // assignment trail (observed variables)
  std::vector<int> chosen;       // true orbits in trail order
  std::vector<size_t> level_start;
  std::vector<int> decision_rows;
  size_t processed = 0;          // chosen[0..processed) have been scanned
  std::vector<std::array<int, 3>> pending;
  size_t pending_pos = 0;
  std::vector<std::pair<int, int>> reason; // per orbit: the pair that forced it false
  std::vector<int> reason_buf;
  size_t reason_pos = 0;
  std::vector<int> conflict;
  size_t conflict_pos = 0;
};

// ---------------------------------------------------------------------------
// Static clauses.

struct Encoder {
  CaDiCaL::Solver &S;
  int next_var;

  int fresh() { return ++next_var; }

  // sum(weights * lits) == 2, weights in {1, 2}; counter for <= 2.
  void exactly_two(const std::vector<std::pair<int, int>> &terms) {
    // at least two: some term true, and a weight-1 term is not alone.
    std::vector<int> all;
    for (auto [l, w] : terms) all.push_back(l);
    S.clause(all);
    for (auto [l, w] : terms) {
      if (w != 1) continue;
      std::vector<int> cl{-l};
      for (auto [l2, w2] : terms)
        if (l2 != l) cl.push_back(l2);
      S.clause(cl);
    }
    // at most two: sequential counter.
    int p1 = 0, p2 = 0;
    for (size_t i = 0; i < terms.size(); i++) {
      auto [x, w] = terms[i];
      if (i > 0) S.clause(-x, w == 2 ? -p1 : -p2);
      if (i + 1 == terms.size()) break;
      int s1 = fresh(), s2 = fresh();
      S.clause(-x, s1);
      if (w == 2) S.clause(-x, s2);
      if (i > 0) {
        S.clause(-p1, s1);
        S.clause(-p2, s2);
        if (w == 1) S.clause(-p1, -x, s2);
      }
      p1 = s1, p2 = s2;
    }
  }
};

void usage() {
  std::fprintf(stderr,
               "usage: no3_cdcl <n> <rot4|rct4> [--seed K] [--own-decisions] [--no-mirror]\n"
               "                [--count C] [--time-limit SECONDS]\n");
  std::exit(2);
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 3) usage();
  int n = std::atoi(argv[1]);
  std::string cls = argv[2];
  uint64_t seed = 1;
  bool own = false, mirror = true;
  int count = 1;
  double time_limit = 0;
  for (int i = 3; i < argc; i++) {
    std::string a = argv[i];
    if (a == "--seed" && i + 1 < argc) seed = std::strtoull(argv[++i], nullptr, 10);
    else if (a == "--own-decisions") own = true;
    else if (a == "--no-mirror") mirror = false;
    else if (a == "--count" && i + 1 < argc) count = std::atoi(argv[++i]);
    else if (a == "--time-limit" && i + 1 < argc) time_limit = std::atof(argv[++i]);
    else usage();
  }
  bool rct4 = cls == "rct4";
  if ((cls != "rot4" && !rct4) || n < 4 || n > 90 || (rct4 ? n % 2 == 0 : n % 2 == 1)) usage();

  auto t0 = std::chrono::steady_clock::now();
  Model M(n, rct4);
  Lines L(M);
  int k = (int)M.orbits.size();

  CaDiCaL::Solver S;
  S.set("seed", (int)(seed & 0x7fffffff));
  Propagator P(M, L, own, seed);
  S.connect_external_propagator(&P);
  for (int v = 1; v <= k; v++) S.add_observed_var(v);
  Encoder E{S, k};

  // Rows and columns: exactly two points. Deduplicate identical constraints.
  std::map<std::vector<std::pair<int, int>>, int> seen;
  for (int i = 0; i < n; i++)
    for (int col = 0; col < 2; col++) {
      std::map<int, int> w;
      for (int j = 0; j < n; j++) {
        int r = col ? j : i, c = col ? i : j;
        int o = M.cell_orbit[r * n + c];
        if (o >= 0) w[o + 1]++;
      }
      std::vector<std::pair<int, int>> terms(w.begin(), w.end());
      if (seen.emplace(terms, 1).second) E.exactly_two(terms);
    }
  // Orbits on a line through two points of another orbit conflict with it.
  // (For rct4 also pairs whose own points are collinear with a third cell of
  // the same pair are excluded; a pair has only two cells so none exist.)
  for (int o = 0; o < k; o++) {
    const auto &pts = M.orbits[o];
    for (size_t i = 0; i < pts.size(); i++)
      for (size_t j = i + 1; j < pts.size(); j++) {
        int line = L.line(pts[i], pts[j]);
        if (line < 0) continue;
        for (uint32_t t = L.start[line]; t < L.start[line + 1]; t++) {
          int q = L.orbit_ids[t];
          if (q != o && q > o) S.clause(-(o + 1), -(q + 1));
          else if (q != o && q < o) S.clause(-(q + 1), -(o + 1));
        }
      }
  }
  if (rct4) {
    std::vector<int> group;
    for (int o : M.diagonal_pairs) group.push_back(o + 1);
    S.clause(group);
    for (size_t i = 0; i < group.size(); i++)
      for (size_t j = i + 1; j < group.size(); j++) S.clause(-group[i], -group[j]);
  }
  if (mirror) {
    // See src/cnf.rs (break_mirror): top-row columns c1 < c2 need c1 + c2 <= n - 1.
    for (int c1 = 0; c1 < n; c1++)
      for (int c2 = c1 + 1; c2 < n; c2++) {
        if (c1 + c2 <= n - 1) continue;
        int a = M.cell_orbit[c1], b = M.cell_orbit[c2];
        if (a >= 0 && b >= 0 && a != b) S.clause(-(a + 1), -(b + 1));
      }
  }
  double setup = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::fprintf(stderr, "c n=%d class=%s orbits=%d setup=%.2fs\n", n, cls.c_str(), k, setup);

  if (time_limit > 0) {
    // CaDiCaL has no wall-clock limit; approximate one with a terminator.
    struct Timer : CaDiCaL::Terminator {
      std::chrono::steady_clock::time_point end;
      bool terminate() override { return std::chrono::steady_clock::now() > end; }
    };
    static Timer timer;
    timer.end = std::chrono::steady_clock::now() + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(time_limit));
    S.connect_terminator(&timer);
  }

  int found = 0;
  while (found < count) {
    int res = S.solve();
    double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (res != 10) {
      std::fprintf(stderr, "c %s after %.1fs, %d solutions\n", res == 20 ? "UNSAT (no more solutions)" : "UNKNOWN", t, found);
      return res == 20 ? 20 : 0;
    }
    std::vector<Cell> pts;
    std::vector<int> block;
    for (int v = 1; v <= k; v++)
      if (S.val(v) > 0) {
        for (Cell p : M.orbits[v - 1]) pts.push_back(p);
        block.push_back(-v);
      }
    std::string why;
    if (!check_solution(n, pts, why)) {
      std::fprintf(stderr, "c INTERNAL ERROR: invalid solution (%s)\n", why.c_str());
      return 1;
    }
    std::printf("%s\n", encode(rct4 ? 'c' : 'o', pts).c_str());
    std::fflush(stdout);
    std::fprintf(stderr, "c solution %d after %.1fs (conflicts %lld, decisions %lld, propagator decisions %llu, external propagations %llu)\n",
                 found + 1, t, (long long)S.get_statistic_value("conflicts"), (long long)S.get_statistic_value("decisions"),
                 (unsigned long long)P.decisions, (unsigned long long)P.propagations);
    found++;
    S.clause(block);
  }
  return 10;
}
