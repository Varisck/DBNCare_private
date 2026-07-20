// Rcpp implementation of the sampling logic of R/dbn_sampling.R.
//
// The R wrapper (dbn.sampling.cpp, in R/dbn_sampling_cpp.R) compiles the
// dbn.fit object down to integer-indexed "plans"; the cores below run the
// sampling loops on a flat value matrix with one row per (observation, time)
// pair and one column per variable: row = obs * (max_time + 1) + t. This is
// the same layout used by the time-series dictionary of dbn.sampling(),
// where the value of a variable at time t for observation obs sits at
// (obs - 1) * (max_time + 1) + t + 1.
//
// Both cores consume the R RNG stream exactly like the pure-R implementation:
//   * gaussian: one rnorm() per prior node, then rnorm(max_time) per
//     transition node, per observation (the values themselves are
//     deterministic given the noise);
//   * discrete: one draw per sampled value, replicating R's
//     sample(x, size = 1, prob = p) internals (FixupProb +
//     ProbSampleNoReplace, with revsort() from the R API).
// For a fixed seed, dbn.sampling() and dbn.sampling.cpp() therefore return
// the same dataset.

#include <R_ext/Random.h>  // unif_rand
#include <R_ext/Utils.h>   // revsort
#include <Rcpp.h>

using namespace Rcpp;

namespace {

// R matrices cannot have more than 2^31 - 1 rows
void check_matrix_size(double n_rows) {
  if (n_rows > 2147483647.0)
    stop(
        "dbn.sampling.cpp: n_samples * (max_time + 1) exceeds the maximum "
        "matrix size");
}

// --------------------------------------------------------------------------
// gaussian: mirrors sample_variable_gaussian()
// --------------------------------------------------------------------------

struct GaussNode {
  int var;   // column written by this node
  int time;  // slice this node writes to (G_0 nodes: 0..g0_span-1; transition
             // nodes: 0, unused - the transition loop supplies the time)
  double intercept;
  double std;                    // residual standard deviation
  std::vector<int> par_var;      // parent columns
  std::vector<int> par_lag;      // parent lags (0 = same slice, k = t-k)
  std::vector<double> par_coef;  // regression coefficients, parents order
};

std::vector<GaussNode> parse_gauss_plan(const List& plan, int n_vars) {
  std::vector<GaussNode> nodes;
  nodes.reserve(plan.size());
  for (R_xlen_t k = 0; k < plan.size(); ++k) {
    List spec = plan[k];
    GaussNode node;
    node.var = as<int>(spec["var"]);
    node.time = as<int>(spec["time"]);
    node.intercept = as<double>(spec["intercept"]);
    node.std = as<double>(spec["std"]);
    node.par_var = as<std::vector<int>>(spec["par_var"]);
    node.par_lag = as<std::vector<int>>(spec["par_lag"]);
    node.par_coef = as<std::vector<double>>(spec["par_coef"]);
    if (node.var < 0 || node.var >= n_vars || node.time < 0 ||
        node.par_var.size() != node.par_lag.size() ||
        node.par_var.size() != node.par_coef.size())
      stop("dbn.sampling.cpp: malformed gaussian sampling plan");
    for (size_t j = 0; j < node.par_var.size(); ++j)
      if (node.par_var[j] < 0 || node.par_var[j] >= n_vars ||
          node.par_lag[j] < 0)
        stop("dbn.sampling.cpp: malformed gaussian sampling plan");
    nodes.push_back(std::move(node));
  }
  return nodes;
}

// intercept + sum coef_j * value(parent_j, t - lag_j). Parents not available
// yet (t - lag < 0, the time_point <= 0 case of get_time_point()) are
// skipped, like in sample_variable_gaussian(). Products are accumulated in
// long double to mirror R's sum().
double gauss_mean(const GaussNode& node, const NumericMatrix& values,
                  R_xlen_t block, int t) {
  long double acc = node.intercept;  // first product: 1 * intercept
  for (size_t j = 0; j < node.par_var.size(); ++j) {
    const int tp = t - node.par_lag[j];
    if (tp >= 0)
      acc += static_cast<double>(values(block + tp, node.par_var[j]) *
                                 node.par_coef[j]);
  }
  return static_cast<double>(acc);
}

// --------------------------------------------------------------------------
// discrete: mirrors sample_variable_discrete() + filter_cpt()
// --------------------------------------------------------------------------

// A CPT is kept as its flat probability vector (column-major, the linear
// order of the prob array, which is also the row order of
// data.frame(as.table(prob))). dims[0] indexes the node's own levels, the
// following dims its parents. Level labels are translated once to canonical
// integer codes per variable, so that the string comparisons of filter_cpt()
// become integer comparisons.
struct DiscNode {
  int var;                    // column written by this node
  int time;                   // slice this node writes to (see GaussNode::time)
  std::vector<int> dims;      // dims[0] = own level count
  std::vector<int> own_code;  // CPT dim-1 index -> canonical code
  std::vector<double> freq;   // CPT probabilities, column-major
  std::vector<int> par_var;   // per parent dim, in CPT dim order
  std::vector<int> par_lag;
  std::vector<std::vector<int>> par_code;  // CPT dim index -> canonical code
  std::vector<std::vector<int>> par_inv;   // canonical code -> CPT dim index
  R_xlen_t n_rows;                         // prod(dims)
};

std::vector<DiscNode> parse_disc_plan(const List& plan, int n_vars,
                                      const IntegerVector& n_levels) {
  std::vector<DiscNode> nodes;
  nodes.reserve(plan.size());
  for (R_xlen_t k = 0; k < plan.size(); ++k) {
    List spec = plan[k];
    DiscNode node;
    node.var = as<int>(spec["var"]);
    node.time = as<int>(spec["time"]);
    node.dims = as<std::vector<int>>(spec["dims"]);
    node.own_code = as<std::vector<int>>(spec["own_map"]);
    node.freq = as<std::vector<double>>(spec["freq"]);
    node.par_var = as<std::vector<int>>(spec["par_var"]);
    node.par_lag = as<std::vector<int>>(spec["par_lag"]);
    List par_map = spec["par_map"];

    if (node.var < 0 || node.var >= n_vars || node.time < 0 ||
        node.dims.empty() || node.dims.size() != node.par_var.size() + 1 ||
        node.par_var.size() != node.par_lag.size() ||
        static_cast<size_t>(par_map.size()) != node.par_var.size())
      stop("dbn.sampling.cpp: malformed discrete sampling plan");

    node.n_rows = 1;
    for (int d : node.dims) {
      if (d <= 0) stop("dbn.sampling.cpp: malformed discrete sampling plan");
      node.n_rows *= d;
    }
    if (static_cast<R_xlen_t>(node.freq.size()) != node.n_rows ||
        static_cast<int>(node.own_code.size()) != node.dims[0])
      stop("dbn.sampling.cpp: malformed discrete sampling plan");

    node.par_code.resize(node.par_var.size());
    node.par_inv.resize(node.par_var.size());
    for (size_t j = 0; j < node.par_var.size(); ++j) {
      if (node.par_var[j] < 0 || node.par_var[j] >= n_vars ||
          node.par_lag[j] < 0)
        stop("dbn.sampling.cpp: malformed discrete sampling plan");

      node.par_code[j] = as<std::vector<int>>(par_map[j]);

      if (static_cast<int>(node.par_code[j].size()) != node.dims[j + 1])
        stop("dbn.sampling.cpp: malformed discrete sampling plan");

      // invert label translation: canonical code -> index in the CPT dim
      node.par_inv[j].assign(n_levels[node.par_var[j]], -1);
      for (size_t i = 0; i < node.par_code[j].size(); ++i) {
        const int code = node.par_code[j][i];
        if (code >= 0 && code < static_cast<int>(node.par_inv[j].size()))
          node.par_inv[j][code] = static_cast<int>(i);
      }
    }
    nodes.push_back(std::move(node));
  }
  return nodes;
}

// Replica of R's sample(x, size = 1, prob = p): FixupProb followed by
// ProbSampleNoReplace (src/main/random.c), using revsort() from the R API so
// that ties between equal probabilities are broken exactly the same way.
// `payload` carries the sampled values alongside the probabilities (R keeps
// 1..n there and indexes x afterwards; the permutation behaves identically).
int sample_one(std::vector<double>& prob, std::vector<int>& payload) {
  const int n = static_cast<int>(prob.size());
  double sum = 0.0;
  int npos = 0;
  for (int i = 0; i < n; ++i) {
    if (!R_FINITE(prob[i]) || prob[i] < 0.0)
      stop("dbn.sampling.cpp: NA in probability vector");
    if (prob[i] > 0.0) ++npos;
    sum += prob[i];
  }
  if (npos == 0) stop("dbn.sampling.cpp: too few positive probabilities");
  for (int i = 0; i < n; ++i) prob[i] /= sum;

  revsort(prob.data(), payload.data(), n);
  const double rT = unif_rand();  // totalmass = 1, one draw per sampled value
  double mass = 0.0;
  int j = 0;
  for (; j < n - 1; ++j) {
    mass += prob[j];
    if (rT <= mass) break;
  }
  return payload[j];
}

// Sample one value of `node` at time t: collect the CPT rows compatible with
// the values of the available parents (filter_cpt()), then draw the node's
// level with probability proportional to Freq (sample_variable_discrete()).
template <typename Matrix>
int sample_disc_value(const DiscNode& node, const Matrix& codes, R_xlen_t block,
                      int t, std::vector<double>& prob_buf,
                      std::vector<int>& code_buf) {
  const int npar = static_cast<int>(node.par_var.size());
  const int n_own = node.dims[0];
  prob_buf.clear();
  code_buf.clear();

  // a parent at lag k is observable iff t - k >= 0 (get_time_point() > 0);
  // unavailable parents are simply not filtered on, like in the R code
  bool all_observed = true;
  for (int j = 0; j < npar; ++j)
    if (t - node.par_lag[j] < 0) {
      all_observed = false;
      break;
    }

  if (all_observed) {
    // every parent fixes its dimension: the candidate rows are the n_own
    // consecutive entries at offset; same order as the filtered data.frame
    R_xlen_t offset = 0;
    R_xlen_t mult = 1;
    bool match = true;
    for (int j = 0; j < npar; ++j) {
      const int code = codes(block + (t - node.par_lag[j]), node.par_var[j]);
      const int idx =
          (code >= 0 && code < static_cast<int>(node.par_inv[j].size()))
              ? node.par_inv[j][code]
              : -1;
      if (idx < 0) {
        match = false;
        break;
      }  // value not among the CPT levels
      offset += mult * idx;
      mult *= node.dims[j + 1];
    }
    if (match) {
      for (int i = 0; i < n_own; ++i) {
        prob_buf.push_back(
            node.freq[i + static_cast<R_xlen_t>(n_own) * offset]);
        code_buf.push_back(node.own_code[i]);
      }
    }
  } else {
    // markov order > 1 at early time points: scan the CPT rows in order and
    // keep those whose observed parent dims match the sampled values
    for (R_xlen_t r = 0; r < node.n_rows; ++r) {
      const int own = static_cast<int>(r % n_own);
      R_xlen_t rest = r / n_own;
      bool keep = true;
      for (int j = 0; j < npar; ++j) {
        const int idx = static_cast<int>(rest % node.dims[j + 1]);
        rest /= node.dims[j + 1];
        if (t - node.par_lag[j] >= 0) {
          const int code =
              codes(block + (t - node.par_lag[j]), node.par_var[j]);
          if (node.par_code[j][idx] != code) {
            keep = false;
            break;
          }
        }
      }
      if (keep) {
        prob_buf.push_back(node.freq[r]);
        code_buf.push_back(node.own_code[own]);
      }
    }
  }

  if (prob_buf.empty())
    stop(
        "dbn.sampling.cpp: no CPT rows match the sampled parents values "
        "(inconsistent levels between the CPTs?)");
  return sample_one(prob_buf, code_buf);
}

}  // namespace

// Sampling core for gaussian DBNs. Returns the value matrix with one row per
// (observation, time) pair and one column per variable.
// [[Rcpp::export]]
NumericMatrix dbn_sample_gaussian_cpp(int n_samples, int max_time, int n_vars,
                                      List plan_0, List plan_t) {
  const std::vector<GaussNode> nodes_0 = parse_gauss_plan(plan_0, n_vars);
  const std::vector<GaussNode> nodes_t = parse_gauss_plan(plan_t, n_vars);
  const int block_len = max_time + 1;
  check_matrix_size(static_cast<double>(n_samples) * block_len);
  NumericMatrix values(n_samples * block_len, n_vars);

  // Number of initial slices covered by the prior network G_0. For a plain
  // (non-extended) G_0 every prior node sits at slice 0, so g0_span = 1 and the
  // transition network runs from t = 1. With an extended G_0 (var_0, var_1,
  // ...) the prior nodes span out over slices 0..g0_span-1  and the transition
  // takes over from t = g0_span.
  int g0_span = 1;
  for (const GaussNode& node : nodes_0)
    if (node.time + 1 > g0_span) g0_span = node.time + 1;

  // per-trajectory noise of the transition nodes, drawn before the time loop
  // exactly like dbn.sampling() does (rnorm(max_time, 0, std) per node)
  std::vector<double> noise(nodes_t.size() * static_cast<size_t>(max_time));

  for (int obs = 0; obs < n_samples; ++obs) {
    const R_xlen_t block = static_cast<R_xlen_t>(obs) * block_len;

    // prior network: one rnorm() per node in topological order, each written to
    // its own slice (block + node.time; slice 0 for a plain G_0)
    for (const GaussNode& node : nodes_0) {
      const double noise_0 = R::rnorm(0.0, node.std);
      if (node.time <= max_time)
        values(block + node.time, node.var) =
            gauss_mean(node, values, block, node.time) + noise_0;
    }

    for (size_t k = 0; k < nodes_t.size(); ++k)
      for (int t = 0; t < max_time; ++t)
        noise[k * max_time + t] = R::rnorm(0.0, nodes_t[k].std);

    // t = g0_span..max_time: transition network in topological order
    for (int t = g0_span; t <= max_time; ++t)
      for (size_t k = 0; k < nodes_t.size(); ++k)
        values(block + t, nodes_t[k].var) =
            gauss_mean(nodes_t[k], values, block, t) +
            noise[k * max_time + (t - 1)];
  }
  return values;
}

// Sampling core for discrete DBNs. Returns the matrix of 0-based canonical
// level codes; the R wrapper maps them back to the level labels.
// [[Rcpp::export]]
IntegerMatrix dbn_sample_discrete_cpp(int n_samples, int max_time, int n_vars,
                                      IntegerVector n_levels, List plan_0,
                                      List plan_t) {
  if (n_levels.size() != n_vars)
    stop("dbn.sampling.cpp: malformed discrete sampling plan");
  const std::vector<DiscNode> nodes_0 =
      parse_disc_plan(plan_0, n_vars, n_levels);
  const std::vector<DiscNode> nodes_t =
      parse_disc_plan(plan_t, n_vars, n_levels);
  const int block_len = max_time + 1;
  check_matrix_size(static_cast<double>(n_samples) * block_len);
  IntegerMatrix codes(n_samples * block_len, n_vars);

  // slices covered by G_0 (1 for a plain G_0; see dbn_sample_gaussian_cpp)
  int g0_span = 1;
  for (const DiscNode& node : nodes_0)
    if (node.time + 1 > g0_span) g0_span = node.time + 1;

  std::vector<double> prob_buf;
  std::vector<int> code_buf;

  for (int obs = 0; obs < n_samples; ++obs) {
    const R_xlen_t block = static_cast<R_xlen_t>(obs) * block_len;
    for (const DiscNode& node : nodes_0)
      if (node.time <= max_time)
        codes(block + node.time, node.var) = sample_disc_value(
            node, codes, block, node.time, prob_buf, code_buf);
    for (int t = g0_span; t <= max_time; ++t)
      for (const DiscNode& node : nodes_t)
        codes(block + t, node.var) =
            sample_disc_value(node, codes, block, t, prob_buf, code_buf);
  }
  return codes;
}

struct MixNode {
  int var;   // column written by this node (continuous)
  int time;  // slice this node writes to (see GaussNode::time)
  std::vector<double> intercept;  // one per discrete-parent combination
  std::vector<double> std;        // residual sd, one per combination
  std::vector<int> cont_par_var;  // continuous parent columns
  std::vector<int> cont_par_lag;  // continuous parent lags (0 = same slice)
  std::vector<std::vector<double>> par_coef;  // [combination][cont parent]

  std::vector<int> disc_dims;     // level count of each discrete parent
  std::vector<int> disc_par_var;  // discrete parent columns, coef-dim order
  std::vector<int> disc_par_lag;
  std::vector<std::vector<int>> par_code;  // coef dim index -> canonical code
  std::vector<std::vector<int>> par_inv;   // canonical code -> coef dim index
  R_xlen_t n_combos;                       // prod(disc_dims), 1 if no disc par
};

std::vector<MixNode> parse_mixed_plan(const List& plan, int n_vars,
                                      const IntegerVector& n_levels) {
  std::vector<MixNode> nodes;
  nodes.reserve(plan.size());

  for (R_xlen_t k = 0; k < plan.size(); ++k) {
    List spec = plan[k];
    MixNode node;

    node.var = as<int>(spec["var"]);
    node.time = as<int>(spec["time"]);
    node.intercept = as<std::vector<double>>(spec["intercept"]);
    node.std = as<std::vector<double>>(spec["std"]);

    node.cont_par_var = as<std::vector<int>>(spec["cont_par_var"]);
    node.cont_par_lag = as<std::vector<int>>(spec["cont_par_lag"]);

    node.disc_dims = as<std::vector<int>>(spec["disc_dims"]);
    node.disc_par_var = as<std::vector<int>>(spec["disc_par_var"]);
    node.disc_par_lag = as<std::vector<int>>(spec["disc_par_lag"]);

    List par_coef = spec["par_coef"];
    List disc_par_map = spec["disc_par_map"];

    if (node.var < 0 || node.var >= n_vars || node.time < 0 ||
        node.cont_par_var.size() != node.cont_par_lag.size() ||
        node.disc_par_var.size() != node.disc_par_lag.size() ||
        node.disc_par_var.size() != node.disc_dims.size() ||
        static_cast<size_t>(disc_par_map.size()) != node.disc_par_var.size())
      stop("dbn.sampling.cpp: malformed mixed sampling plan");

    // number of discrete-parent combinations: the size of the coef tables
    node.n_combos = 1;
    for (int d : node.disc_dims) {
      if (d <= 0) stop("dbn.sampling.cpp: malformed mixed sampling plan");
      node.n_combos *= d;
    }
    if (static_cast<R_xlen_t>(node.intercept.size()) != node.n_combos ||
        static_cast<R_xlen_t>(node.std.size()) != node.n_combos ||
        static_cast<R_xlen_t>(par_coef.size()) != node.n_combos)
      stop("dbn.sampling.cpp: malformed mixed sampling plan");

    // continuous parents: bounds + one coefficient vector per combination
    for (size_t j = 0; j < node.cont_par_var.size(); ++j)
      if (node.cont_par_var[j] < 0 || node.cont_par_var[j] >= n_vars ||
          node.cont_par_lag[j] < 0)
        stop("dbn.sampling.cpp: malformed mixed sampling plan");

    node.par_coef.resize(node.n_combos);
    for (R_xlen_t c = 0; c < node.n_combos; ++c) {
      node.par_coef[c] = as<std::vector<double>>(par_coef[c]);
      if (node.par_coef[c].size() != node.cont_par_var.size())
        stop("dbn.sampling.cpp: malformed mixed sampling plan");
    }

    // discrete parents: invert the label translation, same as parse_disc_plan
    node.par_code.resize(node.disc_par_var.size());
    node.par_inv.resize(node.disc_par_var.size());
    for (size_t j = 0; j < node.disc_par_var.size(); ++j) {
      if (node.disc_par_var[j] < 0 || node.disc_par_var[j] >= n_vars ||
          node.disc_par_lag[j] < 0)
        stop("dbn.sampling.cpp: malformed mixed sampling plan");

      node.par_code[j] = as<std::vector<int>>(disc_par_map[j]);
      if (static_cast<int>(node.par_code[j].size()) != node.disc_dims[j])
        stop("dbn.sampling.cpp: malformed mixed sampling plan");

      node.par_inv[j].assign(n_levels[node.disc_par_var[j]], -1);
      for (size_t i = 0; i < node.par_code[j].size(); ++i) {
        const int code = node.par_code[j][i];
        if (code >= 0 && code < static_cast<int>(node.par_inv[j].size()))
          node.par_inv[j][code] = static_cast<int>(i);
      }
    }

    nodes.push_back(std::move(node));
  }
  return nodes;
}

// Sample one value of a CLG node at time t from the single value matrix
double sample_mixed_value(const MixNode& node, const NumericMatrix& values,
                          R_xlen_t block, int t) {
  R_xlen_t combo = 0;
  R_xlen_t mult = 1;
  for (size_t j = 0; j < node.disc_par_var.size(); ++j) {
    const int tp = t - node.disc_par_lag[j];
    int idx;
    if (tp < 0) {
      // if discrete parent not in time series < markov order
      // sample a discrete parent of the cg at random
      idx = static_cast<int>(unif_rand() * node.disc_dims[j]);
      if (idx >= node.disc_dims[j]) idx = node.disc_dims[j] - 1;
    } else {
      const int code = values(block + tp, node.disc_par_var[j]);
      idx = (code >= 0 && code < static_cast<int>(node.par_inv[j].size()))
                ? node.par_inv[j][code]
                : -1;
    }
    if (idx < 0)
      stop("dbn.sampling.cpp: mixed node parent value not among CLG levels");
    combo += mult * idx;
    mult *= node.disc_dims[j];
  }

  long double acc = node.intercept[combo];  // first product: 1 * intercept
  const std::vector<double>& coef = node.par_coef[combo];
  for (size_t j = 0; j < node.cont_par_var.size(); ++j) {
    const int tp = t - node.cont_par_lag[j];
    // if tp < 0 it does not exist the value in the time series ignore regressor
    if (tp >= 0)
      acc += static_cast<double>(values(block + tp, node.cont_par_var[j]) *
                                 coef[j]);
  }
  return static_cast<double>(acc) + R::rnorm(0.0, node.std[combo]);
}

// Sampling core for mixed (conditional-gaussian) DBNs
// [[Rcpp::export]]
NumericMatrix dbn_sample_mixed_cpp(int n_samples, int max_time, List n_vars,
                                   IntegerVector n_levels, List plan_0,
                                   List plan_t, List node_ordering_0,
                                   List node_ordering_t) {
  const int n_total = as<int>(n_vars["discrete"]) +
                      as<int>(n_vars["gaussian"]) + as<int>(n_vars["mixed"]);
  if (n_levels.size() != n_total)
    stop("dbn.sampling.cpp: malformed mixed sampling plan");

  const std::vector<DiscNode> disc_nodes_0 =
      parse_disc_plan(plan_0["discrete"], n_total, n_levels);
  const std::vector<DiscNode> disc_nodes_t =
      parse_disc_plan(plan_t["discrete"], n_total, n_levels);
  const std::vector<GaussNode> gauss_nodes_0 =
      parse_gauss_plan(plan_0["gaussian"], n_total);
  const std::vector<GaussNode> gauss_nodes_t =
      parse_gauss_plan(plan_t["gaussian"], n_total);
  const std::vector<MixNode> mix_nodes_0 =
      parse_mixed_plan(plan_0["mixed"], n_total, n_levels);
  const std::vector<MixNode> mix_nodes_t =
      parse_mixed_plan(plan_t["mixed"], n_total, n_levels);

  const IntegerVector order_type_0 = node_ordering_0["vector_type"];
  const IntegerVector order_idx_0 = node_ordering_0["index_vec"];
  const IntegerVector order_type_t = node_ordering_t["vector_type"];
  const IntegerVector order_idx_t = node_ordering_t["index_vec"];
  if (order_type_0.size() != order_idx_0.size() ||
      order_type_t.size() != order_idx_t.size())
    stop("dbn.sampling.cpp: malformed mixed node ordering");

  const int block_len = max_time + 1;
  check_matrix_size(static_cast<double>(n_samples) * block_len);
  NumericMatrix values(n_samples * block_len, n_total);

  // slices covered by G_0 (1 for a plain G_0; see dbn_sample_gaussian_cpp)
  int g0_span = 1;
  for (const DiscNode& n : disc_nodes_0)
    if (n.time + 1 > g0_span) g0_span = n.time + 1;
  for (const GaussNode& n : gauss_nodes_0)
    if (n.time + 1 > g0_span) g0_span = n.time + 1;
  for (const MixNode& n : mix_nodes_0)
    if (n.time + 1 > g0_span) g0_span = n.time + 1;

  std::vector<double> prob_buf;
  std::vector<int> code_buf;

  // sample every node of one time slice in topological order. When g0 is true
  // each node is written to its own slice (node.time) instead of the loop's t,
  // so an extended G_0 fans its prior nodes over slices 0..g0_span-1.
  auto run_slice = [&](const IntegerVector& type, const IntegerVector& idx,
                       R_xlen_t block, int t, bool g0,
                       const std::vector<DiscNode>& dn,
                       const std::vector<GaussNode>& gn,
                       const std::vector<MixNode>& mn) {
    for (int k = 0; k < type.size(); ++k) {
      const int i = idx[k];
      switch (type[k]) {
        case 0: {  // discrete
          if (i < 0 || static_cast<size_t>(i) >= dn.size())
            stop("dbn.sampling.cpp: mixed node ordering index out of range");
          const int tt = g0 ? dn[i].time : t;
          if (tt <= max_time)
            values(block + tt, dn[i].var) =
                sample_disc_value(dn[i], values, block, tt, prob_buf, code_buf);
          break;
        }
        case 1: {  // gaussian
          if (i < 0 || static_cast<size_t>(i) >= gn.size())
            stop("dbn.sampling.cpp: mixed node ordering index out of range");
          const int tt = g0 ? gn[i].time : t;
          if (tt <= max_time)
            values(block + tt, gn[i].var) =
                gauss_mean(gn[i], values, block, tt) + R::rnorm(0.0, gn[i].std);
          break;
        }
        case 2: {  // mixed / CLG
          if (i < 0 || static_cast<size_t>(i) >= mn.size())
            stop("dbn.sampling.cpp: mixed node ordering index out of range");
          const int tt = g0 ? mn[i].time : t;
          if (tt <= max_time)
            values(block + tt, mn[i].var) =
                sample_mixed_value(mn[i], values, block, tt);
          break;
        }
        default:
          stop("dbn.sampling.cpp: unknown node type in mixed ordering");
      }
    }
  };

  for (int obs = 0; obs < n_samples; ++obs) {
    const R_xlen_t block = static_cast<R_xlen_t>(obs) * block_len;
    run_slice(order_type_0, order_idx_0, block, 0, true, disc_nodes_0,
              gauss_nodes_0, mix_nodes_0);
    for (int t = g0_span; t <= max_time; ++t)
      run_slice(order_type_t, order_idx_t, block, t, false, disc_nodes_t,
                gauss_nodes_t, mix_nodes_t);
  }
  return values;
}
