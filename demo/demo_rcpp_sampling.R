# ============================================================================
# Verification of dbn.sampling.cpp (Rcpp port of dbn.sampling)
#
# For a hand-designed gaussian DBN and a hand-designed discrete DBN:
#   1. dbn.sampling.cpp() must return the same dataset as dbn.sampling()
#      when run with the same seed (the C++ cores consume the R RNG stream
#      draw for draw);
#   2. hc.dbn() run on the sampled dataset must reconstruct the designed
#      structure (modelstring comparison against the ground truth).
# A timing comparison between the two implementations closes the script.
# ============================================================================

if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(DynamicBayesianNetwork)
}

section <- function(title) cat("\n==", title, "==\n")
check <- function(label, ok) {
  cat(sprintf("[%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label))
  if (!isTRUE(ok)) stop("check failed: ", label, call. = FALSE)
}

# structural equality of two modelstrings: same nodes and same arcs, no matter
# the order in which bnlearn happens to list the parents inside the brackets
same_structure <- function(ms_a, ms_b) {
  isTRUE(all.equal(bnlearn::model2network(ms_a), bnlearn::model2network(ms_b)))
}

# ----------------------------------------------------------------------------
# Hand-designed structure, shared by the gaussian and the discrete model
#
#   G_0:  A_0 -> C_0 <- B_0                 (collider: orientations identifiable)
#   G_t:  A_t-1 -> A_t, B_t-1 -> B_t, C_t-1 -> C_t   (persistence)
#         A_t -> C_t <- B_t                 (collider: orientations identifiable)
# ----------------------------------------------------------------------------
build_structure <- function() {
  dbn <- empty.dbn(dynamic_nodes = c("A", "B", "C"), markov_order = 1)
  dbn <- add.arc.dbn(dbn, from = c("A", "t_0"), to = c("C", "t_0"))
  dbn <- add.arc.dbn(dbn, from = c("B", "t_0"), to = c("C", "t_0"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t-1"), to = c("A", "t"))
  dbn <- add.arc.dbn(dbn, from = c("B", "t-1"), to = c("B", "t"))
  dbn <- add.arc.dbn(dbn, from = c("C", "t-1"), to = c("C", "t"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t"),   to = c("C", "t"))
  dbn <- add.arc.dbn(dbn, from = c("B", "t"),   to = c("C", "t"))
  dbn
}

truth <- build_structure()
truth_ms <- modelstring(truth)
cat("designed G_0:", truth_ms$g_0, "\n")
cat("designed G_t:", truth_ms$g_t, "\n")

# ============================================================================
# 1. gaussian DBN
# ============================================================================
section("gaussian DBN")

# named CPD vector in the format expected by dbn.fit(CPDs = ...):
# (Intercept), coefficients in the stored parents order, Std (res)
make_cpd <- function(dbn, node, intercept, coefs, sd) {
  parents <- get_parent_set(dbn, node)
  stopifnot(setequal(names(coefs), parents))
  values <- c(intercept, as.numeric(coefs[parents]), sd)
  names(values) <- c("(Intercept)", parents, "Std (res)")
  values
}

cpds <- list(
  A_0 = make_cpd(truth, "A_0",  0.0, c(), 1.0),
  B_0 = make_cpd(truth, "B_0",  0.0, c(), 1.0),
  C_0 = make_cpd(truth, "C_0",  1.0, c(A_0 = 0.8, B_0 = -0.7), 0.5),
  A_t = make_cpd(truth, "A_t",  0.0, c(`A_t-1` = 0.6), 0.4),
  B_t = make_cpd(truth, "B_t",  0.2, c(`B_t-1` = 0.5), 0.4),
  C_t = make_cpd(truth, "C_t", -0.3, c(`C_t-1` = 0.4, A_t = 0.7, B_t = -0.6), 0.3)
)
fit_gaussian <- dbn.fit(DBN = truth, CPDs = cpds)

# --- equivalence with dbn.sampling under the same seed ---------------------
n_samples <- 800
horizon <- 8
set.seed(20260610)
reference <- dbn.sampling(fit_gaussian, n_samples, horizon)
set.seed(20260610)
sampled <- dbn.sampling.cpp(fit_gaussian, n_samples, horizon)

if (!identical(reference, sampled)) {
  cat("identical() is FALSE, all.equal():\n")
  print(all.equal(reference, sampled))
}
check("gaussian: same dataset as dbn.sampling under the same seed",
      identical(reference, sampled) || isTRUE(all.equal(reference, sampled)))

# --- structure recovery with hill climbing ---------------------------------
learned_ms <- modelstring(hc.dbn(sampled, markov_order = 1))
cat("hc G_0:", learned_ms$g_0, "\n")
cat("hc G_t:", learned_ms$g_t, "\n")
check("gaussian: hc.dbn reconstructs G_0", same_structure(truth_ms$g_0, learned_ms$g_0))
check("gaussian: hc.dbn reconstructs G_t", same_structure(truth_ms$g_t, learned_ms$g_t))

# ============================================================================
# 2. discrete DBN
# ============================================================================
section("discrete DBN")

# CPT array in the format expected by dbn.fit(CPTs = ...):
# P(node = "yes" | parents) = base_p + sum of bumps[parent] over parents at "yes"
make_cpt <- function(dbn, node, base_p, bumps = c(), levels = c("no", "yes")) {
  parents <- get_parent_set(dbn, node)
  stopifnot(setequal(names(bumps), parents))
  if (length(parents) == 0) {
    probs <- c(1 - base_p, base_p)
  } else {
    # first parent varies fastest, matching the column-major CPT layout
    combos <- expand.grid(rep(list(0:1), length(parents)))
    probs <- numeric(0)
    for (r in seq_len(nrow(combos))) {
      p_yes <- base_p + sum(as.numeric(bumps[parents]) * as.numeric(unlist(combos[r, ])))
      probs <- c(probs, 1 - p_yes, p_yes)
    }
  }
  dn <- c(list(levels), rep(list(levels), length(parents)))
  names(dn) <- c(node, parents)
  array(probs, dim = rep(2L, length(parents) + 1L), dimnames = dn)
}

cpts <- list(
  A_0 = make_cpt(truth, "A_0", 0.45),
  B_0 = make_cpt(truth, "B_0", 0.55),
  C_0 = make_cpt(truth, "C_0", 0.10, c(A_0 = 0.35, B_0 = 0.45)),
  A_t = make_cpt(truth, "A_t", 0.12, c(`A_t-1` = 0.72)),
  B_t = make_cpt(truth, "B_t", 0.15, c(`B_t-1` = 0.65)),
  C_t = make_cpt(truth, "C_t", 0.05, c(`C_t-1` = 0.30, A_t = 0.30, B_t = 0.30))
)
fit_discrete <- dbn.fit(DBN = truth, CPTs = cpts)

# --- equivalence with dbn.sampling under the same seed ---------------------
n_samples_d <- 2000
set.seed(992)
reference_d <- dbn.sampling(fit_discrete, n_samples_d, horizon)
set.seed(992)
sampled_d <- dbn.sampling.cpp(fit_discrete, n_samples_d, horizon)
check("discrete: same dataset as dbn.sampling under the same seed",
      identical(reference_d, sampled_d))

# --- structure recovery with hill climbing ---------------------------------
learned_d_ms <- modelstring(hc.dbn(sampled_d, markov_order = 1))
cat("hc G_0:", learned_d_ms$g_0, "\n")
cat("hc G_t:", learned_d_ms$g_t, "\n")
check("discrete: hc.dbn reconstructs G_0", same_structure(truth_ms$g_0, learned_d_ms$g_0))
check("discrete: hc.dbn reconstructs G_t", same_structure(truth_ms$g_t, learned_d_ms$g_t))

# ============================================================================
# 3. timing
# ============================================================================
section("timing (n_samples = 400, max_time = 10)")

elapsed_r <- system.time(dbn.sampling(fit_gaussian, 400, 10))["elapsed"]
elapsed_cpp <- system.time(dbn.sampling.cpp(fit_gaussian, 400, 10))["elapsed"]
cat(sprintf("gaussian  | dbn.sampling: %7.3fs | dbn.sampling.cpp: %7.3fs | speed-up: x%.0f\n",
            elapsed_r, elapsed_cpp, elapsed_r / max(elapsed_cpp, 1e-3)))

elapsed_r <- system.time(dbn.sampling(fit_discrete, 400, 10))["elapsed"]
elapsed_cpp <- system.time(dbn.sampling.cpp(fit_discrete, 400, 10))["elapsed"]
cat(sprintf("discrete  | dbn.sampling: %7.3fs | dbn.sampling.cpp: %7.3fs | speed-up: x%.0f\n",
            elapsed_r, elapsed_cpp, elapsed_r / max(elapsed_cpp, 1e-3)))

cat("\nall checks passed\n")

# ============================================================================
# 4. profiling
#
# bench::mark() runs each expression repeatedly and reports the distribution of
# the timings (median, min, etc.) together with memory allocation. The two
# expressions are NOT checked for equality here (check = FALSE): under the same
# seed they return the same dataset, but bench::mark() re-runs them many times
# without resetting the seed, so the per-iteration results legitimately differ.
# ============================================================================
section("profiling (bench::mark)")

if (!requireNamespace("bench", quietly = TRUE)) {
  cat("package 'bench' is not installed; skipping profiling section\n")
} else {
  bench_n_samples <- 400
  bench_max_time <- 10
  cat(sprintf("n_samples = %d, max_time = %d\n\n", bench_n_samples, bench_max_time))

  cat("-- gaussian --\n")
  bench_gaussian <- bench::mark(
    dbn.sampling     = dbn.sampling(fit_gaussian, bench_n_samples, bench_max_time),
    dbn.sampling.cpp = dbn.sampling.cpp(fit_gaussian, bench_n_samples, bench_max_time),
    check = FALSE,
    min_iterations = 10
  )
  print(bench_gaussian[, c("expression", "min", "median", "itr/sec", "mem_alloc")])

  cat("\n-- discrete --\n")
  bench_discrete <- bench::mark(
    dbn.sampling     = dbn.sampling(fit_discrete, bench_n_samples, bench_max_time),
    dbn.sampling.cpp = dbn.sampling.cpp(fit_discrete, bench_n_samples, bench_max_time),
    check = FALSE,
    min_iterations = 10
  )
  print(bench_discrete[, c("expression", "min", "median", "itr/sec", "mem_alloc")])

  # speed-up: how many times faster the C++ version is, on the min and the
  # median timings (bench stores them as <bench_time> objects in seconds; row 1
  # is dbn.sampling, row 2 is dbn.sampling.cpp)
  report_speedup <- function(label, res) {
    min_r    <- as.numeric(res$min[1])
    min_cpp  <- as.numeric(res$min[2])
    med_r    <- as.numeric(res$median[1])
    med_cpp  <- as.numeric(res$median[2])
    cat(sprintf("%-9s | speed-up min: x%.0f | speed-up median: x%.0f\n",
                label, min_r / min_cpp, med_r / med_cpp))
  }

  cat("\n-- speed-up (dbn.sampling / dbn.sampling.cpp) --\n")
  report_speedup("gaussian", bench_gaussian)
  report_speedup("discrete", bench_discrete)
}
