library(DynamicBayesianNetwork)

# =============================================================================
# 1. Define a DBN structure manually
# =============================================================================

# Empty DBN with 3 dynamic nodes, Markov order 1
dbn_manual <- empty.dbn(dynamic_nodes = c("A", "B", "C"), markov_order = 1)

# Intra-slice arcs (t -> t) and inter-slice arcs (t-k -> t)
dbn_manual <- add.arc.dbn(dbn_manual, from = c("A", "t"),   to = c("B", "t"))    # A_t   -> B_t
dbn_manual <- add.arc.dbn(dbn_manual, from = c("A", "t-1"), to = c("A", "t"))    # A_t-1 -> A_t
dbn_manual <- add.arc.dbn(dbn_manual, from = c("B", "t-1"), to = c("C", "t"))    # B_t-1 -> C_t

print(dbn_manual)
modelstring.dbn(dbn_manual)   # list(g_0 = "...", g_t = "...")


# =============================================================================
# 2. Generate a random DBN structure
# =============================================================================

set.seed(42)

# Random transition network. With g_0_arcs = TRUE (default), G_0 mirrors the
# intra-slice arcs of G_t. Pass per-lag probabilities for markov_order > 1.
true_dbn <- random.structure.dbn(
  node_names           = c("A", "B", "C", "D"),
  prob_edge_intraslice = 0.4,
  prob_edge_interslice = 0.4,
  markov_order         = 1
)

print(true_dbn)
modelstring.dbn(true_dbn)

# Visualise
plot_g0(true_dbn)
plot_g_transition(true_dbn)


# =============================================================================
# 3. Fit random parameters — discrete (Dirichlet-sampled CPTs)
# =============================================================================

fitted_discrete <- fit_random_dbn(
  generated_dbn             = true_dbn,
  type                      = "discrete",
  max_variables_cardinality = 2,           # binary variables
  fixed_cardinality         = TRUE
)
print(fitted_discrete)


# =============================================================================
# 4. Fit random parameters — continuous (gaussian regressions)
# =============================================================================

# Coefficients ~ U(param.lower, param.upper); same residual sd for every node
fitted_gaussian <- fit_random_dbn(
  generated_dbn = true_dbn,
  type          = "continuous",
  param.lower   = 0.2,
  param.upper   = 0.8,
  sd            = 1
)
print(fitted_gaussian)

# (fit_random_dbn(type = "mixed") -> stop("to be implemented"))


# =============================================================================
# 5. Sampling
# =============================================================================

# Discrete sample: 300 trajectories of length 6
data_discrete <- dbn.sampling(fitted_discrete, n_samples = 300, max_time = 6)
head(data_discrete)

# Gaussian sample: same shape, numeric columns
data_gaussian <- dbn.sampling(fitted_gaussian, n_samples = 300, max_time = 6)
head(data_gaussian)


# =============================================================================
# 6. Forecasting
# =============================================================================

# Last observed time-step of the first sample as the seed.
# dbn.sampling writes Sample_id as "sample1", "sample2", ... — pick whichever
# happens to be first to stay robust across runs.
first_sample <- unique(data_gaussian$Sample_id)[1]
seed_obs <- data_gaussian[data_gaussian$Sample_id == first_sample &
                          data_gaussian$Time      == max(data_gaussian$Time), ]

forecast <- dbn.forecasting(fitted_gaussian, observations = seed_obs, timepoints = 5)
print(forecast)


# =============================================================================
# 7. Structure learning — recover the structure from the sampled data
# =============================================================================

# --- Score-based (Hill Climbing) ---
learned_hc  <- hc.dbn(data_gaussian)

# --- Constraint-based (PC stable) ---
learned_pc  <- pc.stable.dbn(data_gaussian)

# --- Hybrid (Max-Min Hill Climbing) ---
learned_mmhc <- mmhc.dbn(data_gaussian)

# High-level wrapper: auto-selects defaults, handles g0 + transition jointly
learned <- dbn.learn.structure(data_gaussian)


# =============================================================================
# 8. Compare learned vs ground-truth structures (modelstring.dbn)
# =============================================================================

# modelstring.dbn returns a list with the G_0 and G_t modelstrings in the
# bnlearn DSL — easy to eyeball-diff against the ground truth.
truth_ms   <- modelstring.dbn(true_dbn)
hc_ms      <- modelstring.dbn(learned_hc)
pc_ms      <- modelstring.dbn(learned_pc)
mmhc_ms    <- modelstring.dbn(learned_mmhc)
wrapper_ms <- modelstring.dbn(learned)

cat("---- G_0 ----\n")
cat("truth   :", truth_ms$g_0,   "\n")
cat("hc      :", hc_ms$g_0,      "\n")
cat("pc.stab :", pc_ms$g_0,      "\n")
cat("mmhc    :", mmhc_ms$g_0,    "\n")
cat("wrapper :", wrapper_ms$g_0, "\n")

cat("---- G_t ----\n")
cat("truth   :", truth_ms$g_t,   "\n")
cat("hc      :", hc_ms$g_t,      "\n")
cat("pc.stab :", pc_ms$g_t,      "\n")
cat("mmhc    :", mmhc_ms$g_t,    "\n")
cat("wrapper :", wrapper_ms$g_t, "\n")

# Identical structures iff the modelstring lists are equal
identical(truth_ms, hc_ms)


# =============================================================================
# 9. Parameter learning from data
# =============================================================================

# Refit gaussian parameters on the learned structure and compare to the
# generating model using modelstring.dbn.fit (same DSL, on a fitted DBN).
refit <- dbn.fit(DBN = learned_hc, data = data_gaussian)
print(refit)

modelstring.dbn.fit(refit)
modelstring.dbn.fit(fitted_gaussian)


# =============================================================================
# 10. Low-level structure assembly (dbn.build.struct)
# =============================================================================

# dbn.prep.data + individual bnlearn fits + dbn.build.struct
c(data0, dataTR, bl0, blTR, wl0, wlTR) %<-% dbn.prep.data(data_gaussian)

bn_g0 <- bnlearn::hc(data0,  blacklist = bl0)
bn_tr <- bnlearn::hc(dataTR, blacklist = blTR)

dbn_assembled <- dbn.build.struct(PN = bn_g0, TN = bn_tr, markov_order = 1)
print(dbn_assembled)
modelstring.dbn(dbn_assembled)


# =============================================================================
# 11. End-to-end with Markov order 3 (random.dbn wrapper)
# =============================================================================

# random.dbn is the wrapper around random.structure.dbn; the same intra-slice
# / inter-slice / G_0 knobs as before, plus markov_order. Pass a length-3
# vector for prob_edge_interslice to assign a per-lag probability:
# slot 1 = lag t-1, slot 2 = lag t-2, slot 3 = lag t-3 (sparser as lag grows).
set.seed(7)
mo3_dbn <- random.dbn(
  nodes_names          = c("A", "B", "C"),
  is_same              = TRUE,
  prob_edge_intraslice = 0.4,
  prob_edge_interslice = c(0.5, 0.3, 0.15),
  markov_order         = 3
)

print(mo3_dbn)
modelstring.dbn(mo3_dbn)

# Fit random gaussian parameters. fit_random_dbn_g adds a `priors` sub-list
# with marginal intercepts for variable_1, variable_2 (the initial slices
# needed to bootstrap a markov_order > 1 chain).
mo3_fitted <- fit_random_dbn(mo3_dbn, type = "continuous",
                             param.lower = 0.2, param.upper = 0.6, sd = 0.5)

# Sample a longer trajectory so each variable has enough history at every lag
mo3_data <- dbn.sampling(mo3_fitted, n_samples = 5000, max_time = 10)
head(mo3_data)

# Re-learn the structure. hc.dbn (and the other learners) take a markov_order
# argument; pass 3 so the transition network considers lags up to t-3.
mo3_learned <- hc.dbn(mo3_data, markov_order = 3)

# Side-by-side comparison via modelstring.dbn — same DSL as before, but G_t now
# references nodes at t, t-1, t-2, t-3.
mo3_truth_ms   <- modelstring.dbn(mo3_dbn)
mo3_learned_ms <- modelstring.dbn(mo3_learned)

cat("---- MO=3 G_0 ----\n")
cat("truth   :", mo3_truth_ms$g_0,   "\n")
cat("hc      :", mo3_learned_ms$g_0, "\n")

cat("---- MO=3 G_t ----\n")
cat("truth   :", mo3_truth_ms$g_t,   "\n")
cat("hc      :", mo3_learned_ms$g_t, "\n")

# is the transition network identical?
identical(mo3_truth_ms, mo3_learned_ms)
