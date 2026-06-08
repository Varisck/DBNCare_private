# Tests for structure learning in the continuous (Gaussian) case.
#
# Build a Gaussian DBN by hand, sample a long dataset from it, then learn the
# structure back from the data and check that the recovered arcs match the
# ground-truth arcs (G_0 and G_transition).

library(testthat)

# ----- helpers --------------------------------------------------------------

# turn an arcs matrix (or NULL) into a sorted character vector of "from->to"
arc_set <- function(arcs) {
  if (is.null(arcs) || nrow(arcs) == 0) return(character(0))
  sort(paste(arcs[, "from"], arcs[, "to"], sep = "->"))
}

# extract the set of arcs in the recovered DBN scoped to time-0 / transition
arcs_in_dbn <- function(dbn_learned) {
  arc_set(dbn_learned$arcs)
}

# ----- fixture: ground-truth Gaussian DBN -----------------------------------

# Structure
#   G_0           B_0 -> A_0,  C_0 -> A_0
#   G_transition  B_t -> A_t,  B_t -> C_t
#                 A_t-1 -> A_t, B_t-1 -> B_t, C_t-1 -> C_t
build_ground_truth_dbn <- function() {
  d <- empty.dbn(dynamic_nodes = c("A", "B", "C"), markov_order = 1)

  d <- add.arc.dbn(d, from = c("A", "t-1"), to = c("A", "t"))
  d <- add.arc.dbn(d, from = c("B", "t-1"), to = c("B", "t"))
  d <- add.arc.dbn(d, from = c("B", "t"),   to = c("A", "t"))
  d <- add.arc.dbn(d, from = c("C", "t-1"), to = c("C", "t"))
  d <- add.arc.dbn(d, from = c("B", "t"),   to = c("C", "t"))

  d <- add.arc.dbn(d, from = c("B", "t_0"), to = c("A", "t_0"))
  d <- add.arc.dbn(d, from = c("C", "t_0"), to = c("A", "t_0"))

  d
}

# CPDs follow the format documented in learn_param_g_cpds:
#   c("(Intercept)", parents..., "Std (res)")
# parent ordering must equal the parent set order recorded in the DBN.
build_ground_truth_cpds <- function() {
  list(
    "A_0" = c(
      "(Intercept)" = 1,
      "B_0"         = 1,
      "C_0"         = 1,
      "Std (res)"   = .2
    ),
    "B_0" = c(
      "(Intercept)" = .5,
      "Std (res)"   = 1
    ),
    "C_0" = c(
      "(Intercept)" = 0,
      "Std (res)"   = 1
    ),
    "A_t" = c(
      "(Intercept)" = 1,
      "A_t-1"       = -1,
      "B_t"         = 0.5,
      "Std (res)"   = .1
    ),
    "B_t" = c(
      "(Intercept)" = 1,
      "B_t-1"       = 0.5,
      "Std (res)"   = 1
    ),
    "C_t" = c(
      "(Intercept)" = 1,
      "C_t-1"       = -1,
      "B_t"         = 0.5,
      "Std (res)"   = 0.2
    )
  )
}

# ----- shared sample --------------------------------------------------------

dbn_truth <- build_ground_truth_dbn()
cpds_truth <- build_ground_truth_cpds()
fitted_truth <- dbn.fit(DBN = dbn_truth, CPDs = cpds_truth)

# sanity: the auto-detected dbn type must be gaussian
test_that("ground-truth fitted DBN is gaussian", {
  expect_equal(dbn_type(fitted_truth), "gaussian")
})

# sample once and reuse across tests (keeps the suite fast).
# seed is set right before sampling so the dataset is pinned regardless of
# RNG state inherited from earlier test files.
set.seed(123)
sampled_data <- dbn.sampling(fitted_truth, n_samples = 500, max_time = 4)

# sanity: sampled data is numeric (so dataset_type should report gaussian) and
# has the expected columns
test_that("sampled dataset is numeric / gaussian", {
  expect_equal(dataset_type(sampled_data), "gaussian")
  expect_true(all(c("Sample_id", "Time", "A", "B", "C") %in% names(sampled_data)))
  expect_true(all(sapply(sampled_data[, c("A", "B", "C")], is.numeric)))
})

# ----- structure recovery ---------------------------------------------------

expected_arcs <- sort(c(
  "B_0->A_0", "C_0->A_0",
  "B_t->A_t", "B_t->C_t",
  "A_t-1->A_t", "B_t-1->B_t", "C_t-1->C_t"
))

# score-based algorithms (hc, tabu) are the natural fit for continuous data.
# Passing score = NULL lets bnlearn pick its default gaussian score (bic-g),
# which is what works on numeric columns.
test_that("hc recovers the ground-truth gaussian DBN structure", {
  set.seed(1)
  learned <- dbn.learn.structure(sampled_data, algorithm = "hc",
                                       score = NULL)
  expect_s3_class(learned, "dbn")
  expect_setequal(arcs_in_dbn(learned), expected_arcs)
})

test_that("tabu recovers the ground-truth gaussian DBN structure", {
  set.seed(2)
  learned <- dbn.learn.structure(sampled_data, algorithm = "tabu",
                                       score = NULL)
  expect_s3_class(learned, "dbn")
  expect_setequal(arcs_in_dbn(learned), expected_arcs)
})

# hybrid algorithms combine a constraint-based restrict phase with a
# score-based maximize phase. score = NULL / test = NULL lets bnlearn pick
# the gaussian-appropriate defaults for both phases.
test_that("mmhc recovers the ground-truth gaussian DBN structure", {
  set.seed(4)
  learned <- dbn.learn.structure(sampled_data, algorithm = "mmhc",
                                       score = NULL, test = NULL)
  expect_s3_class(learned, "dbn")
  expect_setequal(arcs_in_dbn(learned), expected_arcs)
})

# ----- sanity: parameter learning closes the loop ---------------------------

test_that("parameters re-learned from the gaussian sample are close to truth", {
  set.seed(4)
  learned_struct <- dbn.learn.structure(sampled_data, algorithm = "hc",
                                              score = NULL)
  refit <- dbn.fit(DBN = learned_struct, data = sampled_data)
  expect_equal(class(refit), "dbn.fit")

  # intercepts and standard deviations should be close to the ground-truth
  # values (loose tolerances: this is a finite-sample check, not exact)
  expect_equal(unname(refit$B_0$regs["(Intercept)"]), 0.5, tolerance = 0.2)
  expect_equal(unname(refit$C_0$regs["(Intercept)"]), 0.0, tolerance = 0.2)
  expect_equal(unname(refit$B_0$std),                 1.0, tolerance = 0.2)
  expect_equal(unname(refit$A_t$std),                 0.1, tolerance = 0.1)
  expect_equal(unname(refit$B_t$std),                 1.0, tolerance = 0.2)
  expect_equal(unname(refit$C_t$std),                 0.2, tolerance = 0.1)
})
