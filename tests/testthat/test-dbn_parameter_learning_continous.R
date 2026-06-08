# Tests for parameter learning in the continuous (Gaussian) case.
#
# Build a Gaussian DBN by hand, fit it from CPDs and check the fit matches the
# input CPDs. Then sample a long dataset from the fitted DBN, re-learn the
# parameters from the data and check the recovered regs/std are close to the
# original CPDs.

library(testthat)

# ----- fixture: ground-truth Gaussian DBN -----------------------------------
# (mirrors test-dbn_structure_learning_continous.R)
#
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
      "Std (res)"   = .4
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

dbn_truth  <- build_ground_truth_dbn()
cpds_truth <- build_ground_truth_cpds()

# ----- 1. fit from CPDs and check the fit matches the CPDs ------------------

test_that("dbn.fit(CPDs) produces a gaussian dbn.fit", {
  fitted <- dbn.fit(DBN = dbn_truth, CPDs = cpds_truth)
  expect_equal(class(fitted), "dbn.fit")
  expect_equal(dbn_type(fitted), "gaussian")
  expect_setequal(names(fitted),
                  c("A_0", "B_0", "C_0", "A_t", "B_t", "C_t"))
})

test_that("regs and std stored in dbn.fit match the input CPDs exactly", {
  fitted <- dbn.fit(DBN = dbn_truth, CPDs = cpds_truth)

  for (variable in names(cpds_truth)) {
    cpd       <- cpds_truth[[variable]]
    regs_in   <- cpd[seq_len(length(cpd) - 1)]   # intercept + parents
    std_in    <- cpd[length(cpd)]                # last entry

    # regression coefficients (intercept + parents) stored verbatim
    expect_equal(fitted[[variable]]$regs, regs_in,
                 info = paste("regs mismatch for", variable))
    # std stored verbatim
    expect_equal(fitted[[variable]]$std, std_in,
                 info = paste("std mismatch for", variable))
    # parent set on the fit matches the CPD's parent labels
    expect_equal(fitted[[variable]]$parents,
                 names(regs_in)[-1],
                 info = paste("parents mismatch for", variable))
  }
})

# ----- 2. learn from sampled data and check params are close to truth -------

# pin the sample so the test is reproducible regardless of inherited RNG state
set.seed(123)
sampled_data <- dbn.sampling(
  dbn.fit(DBN = dbn_truth, CPDs = cpds_truth),
  n_samples = 500,
  max_time  = 4
)

test_that("parameters re-learned from gaussian sample are close to truth", {
  refit <- dbn.fit(DBN = dbn_truth, data = sampled_data)
  expect_equal(class(refit), "dbn.fit")
  expect_equal(dbn_type(refit), "gaussian")

  # tolerances are loose: this is a finite-sample MLE check, not exact.
  # std-only nodes (no parents) just have intercept + std.
  for (variable in names(cpds_truth)) {
    cpd     <- cpds_truth[[variable]]
    regs_in <- cpd[seq_len(length(cpd) - 1)]
    std_in  <- unname(cpd[length(cpd)])

    # regression coefficients close to the ground-truth values
    expect_equal(unname(refit[[variable]]$regs),
                 unname(regs_in),
                 tolerance = 0.1,
                 info = paste("regs off for", variable))

    # residual std close to the ground-truth std
    expect_equal(unname(refit[[variable]]$std),
                 std_in,
                 tolerance = 0.05,
                 info = paste("std off for", variable))
  }
})
