# Tests for parameter learning in the mixed (conditional gaussian) case.
#
# Build a mixed DBN by hand (one discrete node, one gaussian node, one
# conditional-gaussian / CLG node) and fit it from a per-node distribution,
# then check the fit dispatched every node to the right subroutine and stored
# the parameters as given.
#
# NOTE: only the distribution path is exercised here. Learning mixed parameters
# from a data.frame is not implemented yet, so there is no data-path test.

library(testthat)

# ----- fixture: ground-truth mixed DBN --------------------------------------
#
# Structure (markov_order = 1)
#   G_0           D_0 -> M_0,  G_0 -> M_0
#   G_transition  D_t-1 -> D_t, G_t-1 -> G_t
#                 D_t   -> M_t, G_t   -> M_t
#
# Node types
#   D : discrete            (parent: D_t-1)
#   G : gaussian            (parent: G_t-1)
#   M : conditional gaussian / CLG
#         continuous parent: G   discrete parent: D
build_ground_truth_dbn <- function() {
  d <- empty.dbn(dynamic_nodes = c("D", "G", "M"), markov_order = 1)

  # transition slice
  d <- add.arc.dbn(d, from = c("D", "t-1"), to = c("D", "t"))
  d <- add.arc.dbn(d, from = c("G", "t-1"), to = c("G", "t"))
  d <- add.arc.dbn(d, from = c("D", "t"),   to = c("M", "t"))
  d <- add.arc.dbn(d, from = c("G", "t"),   to = c("M", "t"))

  # time-0 slice
  d <- add.arc.dbn(d, from = c("D", "t_0"), to = c("M", "t_0"))
  d <- add.arc.dbn(d, from = c("G", "t_0"), to = c("M", "t_0"))

  d
}

# A CLG distribution is a list(coef, sd):
#   coef : array whose first dimension holds the gaussian regressors
#          c("(Intercept)", continuous parents...) and whose remaining
#          dimensions are indexed by the discrete parents' levels.
#   sd   : numeric vector, one residual sd per discrete-parent combination.
build_ground_truth_dist <- function() {
  lv <- c("yes", "no")

  # --- time-0 slice ---

  # D_0: discrete, no parents
  D_0 <- array(c(0.4, 0.6), dim = 2, dimnames = list(D_0 = lv))

  # G_0: gaussian, no parents -> intercept + std
  G_0 <- c("(Intercept)" = 0, "Std (res)" = 1)

  # M_0: CLG, continuous parent G_0, discrete parent D_0
  M_0_coef <- array(
    c(1, 0.5,    # D_0 = "yes": intercept = 1,   slope on G_0 = 0.5
      2, 0.7),   # D_0 = "no" : intercept = 2,   slope on G_0 = 0.7
    dim = c(2, 2),
    dimnames = list(M_0 = c("(Intercept)", "G_0"), D_0 = lv)
  )
  M_0 <- list(coef = M_0_coef, sd = c(0.2, 0.3))

  # --- transition slice ---

  # D_t: discrete, parent D_t-1
  dims_D_t <- list(D_t = lv)
  dims_D_t[["D_t-1"]] <- lv
  D_t <- array(
    c(0.8, 0.2,   # D_t | D_t-1 = "yes"
      0.3, 0.7),  # D_t | D_t-1 = "no"
    dim = c(2, 2),
    dimnames = dims_D_t
  )

  # G_t: gaussian, parent G_t-1
  G_t <- c("(Intercept)" = 1, "G_t-1" = 0.5, "Std (res)" = 1)

  # M_t: CLG, continuous parent G_t, discrete parent D_t
  dims_M_t <- list(M_t = c("(Intercept)", "G_t"))
  dims_M_t[["D_t"]] <- lv
  M_t_coef <- array(
    c(1, -1,     # D_t = "yes": intercept = 1,  slope on G_t = -1
      0, 0.5),   # D_t = "no" : intercept = 0,  slope on G_t = 0.5
    dim = c(2, 2),
    dimnames = dims_M_t
  )
  M_t <- list(coef = M_t_coef, sd = c(0.1, 0.2))

  list(D_0 = D_0, G_0 = G_0, M_0 = M_0,
       D_t = D_t, G_t = G_t, M_t = M_t)
}

dbn_truth  <- build_ground_truth_dbn()
dist_truth <- build_ground_truth_dist()

# ----- 1. fit from the distribution and check the dispatch ------------------

test_that("dbn.fit(distribution = ...) produces a mixed dbn.fit", {
  fitted <- dbn.fit(DBN = dbn_truth, distribution = dist_truth)
  expect_equal(class(fitted), "dbn.fit")
  expect_equal(dbn_type(fitted), "mixed")
  expect_setequal(names(fitted),
                  c("D_0", "G_0", "M_0", "D_t", "G_t", "M_t"))
})

test_that("each node is dispatched to the correct node type / class", {
  fitted <- dbn.fit(DBN = dbn_truth, distribution = dist_truth)

  # discrete nodes -> dnode (carry $prob)
  for (variable in c("D_0", "D_t")) {
    expect_s3_class(fitted[[variable]], "dbn.fit.dnode")
    expect_false(is.null(fitted[[variable]]$prob))
  }

  # gaussian nodes -> gnode (carry $regs / $std)
  for (variable in c("G_0", "G_t")) {
    expect_s3_class(fitted[[variable]], "dbn.fit.gnode")
    expect_false(is.null(fitted[[variable]]$regs))
    expect_false(is.null(fitted[[variable]]$std))
  }

  # CLG nodes -> cgnode (carry $coefficients / $sd / $dlevels)
  for (variable in c("M_0", "M_t")) {
    expect_s3_class(fitted[[variable]], "dbn.fit.cgnode")
    expect_false(is.null(fitted[[variable]]$coefficients))
    expect_false(is.null(fitted[[variable]]$sd))
    expect_false(is.null(fitted[[variable]]$dlevels))
  }
})

# ----- 2. check stored parameters match the input distribution --------------

test_that("discrete CPTs stored on the fit match the input arrays", {
  fitted <- dbn.fit(DBN = dbn_truth, distribution = dist_truth)

  for (variable in c("D_0", "D_t")) {
    cpt <- dist_truth[[variable]]
    # probabilities preserved (aperm only reorders the dimensions)
    expect_equal(sum(fitted[[variable]]$prob), sum(cpt))
    expect_setequal(dimnames(fitted[[variable]]$prob)[[variable]],
                    dimnames(cpt)[[variable]])
  }
})

test_that("gaussian regs/std stored on the fit match the input CPDs", {
  fitted <- dbn.fit(DBN = dbn_truth, distribution = dist_truth)

  for (variable in c("G_0", "G_t")) {
    cpd     <- dist_truth[[variable]]
    regs_in <- cpd[seq_len(length(cpd) - 1)]   # intercept + parents
    std_in  <- cpd[length(cpd)]                # last entry

    expect_equal(fitted[[variable]]$regs, regs_in,
                 info = paste("regs mismatch for", variable))
    expect_equal(fitted[[variable]]$std, std_in,
                 info = paste("std mismatch for", variable))
  }
})

test_that("CLG coefficients/sd/levels stored on the fit match the input", {
  fitted <- dbn.fit(DBN = dbn_truth, distribution = dist_truth)

  for (variable in c("M_0", "M_t")) {
    dist <- dist_truth[[variable]]

    # gaussian regressors (intercept + continuous parents) preserved
    expect_equal(dimnames(fitted[[variable]]$coefficients)[[1]],
                 dimnames(dist$coef)[[1]],
                 info = paste("regressor labels mismatch for", variable))
    # one residual sd per discrete-parent combination
    expect_equal(fitted[[variable]]$sd, dist$sd,
                 info = paste("sd mismatch for", variable))
    # discrete-parent levels recorded
    expect_equal(unname(fitted[[variable]]$dlevels),
                 unname(dimnames(dist$coef)[-1]),
                 info = paste("dlevels mismatch for", variable))
  }
})
