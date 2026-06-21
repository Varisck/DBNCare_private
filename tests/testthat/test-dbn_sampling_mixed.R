# Tests for sampling a mixed (conditional-gaussian) DBN.
#
# Build a markov-order-1 mixed DBN by hand (one discrete, one gaussian, one CLG
# node), fit it from a known distribution, sample a long dataset, then close the
# loop two ways:
#   1. parameter learning on the ground-truth structure recovers the parameters;
#   2. structure learning from the sample recovers the ground-truth arcs.

library(testthat)

# absolute-tolerance closeness check with an informative failure message
expect_close <- function(object, expected, tol, label = "") {
  testthat::expect_true(
    isTRUE(abs(object - expected) < tol),
    info = sprintf("%s: |%.4f - %.4f| = %.4f (tol %.3f)",
                   label, object, expected, abs(object - expected), tol)
  )
}

# turn an arcs matrix (or NULL) into a sorted "from->to" character vector
arc_set <- function(arcs) {
  if (is.null(arcs) || nrow(arcs) == 0) return(character(0))
  sort(paste(arcs[, "from"], arcs[, "to"], sep = "->"))
}

# ----- fixture: ground-truth mixed DBN --------------------------------------
#
# Structure (markov_order = 1)
#   G_0           D_0 -> M_0,  G_0 -> M_0
#   G_transition  D_t-1 -> D_t, G_t-1 -> G_t
#                 D_t   -> M_t, G_t   -> M_t
#
#   D : discrete (yes/no)        G : gaussian        M : CLG (cont G, disc D)
build_ground_truth_dbn <- function() {
  d <- empty.dbn(dynamic_nodes = c("D", "G", "M"), markov_order = 1)

  d <- add.arc.dbn(d, from = c("D", "t-1"), to = c("D", "t"))
  d <- add.arc.dbn(d, from = c("G", "t-1"), to = c("G", "t"))
  d <- add.arc.dbn(d, from = c("D", "t"),   to = c("M", "t"))
  d <- add.arc.dbn(d, from = c("G", "t"),   to = c("M", "t"))

  d <- add.arc.dbn(d, from = c("D", "t_0"), to = c("M", "t_0"))
  d <- add.arc.dbn(d, from = c("G", "t_0"), to = c("M", "t_0"))

  d
}

# CLG generating parameters, keyed by the discrete-parent level, so the same
# numbers drive both the fitted distribution and the recovery checks.
lv <- c("yes", "no")

m0_intercept <- c(yes =  1.0, no = -1.0)
m0_slope     <- c(yes =  2.0, no = -1.0)   # slope on G_0
m0_sd        <- c(yes =  0.4, no =  0.6)

mt_intercept <- c(yes =  2.0, no = -2.0)
mt_slope     <- c(yes =  1.5, no = -0.5)   # slope on G_t
mt_sd        <- c(yes =  0.3, no =  0.5)

build_ground_truth_dist <- function() {
  # --- time-0 slice ---
  D_0 <- array(c(0.5, 0.5), dim = 2, dimnames = list(D_0 = lv))
  G_0 <- c("(Intercept)" = 0, "Std (res)" = 2)        # N(0, 2)
  M_0_coef <- array(
    c(m0_intercept[["yes"]], m0_slope[["yes"]],
      m0_intercept[["no"]],  m0_slope[["no"]]),
    dim = c(2, 2),
    dimnames = list(M_0 = c("(Intercept)", "G_0"), D_0 = lv)
  )
  M_0 <- list(coef = M_0_coef, sd = c(m0_sd[["yes"]], m0_sd[["no"]]))

  # --- transition slice ---
  dims_D_t <- list(D_t = lv); dims_D_t[["D_t-1"]] <- lv
  D_t <- array(
    c(0.8, 0.2,    # D_t | D_t-1 = "yes"
      0.3, 0.7),   # D_t | D_t-1 = "no"
    dim = c(2, 2), dimnames = dims_D_t
  )
  G_t <- c("(Intercept)" = 0, "G_t-1" = 0.6, "Std (res)" = 1)   # AR(1)
  dims_M_t <- list(M_t = c("(Intercept)", "G_t")); dims_M_t[["D_t"]] <- lv
  M_t_coef <- array(
    c(mt_intercept[["yes"]], mt_slope[["yes"]],
      mt_intercept[["no"]],  mt_slope[["no"]]),
    dim = c(2, 2), dimnames = dims_M_t
  )
  M_t <- list(coef = M_t_coef, sd = c(mt_sd[["yes"]], mt_sd[["no"]]))

  list(D_0 = D_0, G_0 = G_0, M_0 = M_0,
       D_t = D_t, G_t = G_t, M_t = M_t)
}

# column index (in expand.grid(dlevels) order) of a fitted CLG node for a given
# discrete-parent level; here there is a single discrete parent, so it is just
# the position of that level.
clg_col <- function(node, d_level) {
  combos <- expand.grid(lapply(node$dlevels, as.character),
                        stringsAsFactors = FALSE)
  which(combos[[1]] == d_level)
}

# ----- shared fit + sample --------------------------------------------------

dbn_truth  <- build_ground_truth_dbn()
dist_truth <- build_ground_truth_dist()
fitted_truth <- dbn.fit(DBN = dbn_truth, distribution = dist_truth)

test_that("ground-truth fitted DBN is mixed", {
  expect_equal(dbn_type(fitted_truth), "mixed")
})

set.seed(424242)
sampled <- dbn.sampling(fitted_truth, n_samples = 2000, max_time = 5)

test_that("sampled mixed dataset has the right shape and column types", {
  expect_equal(dataset_type(sampled), "mixed")
  expect_true(all(c("Sample_id", "Time", "D", "G", "M") %in% names(sampled)))
  expect_true(is.character(sampled$D) || is.factor(sampled$D))
  expect_setequal(unique(as.character(sampled$D)), lv)
  expect_true(is.numeric(sampled$G))
  expect_true(is.numeric(sampled$M))
  expect_equal(nrow(sampled), 2000 * (5 + 1))
})

# ----- 1. parameter learning on the ground-truth structure ------------------

test_that("parameters re-learned from the mixed sample are close to truth", {
  refit <- dbn.fit(DBN = dbn_truth, data = sampled)
  expect_equal(class(refit), "dbn.fit")
  expect_equal(dbn_type(refit), "mixed")

  # discrete transition CPT  P(D_t | D_t-1)
  expect_close(refit$D_t$prob["yes", "yes"], 0.8, 0.05, "P(D_t=yes|yes)")
  expect_close(refit$D_t$prob["no",  "yes"], 0.2, 0.05, "P(D_t=no|yes)")
  expect_close(refit$D_t$prob["yes", "no"],  0.3, 0.05, "P(D_t=yes|no)")
  expect_close(refit$D_t$prob["no",  "no"],  0.7, 0.05, "P(D_t=no|no)")

  # gaussian AR(1) node G_t (regs matched positionally: intercept, lag coef)
  regs_t <- unname(refit$G_t$regs)
  expect_close(regs_t[1], 0.0, 0.1, "G_t intercept")
  expect_close(regs_t[2], 0.6, 0.1, "G_t lag coef")
  expect_close(unname(refit$G_t$std), 1.0, 0.1, "G_t sd")

  # CLG transition node M_t: one (intercept, slope, sd) per D_t level
  for (d in lv) {
    i <- clg_col(refit$M_t, d)
    expect_close(refit$M_t$coefficients["(Intercept)", i], mt_intercept[[d]],
                 0.1, paste0("M_t intercept @ D=", d))
    expect_close(refit$M_t$coefficients["G_t", i], mt_slope[[d]],
                 0.1, paste0("M_t slope @ D=", d))
    expect_close(refit$M_t$sd[i], mt_sd[[d]], 0.05, paste0("M_t sd @ D=", d))
  }

  # CLG time-0 node M_0 (fit on n_samples rows only -> looser tolerance)
  for (d in lv) {
    i <- clg_col(refit$M_0, d)
    expect_close(refit$M_0$coefficients["(Intercept)", i], m0_intercept[[d]],
                 0.15, paste0("M_0 intercept @ D=", d))
    expect_close(refit$M_0$coefficients["G_0", i], m0_slope[[d]],
                 0.15, paste0("M_0 slope @ D=", d))
    expect_close(refit$M_0$sd[i], m0_sd[[d]], 0.1, paste0("M_0 sd @ D=", d))
  }
})

# ----- 2. structure learning from the sample --------------------------------

expected_arcs <- sort(c(
  "D_0->M_0",   "G_0->M_0",
  "D_t-1->D_t", "G_t-1->G_t",
  "D_t->M_t",   "G_t->M_t"
))

# hc with the default mixed score (bic-cg). The continuous->discrete arcs are
# structurally forbidden in a CG network, so the D->M edges are forced into M,
# and the v-structures D->M<-G make the G->M directions identifiable too.
test_that("hc recovers the ground-truth mixed DBN structure", {
  set.seed(1)
  learned <- dbn.learn.structure(sampled, algorithm = "hc", score = NULL)
  expect_s3_class(learned, "dbn")
  expect_setequal(arc_set(learned$arcs), expected_arcs)
})

test_that("tabu recovers the ground-truth mixed DBN structure", {
  set.seed(2)
  learned <- dbn.learn.structure(sampled, algorithm = "tabu", score = NULL)
  expect_s3_class(learned, "dbn")
  expect_setequal(arc_set(learned$arcs), expected_arcs)
})

# ----- 3. markov order 2 ----------------------------------------------------
#
# Order-2 dependence lives in the discrete and gaussian self-loops
# (D_t-1, D_t-2 -> D_t and G_t-1, G_t-2 -> G_t); the CLG node M keeps
# contemporaneous parents (G_t, D_t). This exercises the unavailable-parent
# paths at t = 1 (the discrete CPT pooled over the missing D_t-2, the skipped
# G_t-2 regressor) while staying clear of a lagged discrete CLG parent, which
# the mixed sampler does not support.
build_mo2_dbn <- function() {
  d <- empty.dbn(dynamic_nodes = c("D", "G", "M"), markov_order = 2)

  d <- add.arc.dbn(d, from = c("D", "t-1"), to = c("D", "t"))
  d <- add.arc.dbn(d, from = c("D", "t-2"), to = c("D", "t"))
  d <- add.arc.dbn(d, from = c("G", "t-1"), to = c("G", "t"))
  d <- add.arc.dbn(d, from = c("G", "t-2"), to = c("G", "t"))
  d <- add.arc.dbn(d, from = c("G", "t"),   to = c("M", "t"))
  d <- add.arc.dbn(d, from = c("D", "t"),   to = c("M", "t"))

  d <- add.arc.dbn(d, from = c("G", "t_0"), to = c("M", "t_0"))
  d <- add.arc.dbn(d, from = c("D", "t_0"), to = c("M", "t_0"))

  d
}

build_mo2_dist <- function() {
  D_0 <- array(c(0.5, 0.5), dim = 2, dimnames = list(D_0 = lv))
  G_0 <- c("(Intercept)" = 0, "Std (res)" = 2)
  M_0_coef <- array(
    c(m0_intercept[["yes"]], m0_slope[["yes"]],
      m0_intercept[["no"]],  m0_slope[["no"]]),
    dim = c(2, 2), dimnames = list(M_0 = c("(Intercept)", "G_0"), D_0 = lv)
  )
  M_0 <- list(coef = M_0_coef, sd = c(m0_sd[["yes"]], m0_sd[["no"]]))

  # D_t | D_t-1, D_t-2 : depends on both lags (dims D_t x D_t-1 x D_t-2)
  dims_D_t <- list(D_t = lv); dims_D_t[["D_t-1"]] <- lv; dims_D_t[["D_t-2"]] <- lv
  D_t <- array(
    c(0.9, 0.1,    # D_t-1 = yes, D_t-2 = yes
      0.6, 0.4,    # D_t-1 = no,  D_t-2 = yes
      0.4, 0.6,    # D_t-1 = yes, D_t-2 = no
      0.1, 0.9),   # D_t-1 = no,  D_t-2 = no
    dim = c(2, 2, 2), dimnames = dims_D_t
  )
  G_t <- c("(Intercept)" = 0, "G_t-1" = 0.5, "G_t-2" = 0.3, "Std (res)" = 1)
  dims_M_t <- list(M_t = c("(Intercept)", "G_t")); dims_M_t[["D_t"]] <- lv
  M_t_coef <- array(
    c(mt_intercept[["yes"]], mt_slope[["yes"]],
      mt_intercept[["no"]],  mt_slope[["no"]]),
    dim = c(2, 2), dimnames = dims_M_t
  )
  M_t <- list(coef = M_t_coef, sd = c(mt_sd[["yes"]], mt_sd[["no"]]))

  list(D_0 = D_0, G_0 = G_0, M_0 = M_0, D_t = D_t, G_t = G_t, M_t = M_t)
}

dbn_mo2     <- build_mo2_dbn()
fitted_mo2  <- dbn.fit(DBN = dbn_mo2, distribution = build_mo2_dist())

set.seed(20260621)
sampled_mo2 <- dbn.sampling(fitted_mo2, n_samples = 2000, max_time = 6)

test_that("markov-order-2 mixed sampling runs and is well-formed", {
  expect_equal(dbn_type(fitted_mo2), "mixed")
  expect_equal(dataset_type(sampled_mo2), "mixed")
  expect_true(all(c("Sample_id", "Time", "D", "G", "M") %in% names(sampled_mo2)))
  expect_setequal(unique(as.character(sampled_mo2$D)), lv)
  expect_true(is.numeric(sampled_mo2$G) && all(is.finite(sampled_mo2$G)))
  expect_true(is.numeric(sampled_mo2$M) && all(is.finite(sampled_mo2$M)))
  expect_equal(nrow(sampled_mo2), 2000 * (6 + 1))
})

test_that("parameters re-learned from the markov-order-2 mixed sample match truth", {
  refit <- dbn.fit(DBN = dbn_mo2, data = sampled_mo2)
  expect_equal(dbn_type(refit), "mixed")

  # discrete order-2 CPT  P(D_t | D_t-1, D_t-2)
  expect_close(refit$D_t$prob["yes", "yes", "yes"], 0.9, 0.05, "P(yes|yes,yes)")
  expect_close(refit$D_t$prob["yes", "no",  "yes"], 0.6, 0.05, "P(yes|no,yes)")
  expect_close(refit$D_t$prob["yes", "yes", "no"],  0.4, 0.05, "P(yes|yes,no)")
  expect_close(refit$D_t$prob["yes", "no",  "no"],  0.1, 0.05, "P(yes|no,no)")

  # gaussian AR(2): intercept ~ 0, sd ~ 1, lag coefs ~ {0.5, 0.3} (order-agnostic)
  regs_t <- unname(refit$G_t$regs)
  expect_close(regs_t[1], 0.0, 0.1, "G_t intercept")
  expect_close(max(regs_t[-1]), 0.5, 0.1, "G_t larger lag coef")
  expect_close(min(regs_t[-1]), 0.3, 0.1, "G_t smaller lag coef")
  expect_close(unname(refit$G_t$std), 1.0, 0.1, "G_t sd")

  # CLG node M_t: one (intercept, slope, sd) per D_t level
  for (d in lv) {
    i <- clg_col(refit$M_t, d)
    expect_close(refit$M_t$coefficients["(Intercept)", i], mt_intercept[[d]],
                 0.1, paste0("M_t intercept @ D=", d))
    expect_close(refit$M_t$coefficients["G_t", i], mt_slope[[d]],
                 0.1, paste0("M_t slope @ D=", d))
    expect_close(refit$M_t$sd[i], mt_sd[[d]], 0.05, paste0("M_t sd @ D=", d))
  }
})
