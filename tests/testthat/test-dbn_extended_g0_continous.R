# Tests for the "extended G_0" behaviour in the continuous (Gaussian) case.
#
# When markov_order > 1 and extend_g0 = TRUE (the default), the prior network
# G_0 no longer holds just the time-0 slice: it spans every initial slice
# var_0, var_1, ..., var_{markov_order-1}. Both structure learning and
# parameter learning must therefore learn / fit those extra nodes.
#
# This file builds a markov_order = 2 ground-truth Gaussian DBN, samples from
# it, and checks that:
#   1. dbn.prep.data exposes the extended G_0 frame (and the matching blacklist)
#      only when extend_g0 = TRUE;
#   2. structure learning recovers the expected extended-G_0 arcs;
#   3. parameter learning fits the extended-G_0 nodes with parameters close to
#      the ones that generated the data.

library(testthat)

# ----- helpers --------------------------------------------------------------

arc_set <- function(arcs) {
  if (is.null(arcs) || nrow(arcs) == 0) return(character(0))
  sort(paste(arcs[, "from"], arcs[, "to"], sep = "->"))
}

# G_0 arcs of a learned DBN are the ones whose head is an initial-slice node
# (var_0, var_1, ...); transition arcs point into var_t.
g0_arc_set <- function(dbn_learned) {
  arcs <- dbn_learned$arcs
  arc_set(arcs[grepl("_[0-9]+$", arcs[, "to"]), , drop = FALSE])
}

# ----- fixture: ground-truth markov_order = 2 Gaussian DBN -------------------
#
# Structure
#   G_0           B_0 -> A_0,  C_0 -> A_0          (A_0 is a collider)
#   G_transition  A_t-1 -> A_t, A_t-2 -> A_t, B_t-1 -> A_t
#                 B_t-1 -> B_t, C_t-1 -> C_t       (no intra-slice _t arcs)
#
# With this transition, the first two sampled slices (0 and 1) depend on each
# other only through forward-in-time lag-1 arcs (the lag-2 self arc on A is
# zeroed at t = 1), so the joint over slices {0, 1} - i.e. the extended G_0 -
# is exactly:
#   B_0 -> A_0, C_0 -> A_0           (time-0 prior)
#   A_0 -> A_1, B_0 -> A_1           (A_t-1, B_t-1 into A_t at t = 1)
#   B_0 -> B_1, C_0 -> C_1           (self lag-1 on B, C)
build_ground_truth_dbn <- function() {
  d <- empty.dbn(dynamic_nodes = c("A", "B", "C"), markov_order = 2)

  d <- add.arc.dbn(d, from = c("B", "t_0"), to = c("A", "t_0"))
  d <- add.arc.dbn(d, from = c("C", "t_0"), to = c("A", "t_0"))

  d <- add.arc.dbn(d, from = c("A", "t-1"), to = c("A", "t"))
  d <- add.arc.dbn(d, from = c("A", "t-2"), to = c("A", "t"))
  d <- add.arc.dbn(d, from = c("B", "t-1"), to = c("A", "t"))
  d <- add.arc.dbn(d, from = c("B", "t-1"), to = c("B", "t"))
  d <- add.arc.dbn(d, from = c("C", "t-1"), to = c("C", "t"))

  d
}

build_ground_truth_cpds <- function() {
  list(
    "A_0" = c("(Intercept)" = 0, "B_0" = 1,   "C_0" = 1,   "Std (res)" = 0.3),
    "B_0" = c("(Intercept)" = 0, "Std (res)" = 1),
    "C_0" = c("(Intercept)" = 0, "Std (res)" = 1),
    "A_t" = c("(Intercept)" = 0, "A_t-1" = 0.6, "A_t-2" = 0.3, "B_t-1" = 0.5,
              "Std (res)" = 0.2),
    "B_t" = c("(Intercept)" = 0, "B_t-1" = 0.7, "Std (res)" = 0.5),
    "C_t" = c("(Intercept)" = 0, "C_t-1" = 0.7, "Std (res)" = 0.5)
  )
}

dbn_truth    <- build_ground_truth_dbn()
cpds_truth   <- build_ground_truth_cpds()
fitted_truth <- dbn.fit(DBN = dbn_truth, distribution = cpds_truth)

# sample once and reuse across tests; seed pinned right before sampling so the
# dataset does not depend on RNG state inherited from earlier test files.
set.seed(123)
sampled_data <- dbn.sampling(fitted_truth, n_samples = 1500, max_time = 4)

# arcs of the extended G_0 we expect to recover (see fixture comment)
expected_g0_arcs <- sort(c(
  "B_0->A_0", "C_0->A_0",
  "A_0->A_1", "B_0->A_1",
  "B_0->B_1", "C_0->C_1"
))

# ----- 1. dbn.prep.data exposes the extended G_0 frame ----------------------

test_that("extend_g0 = TRUE pivots the initial slices into the G_0 frame", {
  prep_ext <- dbn.prep.data(sampled_data, markov_order = 2, extend_g0 = TRUE)
  df_0_ext <- prep_ext[[1]]

  # G_0 now spans slices 0 and 1 for every variable
  expect_setequal(names(df_0_ext),
                  c("A_0", "A_1", "B_0", "B_1", "C_0", "C_1"))
})

test_that("extend_g0 = FALSE keeps only the time-0 slice in the G_0 frame", {
  prep_flat <- dbn.prep.data(sampled_data, markov_order = 2, extend_g0 = FALSE)
  df_0_flat <- prep_flat[[1]]

  expect_setequal(names(df_0_flat), c("A_0", "B_0", "C_0"))
})

test_that("extended G_0 blacklist forbids backward-in-time arcs only", {
  prep_ext <- dbn.prep.data(sampled_data, markov_order = 2, extend_g0 = TRUE)
  bl_0 <- arc_set(prep_ext[[3]])

  # arcs running back in time (slice 1 -> slice 0) must be blacklisted
  expect_true("A_1->A_0" %in% bl_0)
  expect_true("B_1->B_0" %in% bl_0)
  expect_true("C_1->C_0" %in% bl_0)

  # forward-in-time arcs (slice 0 -> slice 1) must remain allowed
  expect_false("A_0->A_1" %in% bl_0)
  expect_false("B_0->B_1" %in% bl_0)
  expect_false("C_0->C_1" %in% bl_0)
})

# ----- 2. structure learning recovers the extended G_0 ----------------------

test_that("hc recovers the extended G_0 structure (markov_order = 2)", {
  set.seed(1)
  learned <- dbn.learn.structure(sampled_data, markov_order = 2,
                                 algorithm = "hc", score = NULL)
  expect_s3_class(learned, "dbn")

  # the prior network must include the time-1 slice nodes
  g0_net <- get.g0.net(learned)
  expect_setequal(names(g0_net$nodes),
                  c("A_0", "A_1", "B_0", "B_1", "C_0", "C_1"))

  # and its arcs must match the ground-truth extended-G_0 arcs
  expect_setequal(g0_arc_set(learned), expected_g0_arcs)
})

test_that("tabu recovers the extended G_0 structure (markov_order = 2)", {
  set.seed(2)
  learned <- dbn.learn.structure(sampled_data, markov_order = 2,
                                 algorithm = "tabu", score = NULL)
  expect_s3_class(learned, "dbn")
  expect_setequal(g0_arc_set(learned), expected_g0_arcs)
})

# ----- 3. parameter learning fits the extended G_0 --------------------------

test_that("parameters of the extended G_0 are fitted close to the truth", {
  set.seed(1)
  learned <- dbn.learn.structure(sampled_data, markov_order = 2,
                                 algorithm = "hc", score = NULL)

  refit <- dbn.fit(DBN = learned, data = sampled_data)
  expect_equal(class(refit), "dbn.fit")
  expect_equal(dbn_type(refit), "gaussian")

  # every initial-slice node must be fitted, including the time-1 ones
  expect_true(all(c("A_0", "A_1", "B_0", "B_1", "C_0", "C_1") %in% names(refit)))

  # time-0 prior: A_0 = B_0 + C_0
  expect_equal(unname(refit[["A_0"]]$regs[["B_0"]]), 1.0, tolerance = 0.1)
  expect_equal(unname(refit[["A_0"]]$regs[["C_0"]]), 1.0, tolerance = 0.1)
  expect_equal(unname(refit[["A_0"]]$std),           0.3, tolerance = 0.1)

  # time-1 slice reflects the lag-1 transition dynamics (lag-2 zeroed at t = 1)
  expect_equal(unname(refit[["A_1"]]$regs[["A_0"]]), 0.6, tolerance = 0.1)
  expect_equal(unname(refit[["A_1"]]$regs[["B_0"]]), 0.5, tolerance = 0.1)
  expect_equal(unname(refit[["A_1"]]$std),           0.2, tolerance = 0.1)

  expect_equal(unname(refit[["B_1"]]$regs[["B_0"]]), 0.7, tolerance = 0.1)
  expect_equal(unname(refit[["C_1"]]$regs[["C_0"]]), 0.7, tolerance = 0.1)
})
