# Tests for get_unrolled_dbn across the three network types (discrete, gaussian,
# mixed/conditional-gaussian), with markov_order > 1, both values of
# extend_with_transition, and edge cases where a lagged parent must be dropped
# (including a node whose single parent is dropped).
#
# All fixtures are built from an exact `distribution =` input so the unrolled
# node distributions can be asserted against known values.

library(testthat)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

INT <- "(Intercept)"

# Compare an unrolled discrete node CPT `actual` to an `expected` array. `map`
# renames the expected array's dimnames-names to the unrolled node's names; the
# comparison is invariant to dimension order and to within-dimension level order.
same_cpt <- function(actual, expected, map) {
  names(dimnames(expected)) <- unname(map[names(dimnames(expected))])
  idx <- dimnames(expected)                       # dim -> levels (unrolled names)
  actual <- aperm(actual, names(idx))             # align dimension order
  actual <- do.call(`[`, c(list(actual), idx, list(drop = FALSE)))  # align levels
  # compare values only (actual is a 'table', expected a plain array)
  isTRUE(all.equal(as.vector(actual), as.vector(expected), tolerance = 1e-8))
}

# coefficient / sd for a single discrete-parent configuration of a fitted CG node
# (works for both dbn.fit.cgnode and bn.fit.cgnode: their $coefficients columns
# are both in expand.grid($dlevels) order).
cg_lookup <- function(node, cfg) {
  combos <- expand.grid(node$dlevels, stringsAsFactors = FALSE)
  keep <- rep(TRUE, nrow(combos))
  for (p in names(cfg)) keep <- keep & (combos[[p]] == cfg[[p]])
  j <- which(keep)
  cf <- node$coefficients[, j]
  names(cf) <- rownames(node$coefficients)   # preserve reg names even for a 1-row matrix
  list(coef = cf, sd = unname(node$sd[j]))
}

# ===========================================================================
# DISCRETE
# ===========================================================================

build_discrete_dbn <- function() {
  lv <- c("y", "n")
  d <- empty.dbn(dynamic_nodes = c("A", "B"), markov_order = 2)
  # extended G_0 (A_1 | A_0, B_1 | B_0)
  d <- add.arc.dbn(d, c("A", "t_0"), c("A", "t_1"))
  d <- add.arc.dbn(d, c("B", "t_0"), c("B", "t_1"))
  # transition: A has two lagged parents, B has a single lagged parent
  d <- add.arc.dbn(d, c("A", "t-1"), c("A", "t"))
  d <- add.arc.dbn(d, c("A", "t-2"), c("A", "t"))
  d <- add.arc.dbn(d, c("B", "t-2"), c("B", "t"))

  A_0 <- array(c(.3, .7), 2, dimnames = list(A_0 = lv))
  A_1 <- array(c(.6, .4, .55, .45), c(2, 2), dimnames = list(A_1 = lv, A_0 = lv))
  B_0 <- array(c(.2, .8), 2, dimnames = list(B_0 = lv))
  B_1 <- array(c(.5, .5, .35, .65), c(2, 2), dimnames = list(B_1 = lv, B_0 = lv))
  dA <- list(A_t = lv); dA[["A_t-1"]] <- lv; dA[["A_t-2"]] <- lv
  A_t <- array(c(.8, .2, .3, .7, .65, .35, .1, .9), c(2, 2, 2), dimnames = dA)
  dB <- list(B_t = lv); dB[["B_t-2"]] <- lv
  B_t <- array(c(.9, .1, .25, .75), c(2, 2), dimnames = dB)

  fit <- dbn.fit(d, distribution = list(A_0 = A_0, A_1 = A_1, B_0 = B_0,
                                        B_1 = B_1, A_t = A_t, B_t = B_t))
  list(fit = fit, A_t = A_t, B_t = B_t, A_1 = A_1, B_1 = B_1)
}

test_that("discrete mo=2, extend_with_transition=FALSE: extended G_0 + full-history transition", {
  f <- build_discrete_dbn()
  u <- get_unrolled_dbn(f$fit, 3)

  expect_setequal(bnlearn::nodes(u),
                  c("A_0", "A_1", "B_0", "B_1", "A_2", "B_2", "A_3", "B_3"))

  # extended G_0 copied verbatim
  expect_true(same_cpt(u$A_1$prob, f$A_1, c("A_1" = "A_1", "A_0" = "A_0")))
  expect_true(same_cpt(u$B_1$prob, f$B_1, c("B_1" = "B_1", "B_0" = "B_0")))

  # transition slices carry the full Markov history and the verbatim CPT
  expect_setequal(u$A_2$parents, c("A_1", "A_0"))
  expect_true(same_cpt(u$A_2$prob, f$A_t,
                       c("A_t" = "A_2", "A_t-1" = "A_1", "A_t-2" = "A_0")))
  expect_true(same_cpt(u$A_3$prob, f$A_t,
                       c("A_t" = "A_3", "A_t-1" = "A_2", "A_t-2" = "A_1")))
  expect_setequal(u$B_2$parents, "B_0")
  expect_true(same_cpt(u$B_2$prob, f$B_t, c("B_t" = "B_2", "B_t-2" = "B_0")))
})

test_that("discrete mo=2, extend_with_transition=TRUE: drop missing lag = marginalize; single dropped parent = root", {
  f <- build_discrete_dbn()
  u <- get_unrolled_dbn(f$fit, 3, extend_with_transition = TRUE)

  expect_setequal(bnlearn::nodes(u),
                  c("A_0", "B_0", "A_1", "B_1", "A_2", "B_2", "A_3", "B_3"))

  # A_1: A_t-2 has no history -> marginalized (normalized) out
  expect_setequal(u$A_1$parents, "A_0")
  expA1 <- apply(f$A_t, c("A_t", "A_t-1"), mean)
  expect_true(same_cpt(u$A_1$prob, expA1, c("A_t" = "A_1", "A_t-1" = "A_0")))
  expect_true(all(abs(colSums(u$A_1$prob) - 1) < 1e-9))

  # EDGE: B's single parent B_t-2 is dropped -> B_1 becomes a root (marginal CPT)
  expect_length(u$B_1$parents, 0)
  expB1 <- apply(f$B_t, "B_t", mean)
  expect_equal(as.vector(u$B_1$prob), as.vector(expB1))
  expect_equal(sum(u$B_1$prob), 1)

  # from slice 2 on the full history is available again
  expect_true(same_cpt(u$A_2$prob, f$A_t,
                       c("A_t" = "A_2", "A_t-1" = "A_1", "A_t-2" = "A_0")))
})

test_that("discrete mo=3: transition drops several missing lags (marginalized)", {
  lv <- c("y", "n")
  d <- empty.dbn(dynamic_nodes = "A", markov_order = 3)
  d <- add.arc.dbn(d, c("A", "t-1"), c("A", "t"))
  d <- add.arc.dbn(d, c("A", "t-2"), c("A", "t"))
  d <- add.arc.dbn(d, c("A", "t-3"), c("A", "t"))
  A_0 <- array(c(.5, .5), 2, dimnames = list(A_0 = lv))
  dA <- list(A_t = lv); dA[["A_t-1"]] <- lv; dA[["A_t-2"]] <- lv; dA[["A_t-3"]] <- lv
  set.seed(1)
  vals <- as.vector(vapply(1:8, function(i) { p <- runif(1, .1, .9); c(p, 1 - p) },
                           numeric(2)))
  A_t <- array(vals, c(2, 2, 2, 2), dimnames = dA)
  fit <- dbn.fit(d, distribution = list(A_0 = A_0, A_t = A_t))

  u <- get_unrolled_dbn(fit, 4, extend_with_transition = TRUE)
  expect_setequal(bnlearn::nodes(u), paste0("A_", 0:4))

  # slice 1: drop A_t-2 and A_t-3 -> marginalize over both
  expect_setequal(u$A_1$parents, "A_0")
  exp1 <- apply(A_t, c("A_t", "A_t-1"), mean)
  expect_true(same_cpt(u$A_1$prob, exp1, c("A_t" = "A_1", "A_t-1" = "A_0")))
  # slice 2: drop A_t-3 only -> marginalize over it
  expect_setequal(u$A_2$parents, c("A_1", "A_0"))
  exp2 <- apply(A_t, c("A_t", "A_t-1", "A_t-2"), mean)
  expect_true(same_cpt(u$A_2$prob, exp2,
                       c("A_t" = "A_2", "A_t-1" = "A_1", "A_t-2" = "A_0")))
  # slice 3: full history
  expect_setequal(u$A_3$parents, c("A_2", "A_1", "A_0"))
  expect_true(same_cpt(u$A_3$prob, A_t,
                       c("A_t" = "A_3", "A_t-1" = "A_2", "A_t-2" = "A_1", "A_t-3" = "A_0")))
})

# ===========================================================================
# GAUSSIAN
# ===========================================================================

build_gaussian_dbn <- function() {
  d <- empty.dbn(dynamic_nodes = c("A", "B"), markov_order = 2)
  d <- add.arc.dbn(d, c("A", "t_0"), c("A", "t_1"))
  d <- add.arc.dbn(d, c("B", "t_0"), c("B", "t_1"))
  d <- add.arc.dbn(d, c("A", "t-1"), c("A", "t"))
  d <- add.arc.dbn(d, c("A", "t-2"), c("A", "t"))
  d <- add.arc.dbn(d, c("B", "t-2"), c("B", "t"))
  cpds <- list(
    A_0 = c("(Intercept)" = 0,   "Std (res)" = 1),
    A_1 = c("(Intercept)" = 0.2, "A_0" = 0.5, "Std (res)" = 0.4),
    B_0 = c("(Intercept)" = 1,   "Std (res)" = 1),
    B_1 = c("(Intercept)" = 0.3, "B_0" = 0.8, "Std (res)" = 0.5),
    A_t = c("(Intercept)" = 1,   "A_t-1" = 0.6, "A_t-2" = 0.3, "Std (res)" = 0.2),
    B_t = c("(Intercept)" = 2,   "B_t-2" = 0.7, "Std (res)" = 0.5)
  )
  dbn.fit(d, distribution = cpds)
}

test_that("gaussian mo=2, FALSE: extended G_0 + full-history transition coefficients", {
  fit <- build_gaussian_dbn()
  u <- get_unrolled_dbn(fit, 3)

  expect_setequal(bnlearn::nodes(u),
                  c("A_0", "A_1", "B_0", "B_1", "A_2", "B_2", "A_3", "B_3"))

  # extended G_0 verbatim
  expect_equal(u$A_1$coefficients[c("(Intercept)", "A_0")],
               c("(Intercept)" = 0.2, "A_0" = 0.5))
  expect_equal(unname(u$A_1$sd), 0.4)

  # full-history transition: A_t-1 -> A_1, A_t-2 -> A_0
  expect_setequal(u$A_2$parents, c("A_1", "A_0"))
  expect_equal(unname(u$A_2$coefficients[["(Intercept)"]]), 1)
  expect_equal(unname(u$A_2$coefficients[["A_1"]]), 0.6)
  expect_equal(unname(u$A_2$coefficients[["A_0"]]), 0.3)
  expect_equal(unname(u$A_2$sd), 0.2)
})

test_that("gaussian mo=2, TRUE: drop the missing lag's regressor; single dropped parent = intercept-only root", {
  fit <- build_gaussian_dbn()
  u <- get_unrolled_dbn(fit, 3, extend_with_transition = TRUE)

  expect_setequal(bnlearn::nodes(u),
                  c("A_0", "B_0", "A_1", "B_1", "A_2", "B_2", "A_3", "B_3"))

  # A_1: A_t-2 dropped -> its regressor removed; intercept, A_t-1 coef, sd unchanged
  expect_setequal(u$A_1$parents, "A_0")
  expect_setequal(names(u$A_1$coefficients), c("(Intercept)", "A_0"))
  expect_equal(unname(u$A_1$coefficients[["(Intercept)"]]), 1)
  expect_equal(unname(u$A_1$coefficients[["A_0"]]), 0.6)   # from A_t-1
  expect_equal(unname(u$A_1$sd), 0.2)

  # EDGE: B's single parent dropped -> intercept-only gaussian root, sd unchanged
  expect_length(u$B_1$parents, 0)
  expect_equal(names(u$B_1$coefficients), "(Intercept)")
  expect_equal(unname(u$B_1$coefficients[["(Intercept)"]]), 2)
  expect_equal(unname(u$B_1$sd), 0.5)
})

# ===========================================================================
# MIXED (conditional gaussian)
# ===========================================================================

build_mixed_dbn <- function() {
  lv <- c("hi", "lo")
  d <- empty.dbn(dynamic_nodes = c("D", "G", "M"), markov_order = 2)
  # extended G_0 (slice-1 nodes)
  d <- add.arc.dbn(d, c("D", "t_0"), c("D", "t_1"))
  d <- add.arc.dbn(d, c("G", "t_0"), c("G", "t_1"))
  d <- add.arc.dbn(d, c("G", "t_0"), c("M", "t_0")); d <- add.arc.dbn(d, c("D", "t_0"), c("M", "t_0"))
  d <- add.arc.dbn(d, c("G", "t_1"), c("M", "t_1")); d <- add.arc.dbn(d, c("D", "t_1"), c("M", "t_1"))
  # transition
  d <- add.arc.dbn(d, c("D", "t-1"), c("D", "t"))
  d <- add.arc.dbn(d, c("G", "t-1"), c("G", "t"))
  d <- add.arc.dbn(d, c("G", "t"),   c("M", "t"))
  d <- add.arc.dbn(d, c("D", "t"),   c("M", "t"))
  d <- add.arc.dbn(d, c("G", "t-2"), c("M", "t"))   # gaussian lag-2 -> droppable
  d <- add.arc.dbn(d, c("D", "t-2"), c("M", "t"))   # discrete lag-2 -> droppable

  D_0 <- array(c(.4, .6), 2, dimnames = list(D_0 = lv))
  D_1 <- array(c(.7, .3, .2, .8), c(2, 2), dimnames = list(D_1 = lv, D_0 = lv))
  G_0 <- c("(Intercept)" = 0, "Std (res)" = 1)
  G_1 <- c("(Intercept)" = 0.1, "G_0" = 0.5, "Std (res)" = 0.7)
  M0c <- array(c(1, 0.5, 2, 0.7), c(2, 2),
               dimnames = list(M_0 = c("(Intercept)", "G_0"), D_0 = lv))
  M_0 <- list(coef = M0c, sd = c(.2, .3))
  M1c <- array(c(1, -1, 0, 0.5), c(2, 2),
               dimnames = list(M_1 = c("(Intercept)", "G_1"), D_1 = lv))
  M_1 <- list(coef = M1c, sd = c(.15, .25))
  D_t <- array(c(.8, .2, .3, .7), c(2, 2), dimnames = list(D_t = lv, `D_t-1` = lv))
  G_t <- c("(Intercept)" = 1, "G_t-1" = 0.5, "Std (res)" = 1)
  Mtc <- array(c(1, -1, 0.2,  0, 0.5, 0.1,  3, -2, 0.3,  -1, 1.5, 0.4),
               dim = c(3, 2, 2),
               dimnames = list(M_t = c("(Intercept)", "G_t", "G_t-2"),
                               D_t = lv, `D_t-2` = lv))
  M_t <- list(coef = Mtc, sd = c(.1, .2, .3, .4))

  fit <- dbn.fit(d, distribution = list(D_0 = D_0, D_1 = D_1, G_0 = G_0, G_1 = G_1,
                                        M_0 = M_0, M_1 = M_1, D_t = D_t, G_t = G_t, M_t = M_t))
  list(fit = fit, lv = lv, M1c = M1c, M_1_sd = c(.15, .25))
}

test_that("mixed mo=2, FALSE: verbatim CLG G_0 node and full-history CLG transition node", {
  f <- build_mixed_dbn()
  u <- suppressWarnings(get_unrolled_dbn(f$fit, 3))
  fit <- f$fit; lv <- f$lv

  expect_true("bn.fit.cgnode" %in% class(u$M_1))
  expect_true("bn.fit.cgnode" %in% class(u$M_2))

  # --- verbatim G_0 CLG node M_1 (parents G_1, D_1) matches the input CLG ---
  expect_setequal(u$M_1$parents, c("G_1", "D_1"))
  for (lvl in lv) {
    got <- cg_lookup(u$M_1, list(D_1 = lvl))
    exp_coef <- f$M1c[, lvl]                       # c("(Intercept)", "G_1")
    expect_equal(unname(got$coef[c("(Intercept)", "G_1")]),
                 unname(exp_coef[c("(Intercept)", "G_1")]),
                 info = paste("M_1 coef @ D_1 =", lvl))
  }
  expect_equal(sort(unname(u$M_1$sd)), sort(f$M_1_sd))

  # --- full-history transition M_2: G_t->G_2, D_t->D_2, G_t-2->G_0, D_t-2->D_0 ---
  expect_setequal(u$M_2$parents, c("G_2", "D_2", "G_0", "D_0"))
  s_combos <- expand.grid(fit$M_t$dlevels, stringsAsFactors = FALSE)  # D_t, D_t-2
  for (i in seq_len(nrow(s_combos))) {
    src <- cg_lookup(fit$M_t, list("D_t" = s_combos[i, "D_t"],
                                   "D_t-2" = s_combos[i, "D_t-2"]))
    got <- cg_lookup(u$M_2, list(D_2 = s_combos[i, "D_t"],
                                 D_0 = s_combos[i, "D_t-2"]))
    # rename source gaussian coefficients: G_t -> G_2, G_t-2 -> G_0
    expect_equal(unname(got$coef[["(Intercept)"]]), unname(src$coef[["(Intercept)"]]))
    expect_equal(unname(got$coef[["G_2"]]),         unname(src$coef[["G_t"]]))
    expect_equal(unname(got$coef[["G_0"]]),         unname(src$coef[["G_t-2"]]))
    expect_equal(got$sd, src$sd)
  }
})

test_that("mixed mo=2, TRUE: drop a gaussian lagged parent (regressor) and a discrete lagged parent (dimension)", {
  f <- build_mixed_dbn()
  u <- suppressWarnings(get_unrolled_dbn(f$fit, 3, extend_with_transition = TRUE))
  fit <- f$fit; lv <- f$lv

  # M_1 from the transition: G_t-2 -> G_-1 (drop gaussian regressor),
  # D_t-2 -> D_-1 (drop discrete dimension, fixed to its first level)
  expect_setequal(u$M_1$parents, c("G_1", "D_1"))
  expect_true("bn.fit.cgnode" %in% class(u$M_1))          # still CG (keeps D_1)
  expect_setequal(rownames(u$M_1$coefficients), c("(Intercept)", "G_1"))  # G_t-2 row gone

  first_d2 <- fit$M_t$dlevels[["D_t-2"]][1]               # first level of the dropped parent
  for (lvl in lv) {
    # expected: source M_t at (D_t = lvl, D_t-2 = first level), G_t-2 row dropped, G_t -> G_1
    src <- cg_lookup(fit$M_t, list("D_t" = lvl, "D_t-2" = first_d2))
    got <- cg_lookup(u$M_1, list(D_1 = lvl))
    expect_equal(unname(got$coef[["(Intercept)"]]), unname(src$coef[["(Intercept)"]]),
                 info = paste("M_1 intercept @ D_1 =", lvl))
    expect_equal(unname(got$coef[["G_1"]]), unname(src$coef[["G_t"]]),
                 info = paste("M_1 G_1 slope @ D_1 =", lvl))
    expect_equal(got$sd, src$sd, info = paste("M_1 sd @ D_1 =", lvl))
  }
})

test_that("mixed edge: a CG node whose single (discrete) parent is dropped degrades to a gaussian node", {
  lv <- c("hi", "lo")
  d <- empty.dbn(dynamic_nodes = c("D", "N"), markov_order = 2)
  d <- add.arc.dbn(d, c("D", "t-1"), c("D", "t"))
  d <- add.arc.dbn(d, c("D", "t-2"), c("N", "t"))   # N: single discrete lagged parent

  D_0 <- array(c(.5, .5), 2, dimnames = list(D_0 = lv))
  N_0 <- c("(Intercept)" = 0, "Std (res)" = 1)          # gaussian root at time 0
  D_t <- array(c(.8, .2, .3, .7), c(2, 2), dimnames = list(D_t = lv, `D_t-1` = lv))
  # N_t is CLG with only a discrete parent: one (intercept) coefficient per D_t-2 level
  Nc <- array(c(5, 10), dim = c(1, 2),
              dimnames = list(N_t = "(Intercept)", `D_t-2` = lv))
  N_t <- list(coef = Nc, sd = c(0.5, 0.9))
  fit <- dbn.fit(d, distribution = list(D_0 = D_0, N_0 = N_0, D_t = D_t, N_t = N_t))
  expect_equal(dbn_type(fit), "mixed")

  u <- suppressWarnings(get_unrolled_dbn(fit, 3, extend_with_transition = TRUE))

  # N_1: its only parent D_t-2 has no history -> dropped, N_1 becomes a root
  expect_length(u$N_1$parents, 0)
  expect_true("bn.fit.gnode" %in% class(u$N_1))          # no discrete parent left -> gaussian
  # values come from fixing D_t-2 to its first level
  first_d <- fit$N_t$dlevels[["D_t-2"]][1]
  src <- cg_lookup(fit$N_t, setNames(list(first_d), "D_t-2"))
  expect_equal(unname(u$N_1$coefficients[["(Intercept)"]]), unname(src$coef[["(Intercept)"]]))
  expect_equal(unname(u$N_1$sd), src$sd)

  # once history exists (slice 2) the discrete parent is kept and N is CG again
  expect_setequal(u$N_2$parents, "D_0")
  expect_true("bn.fit.cgnode" %in% class(u$N_2))
})

# ===========================================================================
# truncation + input validation
# ===========================================================================

test_that("slices <= last G_0 slice returns only the G_0 slices (no transition)", {
  f <- build_discrete_dbn()                 # extended G_0 spans slices 0 and 1
  u <- get_unrolled_dbn(f$fit, 1)
  expect_setequal(bnlearn::nodes(u), c("A_0", "A_1", "B_0", "B_1"))
  expect_length(grep("_2$|_3$", bnlearn::nodes(u)), 0)
})

test_that("input validation", {
  f <- build_discrete_dbn()
  expect_error(get_unrolled_dbn(f$fit, "3"))
  expect_error(get_unrolled_dbn(f$fit, 0))
  expect_error(get_unrolled_dbn(c(), 3))
  expect_error(get_unrolled_dbn(f$fit, 3, extend_with_transition = "yes"))
})
