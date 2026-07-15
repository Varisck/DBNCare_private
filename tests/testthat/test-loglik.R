# Tests for logLik.dbn.fit — the plug-in (predictive) log-likelihood of data
# under a fitted DBN.
#
# The primary oracle is bnlearn: for fixed-length trajectories and markov_order
# 1, logLik(fit, data, component = "both") must equal bnlearn::logLik() on the
# fully unrolled bn.fit (get_unrolled_dbn) over the same data reshaped wide.
# All fixtures are fit from an exact `distribution =` so parameters are known
# and free of unidentifiable NAs.

library(testthat)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# Reshape long data into the wide, one-row-per-trajectory frame the unrolled
# bn.fit expects (columns var_0 .. var_slices), coercing each node column to the
# type bnlearn needs: factors (with the fitted CPT levels) for discrete nodes,
# numeric otherwise.
wide_for_unrolled <- function(data, unrolled, slices) {
  w <- build_g0_df(data, markov_order = slices + 1)
  ns <- bnlearn::nodes(unrolled)
  for (nd in ns) {
    node <- unrolled[[nd]]
    if (inherits(node, "bn.fit.dnode"))
      w[[nd]] <- factor(as.character(w[[nd]]), levels = dimnames(node$prob)[[nd]])
    else
      w[[nd]] <- as.numeric(w[[nd]])
  }
  w[, ns, drop = FALSE]
}

# ---- gaussian fixture: two independent AR(1) chains ----
build_gauss_fit <- function() {
  d <- empty.dbn(dynamic_nodes = c("A", "B"), markov_order = 1)
  d <- add.arc.dbn(d, c("A", "t-1"), c("A", "t"))
  d <- add.arc.dbn(d, c("B", "t-1"), c("B", "t"))
  cpds <- list(
    A_0 = c("(Intercept)" = 0,    "Std (res)" = 1),
    B_0 = c("(Intercept)" = 0,    "Std (res)" = 1),
    A_t = c("(Intercept)" = 0.5,  "A_t-1" = 0.6, "Std (res)" = 0.8),
    B_t = c("(Intercept)" = -0.2, "B_t-1" = 0.4, "Std (res)" = 1.2)
  )
  dbn.fit(d, distribution = cpds)
}
gen_gauss_data <- function(N = 25, Tmax = 4, seed = 11) {
  set.seed(seed)
  n <- N * (Tmax + 1)
  data.frame(Sample_id = rep(sprintf("s%02d", seq_len(N)), each = Tmax + 1),
             Time = rep(0:Tmax, N),
             A = rnorm(n), B = rnorm(n))
}

# ---- discrete fixture: two independent Markov chains ----
build_disc_fit <- function() {
  lv <- c("y", "n")
  d <- empty.dbn(dynamic_nodes = c("A", "B"), markov_order = 1)
  d <- add.arc.dbn(d, c("A", "t-1"), c("A", "t"))
  d <- add.arc.dbn(d, c("B", "t-1"), c("B", "t"))
  A_0 <- array(c(.5, .5), 2, dimnames = list(A_0 = lv))
  B_0 <- array(c(.3, .7), 2, dimnames = list(B_0 = lv))
  dA <- list(A_t = lv); dA[["A_t-1"]] <- lv
  A_t <- array(c(.8, .2, .4, .6), c(2, 2), dimnames = dA)
  dB <- list(B_t = lv); dB[["B_t-1"]] <- lv
  B_t <- array(c(.6, .4, .35, .65), c(2, 2), dimnames = dB)
  dbn.fit(d, distribution = list(A_0 = A_0, B_0 = B_0, A_t = A_t, B_t = B_t))
}
gen_disc_data <- function(N = 30, Tmax = 4, seed = 7) {
  set.seed(seed)
  n <- N * (Tmax + 1)
  data.frame(Sample_id = rep(sprintf("s%02d", seq_len(N)), each = Tmax + 1),
             Time = rep(0:Tmax, N),
             A = sample(c("y", "n"), n, replace = TRUE),
             B = sample(c("y", "n"), n, replace = TRUE),
             stringsAsFactors = FALSE)
}

# ---- mixed / CLG fixture (markov_order 1): D discrete, G gaussian, M = CLG ----
build_mixed_fit <- function() {
  lv <- c("yes", "no")
  d <- empty.dbn(dynamic_nodes = c("D", "G", "M"), markov_order = 1)
  d <- add.arc.dbn(d, c("D", "t-1"), c("D", "t"))
  d <- add.arc.dbn(d, c("G", "t-1"), c("G", "t"))
  d <- add.arc.dbn(d, c("D", "t"),   c("M", "t"))
  d <- add.arc.dbn(d, c("G", "t"),   c("M", "t"))
  d <- add.arc.dbn(d, c("D", "t_0"), c("M", "t_0"))
  d <- add.arc.dbn(d, c("G", "t_0"), c("M", "t_0"))
  D_0 <- array(c(.4, .6), 2, dimnames = list(D_0 = lv))
  G_0 <- c("(Intercept)" = 0, "Std (res)" = 1)
  M0c <- array(c(1, .5, 2, .7), c(2, 2),
               dimnames = list(M_0 = c("(Intercept)", "G_0"), D_0 = lv))
  M_0 <- list(coef = M0c, sd = c(.2, .3))
  dDt <- list(D_t = lv); dDt[["D_t-1"]] <- lv
  D_t <- array(c(.8, .2, .3, .7), c(2, 2), dimnames = dDt)
  G_t <- c("(Intercept)" = 1, "G_t-1" = 0.5, "Std (res)" = 1)
  dMt <- list(M_t = c("(Intercept)", "G_t")); dMt[["D_t"]] <- lv
  Mtc <- array(c(1, -1, 0, .5), c(2, 2), dimnames = dMt)
  M_t <- list(coef = Mtc, sd = c(.1, .2))
  dbn.fit(d, distribution = list(D_0 = D_0, G_0 = G_0, M_0 = M_0,
                                 D_t = D_t, G_t = G_t, M_t = M_t))
}
gen_mixed_data <- function(N = 30, Tmax = 4, seed = 3) {
  set.seed(seed)
  n <- N * (Tmax + 1)
  data.frame(Sample_id = rep(sprintf("s%02d", seq_len(N)), each = Tmax + 1),
             Time = rep(0:Tmax, N),
             D = sample(c("yes", "no"), n, replace = TRUE),
             G = rnorm(n), M = rnorm(n),
             stringsAsFactors = FALSE)
}

# ===========================================================================
# cross-check against bnlearn on the unrolled network
# ===========================================================================

test_that("gaussian: logLik equals bnlearn::logLik on the unrolled net", {
  fit <- build_gauss_fit(); data <- gen_gauss_data(); Tmax <- 4
  u <- get_unrolled_dbn(fit, Tmax)
  w <- wide_for_unrolled(data, u, Tmax)
  expect_equal(as.numeric(logLik(fit, data)),
               as.numeric(logLik(u, w)), tolerance = 1e-6)
})

test_that("discrete: logLik equals bnlearn::logLik on the unrolled net", {
  fit <- build_disc_fit(); data <- gen_disc_data(); Tmax <- 4
  u <- get_unrolled_dbn(fit, Tmax)
  w <- wide_for_unrolled(data, u, Tmax)
  expect_equal(as.numeric(logLik(fit, data)),
               as.numeric(logLik(u, w)), tolerance = 1e-6)
})

test_that("mixed / CLG: logLik equals bnlearn::logLik on the unrolled net", {
  fit <- build_mixed_fit(); data <- gen_mixed_data(); Tmax <- 4
  u <- suppressWarnings(get_unrolled_dbn(fit, Tmax))
  w <- wide_for_unrolled(data, u, Tmax)
  expect_equal(as.numeric(logLik(fit, data)),
               as.numeric(logLik(u, w)), tolerance = 1e-6)
})

# ===========================================================================
# component decomposition, attributes, AIC/BIC, by.sample, node filter, na.rm
# ===========================================================================

test_that("component 'both' decomposes into 'initial' + 'transition'", {
  fit <- build_gauss_fit(); data <- gen_gauss_data()
  both  <- as.numeric(logLik(fit, data, component = "both"))
  init  <- as.numeric(logLik(fit, data, component = "initial"))
  trans <- as.numeric(logLik(fit, data, component = "transition"))
  expect_equal(both, init + trans, tolerance = 1e-9)
})

test_that("logLik carries nobs/df so AIC and BIC follow their formulas", {
  fit <- build_disc_fit(); data <- gen_disc_data()
  ll <- logLik(fit, data)
  expect_s3_class(ll, "logLik")
  k    <- attr(ll, "df")
  nobs <- attr(ll, "nobs")
  expect_true(is.finite(k) && k > 0)
  expect_true(is.finite(nobs) && nobs > 0)
  expect_equal(AIC(ll), -2 * as.numeric(ll) + 2 * k, tolerance = 1e-9)
  expect_equal(BIC(ll), -2 * as.numeric(ll) + log(nobs) * k, tolerance = 1e-9)
})

test_that("df equals the number of free parameters (2-level chain fixture)", {
  fit <- build_disc_fit(); data <- gen_disc_data()
  # A_0,B_0: 1 free param each (2 levels, no parents) = 2
  # A_t,B_t: (2-1)*2 configs = 2 each                 = 4
  expect_equal(attr(logLik(fit, data), "df"), 6)
})

test_that("by.sample sums to the total and is one value per trajectory", {
  fit <- build_gauss_fit(); data <- gen_gauss_data(N = 25)
  bs <- logLik(fit, data, by.sample = TRUE)
  expect_length(bs, length(unique(data$Sample_id)))
  expect_equal(sum(bs), as.numeric(logLik(fit, data)), tolerance = 1e-9)
})

test_that("the node filter partitions the total across variables", {
  fit <- build_gauss_fit(); data <- gen_gauss_data()
  llA <- as.numeric(logLik(fit, data, nodes = "A"))
  llB <- as.numeric(logLik(fit, data, nodes = "B"))
  expect_equal(llA + llB, as.numeric(logLik(fit, data)), tolerance = 1e-9)
})

test_that("held-out predictive log-likelihood is finite on unseen data", {
  fit <- build_gauss_fit()
  test <- gen_gauss_data(N = 10, seed = 99)   # different sample than any train set
  expect_true(is.finite(as.numeric(logLik(fit, test, component = "transition"))))
})

test_that("unscorable observations give NA, dropped with na.rm = TRUE", {
  fit <- build_disc_fit()
  data <- gen_disc_data()
  data$A[data$Time == 2][1] <- "z"            # level never seen in training
  expect_true(is.na(as.numeric(logLik(fit, data))))
  expect_true(is.finite(as.numeric(logLik(fit, data, na.rm = TRUE))))
})

test_that("input validation", {
  fit <- build_gauss_fit()
  expect_error(logLik.dbn.fit(list(), data.frame()), "dbn.fit")
  expect_error(logLik(fit, "not a data.frame"), "data.frame")
  expect_error(logLik(fit, gen_gauss_data(), nodes = "Z"), "unrecognized")
})
