# Tests for parameter learning in the mixed (conditional gaussian) case.

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

build_two_discrete_parents_dbn <- function() {
  d <- empty.dbn(dynamic_nodes = c("D1", "D2", "N"), markov_order = 1)
  d <- add.arc.dbn(d, from = c("D1", "t-1"), to = c("D1", "t"))
  d <- add.arc.dbn(d, from = c("D2", "t-1"), to = c("D2", "t"))
  d <- add.arc.dbn(d, from = c("D1", "t"),   to = c("N",  "t"))
  d <- add.arc.dbn(d, from = c("D2", "t"),   to = c("N",  "t"))
  d <- add.arc.dbn(d, from = c("D1", "t_0"), to = c("N",  "t_0"))
  d <- add.arc.dbn(d, from = c("D2", "t_0"), to = c("N",  "t_0"))
  d
}

build_two_discrete_parents_dist <- function() {
  lv1 <- c("yes", "no")
  lv2 <- c("a", "b", "c")

  D1_0 <- array(c(0.5, 0.5), dim = 2, dimnames = list(D1_0 = lv1))
  D2_0 <- array(c(0.3, 0.4, 0.3), dim = 3, dimnames = list(D2_0 = lv2))

  # N_0: CLG, no continuous parent, discrete parents D1_0 and D2_0
  # coef first dim = intercept only; remaining dims = D1_0 x D2_0 (2 x 3 = 6 combos)
  N_0_coef <- array(
    c(1, 2, 3, 4, 5, 6),   # one intercept per (D1_0, D2_0) combo
    dim = c(1, 2, 3),
    dimnames = list(N_0 = c("(Intercept)"), D1_0 = lv1, D2_0 = lv2)
  )
  N_0 <- list(coef = N_0_coef, sd = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))

  dims_D1_t <- list(D1_t = lv1); dims_D1_t[["D1_t-1"]] <- lv1
  D1_t <- array(c(0.8, 0.2, 0.3, 0.7), dim = c(2, 2), dimnames = dims_D1_t)

  dims_D2_t <- list(D2_t = lv2); dims_D2_t[["D2_t-1"]] <- lv2
  D2_t <- array(
    c(0.7, 0.2, 0.1,
      0.1, 0.8, 0.1,
      0.1, 0.2, 0.7),
    dim = c(3, 3), dimnames = dims_D2_t
  )

  dims_N_t <- list(N_t = c("(Intercept)"))
  dims_N_t[["D1_t"]] <- lv1; dims_N_t[["D2_t"]] <- lv2
  N_t_coef <- array(
    c(10, 20, 30, 40, 50, 60),
    dim = c(1, 2, 3),
    dimnames = dims_N_t
  )
  N_t <- list(coef = N_t_coef, sd = c(1, 2, 3, 4, 5, 6))

  list(D1_0 = D1_0, D2_0 = D2_0, N_0 = N_0,
       D1_t = D1_t, D2_t = D2_t, N_t = N_t)
}

test_that("CLG node with two discrete parents and no continuous parent fits correctly", {
  dbn  <- build_two_discrete_parents_dbn()
  dist <- build_two_discrete_parents_dist()

  fitted <- dbn.fit(DBN = dbn, distribution = dist)

  expect_equal(class(fitted), "dbn.fit")
  expect_equal(dbn_type(fitted), "mixed")

  for (variable in c("N_0", "N_t")) {
    expect_s3_class(fitted[[variable]], "dbn.fit.cgnode")

    # only the intercept row — no continuous-parent regressors
    expect_equal(dimnames(fitted[[variable]]$coefficients)[[1]],
                 c("(Intercept)"),
                 info = paste("regressor labels mismatch for", variable))

    # 2 * 3 = 6 sd values, one per discrete-parent combination
    expect_length(fitted[[variable]]$sd, 6)
    expect_equal(fitted[[variable]]$sd, dist[[variable]]$sd,
                 info = paste("sd mismatch for", variable))

    # dlevels carries the discrete parents' level vectors
    expect_equal(length(fitted[[variable]]$dlevels), 2,
                 info = paste("expected 2 discrete parents in dlevels for", variable))
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

# ----- 4. learn mixed parameters from DATA (no sampling) ---------------------
# Structure (markov_order = 1)
#   D1_t-1 -> D1_t      D2_t-1 -> D2_t      G_t-1 -> G_t       (self-loops)
#   G_t -> M_t          D1_t -> M_t         D2_t -> M_t        (M_t is CLG)
#   G_0 -> M_0          D1_0 -> M_0         D2_0 -> M_0        (M_0 is CLG)
#
# Node types
#   D1, D2 : discrete
#   G      : gaussian
#   M      : conditional gaussian, 1 continuous + 2 discrete parents
#              M | (G, D1, D2) ~ N(b0[D1,D2] + b1[D1,D2] * G,  sd[D1,D2])

build_clg_data_dbn <- function() {
  d <- empty.dbn(dynamic_nodes = c("D1", "D2", "G", "M"), markov_order = 1)

  # transition self-loops
  d <- add.arc.dbn(d, from = c("D1", "t-1"), to = c("D1", "t"))
  d <- add.arc.dbn(d, from = c("D2", "t-1"), to = c("D2", "t"))
  d <- add.arc.dbn(d, from = c("G",  "t-1"), to = c("G",  "t"))

  # CLG node M_t: continuous parent G_t, discrete parents D1_t, D2_t
  d <- add.arc.dbn(d, from = c("G",  "t"),   to = c("M",  "t"))
  d <- add.arc.dbn(d, from = c("D1", "t"),   to = c("M",  "t"))
  d <- add.arc.dbn(d, from = c("D2", "t"),   to = c("M",  "t"))

  # time-0 slice: M_0 has the same parents
  d <- add.arc.dbn(d, from = c("G",  "t_0"), to = c("M",  "t_0"))
  d <- add.arc.dbn(d, from = c("D1", "t_0"), to = c("M",  "t_0"))
  d <- add.arc.dbn(d, from = c("D2", "t_0"), to = c("M",  "t_0"))

  d
}

# CLG generating parameters, keyed by paste(D1, D2, sep = ".")
clg_lv1 <- c("low", "high")
clg_lv2 <- c("x", "y", "z")
clg_b0 <- c("low.x" =  1.0, "low.y" = -1.0, "low.z" =  2.0,
            "high.x" =  0.5, "high.y" =  3.0, "high.z" = -2.0)
clg_b1 <- c("low.x" =  0.5, "low.y" =  1.2, "low.z" = -0.5,
            "high.x" =  2.0, "high.y" =  0.8, "high.z" = -1.0)
clg_sd <- c("low.x" =  0.30, "low.y" =  0.40, "low.z" =  0.20,
            "high.x" =  0.50, "high.y" =  0.30, "high.z" =  0.25)


clg_d1_init  <- c(low = 0.5, high = 0.5)
clg_d1_trans <- rbind(low  = c(low = 0.75, high = 0.25),
                      high = c(low = 0.35, high = 0.65))

clg_d2_init  <- c(x = 1 / 3, y = 1 / 3, z = 1 / 3)
clg_d2_trans <- rbind(x = c(x = 0.6, y = 0.3, z = 0.1),
                      y = c(x = 0.2, y = 0.6, z = 0.2),
                      z = c(x = 0.1, y = 0.3, z = 0.6))

# One Markov step for a whole cohort at once: prev / return value are length-N
# state vectors; sampling is batched per from-state so it stays vectorised.
markov_step <- function(prev, levels, transmat) {
  nxt <- character(length(prev))
  for (s in levels) {
    idx <- which(prev == s)
    if (length(idx) > 0)
      nxt[idx] <- sample(levels, length(idx), replace = TRUE,
                         prob = transmat[s, levels])
  }
  nxt
}

# Generate a long-format data.frame straight from the generating process
generate_clg_data <- function(N, Tmax, seed = 20240618) {
  set.seed(seed)
  n_time <- Tmax + 1
  total  <- N * n_time

  Sample_id <- rep(sprintf("s%04d", seq_len(N)), each = n_time)
  Time      <- rep(0:Tmax, times = N)
  # rows with Time == k sit at positions k+1, k+1+n_time, ... (one per sample)
  pos_at <- function(k) seq(k + 1, total, by = n_time)

  D1 <- character(total)
  D2 <- character(total)
  G  <- numeric(total)

  # initial slice (Time 0)
  d1_prev <- sample(clg_lv1, N, replace = TRUE, prob = clg_d1_init[clg_lv1])
  d2_prev <- sample(clg_lv2, N, replace = TRUE, prob = clg_d2_init[clg_lv2])
  g_prev  <- rnorm(N, mean = 0, sd = 2)
  D1[pos_at(0)] <- d1_prev
  D2[pos_at(0)] <- d2_prev
  G[pos_at(0)]  <- g_prev

  for (k in seq_len(Tmax)) {
    d1_prev <- markov_step(d1_prev, clg_lv1, clg_d1_trans)
    d2_prev <- markov_step(d2_prev, clg_lv2, clg_d2_trans)
    g_prev  <- 0.6 * g_prev + rnorm(N, mean = 0, sd = 1.5)
    D1[pos_at(k)] <- d1_prev
    D2[pos_at(k)] <- d2_prev
    G[pos_at(k)]  <- g_prev
  }

  key <- paste(D1, D2, sep = ".")
  M <- unname(clg_b0[key] + clg_b1[key] * G +
                rnorm(total, mean = 0, sd = clg_sd[key]))

  data.frame(Sample_id = Sample_id, Time = Time,
             D1 = D1, D2 = D2, G = G, M = M,
             stringsAsFactors = FALSE, check.names = FALSE)
}

clg_dbn  <- build_clg_data_dbn()
clg_data <- generate_clg_data(N = 1500, Tmax = 6)
clg_fit  <- dbn.fit(DBN = clg_dbn, data = clg_data)

# absolute-tolerance closeness check with an informative failure message
expect_close <- function(object, expected, tol, label = "") {
  diff <- abs(object - expected)
  testthat::expect_true(
    isTRUE(diff < tol),
    info = sprintf("%s: |%.4f - %.4f| = %.4f (tol %.3f)",
                   label, object, expected, diff, tol)
  )
}

# Check a fitted CLG node against the generating look-up tables
check_clg_node <- function(node, g_name, d1_name, d2_name,
                           tol_coef, tol_sd) {
  combos <- expand.grid(lapply(node$dlevels, as.character),
                        stringsAsFactors = FALSE)
  coefs  <- node$coefficients

  testthat::expect_equal(dimnames(coefs)[[1]], c("(Intercept)", g_name))
  testthat::expect_length(node$sd, nrow(combos))

  for (i in seq_len(nrow(combos))) {
    key <- paste(combos[i, d1_name], combos[i, d2_name], sep = ".")
    expect_close(coefs["(Intercept)", i], clg_b0[[key]], tol_coef,
                 label = paste0(node$node, " intercept @ ", key))
    expect_close(coefs[g_name, i],        clg_b1[[key]], tol_coef,
                 label = paste0(node$node, " slope @ ", key))
    expect_close(node$sd[i],              clg_sd[[key]], tol_sd,
                 label = paste0(node$node, " sd @ ", key))
  }
}

test_that("mixed data path produces a mixed dbn.fit with the right node classes", {
  expect_equal(class(clg_fit), "dbn.fit")
  expect_equal(dbn_type(clg_fit), "mixed")
  expect_setequal(names(clg_fit),
                  c("D1_0", "D2_0", "G_0", "M_0",
                    "D1_t", "D2_t", "G_t", "M_t"))

  for (v in c("D1_0", "D2_0", "D1_t", "D2_t"))
    expect_s3_class(clg_fit[[v]], "dbn.fit.dnode")
  for (v in c("G_0", "G_t"))
    expect_s3_class(clg_fit[[v]], "dbn.fit.gnode")
  for (v in c("M_0", "M_t"))
    expect_s3_class(clg_fit[[v]], "dbn.fit.cgnode")
})

test_that("CLG coefficients and sd learned from data match the generating process", {
  # M_t is fit on the transition slice (lots of data) -> tight tolerance.
  check_clg_node(clg_fit[["M_t"]], g_name = "G_t",
                 d1_name = "D1_t", d2_name = "D2_t",
                 tol_coef = 0.1, tol_sd = 0.05)
  # M_0 is fit on the time-0 slice only (N rows) -> looser tolerance.
  check_clg_node(clg_fit[["M_0"]], g_name = "G_0",
                 d1_name = "D1_0", d2_name = "D2_0",
                 tol_coef = 0.15, tol_sd = 0.1)
})

test_that("gaussian node learned from data matches the AR(1) generating process", {
  # G_t : intercept 0, lag coefficient 0.6, innovation sd 1.5
  regs_t <- unname(clg_fit[["G_t"]]$regs)
  expect_close(regs_t[1], 0.0, 0.1,           label = "G_t intercept")
  expect_close(regs_t[2], 0.6, 0.1,           label = "G_t lag coef")
  expect_close(unname(clg_fit[["G_t"]]$std), 1.5, 0.1, label = "G_t sd")

  # G_0 : intercept only -> intercept 0, marginal sd 2
  regs_0 <- unname(clg_fit[["G_0"]]$regs)
  expect_close(regs_0[1], 0.0, 0.15,          label = "G_0 intercept")
  expect_close(unname(clg_fit[["G_0"]]$std), 2.0, 0.15, label = "G_0 sd")
})

# Verify a learned transition CPT P(var_t | var_t-1) against the generating
# transition matrix. The CPT array has dims [var_t] x [var_t-1] (both sorted
# levels), so prob[to, from] = P(var_t = to | var_t-1 = from).
check_transition_cpt <- function(prob, levels, transmat, tol = 0.05) {
  var <- names(dimnames(prob))[1]
  for (from in levels) for (to in levels)
    expect_close(prob[to, from], transmat[from, to], tol,
                 label = sprintf("P(%s=%s | _t-1=%s)", var, to, from))
}

test_that("discrete nodes learned from data recover the Markov chain", {
  # time-0 slice ~ the (uniform) initial distributions
  expect_true(all(abs(clg_fit[["D1_0"]]$prob - 1 / 2) < 0.05))
  expect_true(all(abs(clg_fit[["D2_0"]]$prob - 1 / 3) < 0.05))

  # transition slice ~ the generating transition matrices (the self-loops)
  check_transition_cpt(clg_fit[["D1_t"]]$prob, clg_lv1, clg_d1_trans)
  check_transition_cpt(clg_fit[["D2_t"]]$prob, clg_lv2, clg_d2_trans)

  # every conditional distribution (column) is a proper distribution
  expect_true(all(abs(colSums(clg_fit[["D1_t"]]$prob) - 1) < 1e-8))
  expect_true(all(abs(colSums(clg_fit[["D2_t"]]$prob) - 1) < 1e-8))
})
