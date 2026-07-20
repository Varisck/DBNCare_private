# dbn.sampling() must replicate dbn.sampling.R() draw for draw: the Rcpp
# cores consume the R RNG stream in the same order as the pure-R
# implementation, so with the same seed the two functions must return the
# same dataset.

# named CPD vector in the format expected by dbn.fit(distribution = ...)
cpd_for <- function(dbn, node, intercept, coefs, sd) {
  parents <- get_parent_set(dbn, node)
  stopifnot(setequal(names(coefs), parents))
  values <- c(intercept, as.numeric(coefs[parents]), sd)
  names(values) <- c("(Intercept)", parents, "Std (res)")
  values
}

# binary CPT array in the format expected by dbn.fit(distribution = ...):
# P(node = "yes" | parents) = base_p + sum of bumps[parent] over parents at "yes"
cpt_for <- function(dbn, node, base_p, bumps = c(), levels = c("no", "yes")) {
  parents <- get_parent_set(dbn, node)
  stopifnot(setequal(names(bumps), parents))
  if (length(parents) == 0) {
    probs <- c(1 - base_p, base_p)
  } else {
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

# markov order 1: G_0 collider + persistence and collider in G_t
structure_mo1 <- function() {
  dbn <- empty.dbn(dynamic_nodes = c("A", "B", "C"), markov_order = 1)
  dbn <- add.arc.dbn(dbn, from = c("A", "t_0"), to = c("C", "t_0"))
  dbn <- add.arc.dbn(dbn, from = c("B", "t_0"), to = c("C", "t_0"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t-1"), to = c("A", "t"))
  dbn <- add.arc.dbn(dbn, from = c("B", "t-1"), to = c("B", "t"))
  dbn <- add.arc.dbn(dbn, from = c("C", "t-1"), to = c("C", "t"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t"), to = c("C", "t"))
  dbn <- add.arc.dbn(dbn, from = c("B", "t"), to = c("C", "t"))
  dbn
}

# markov order 2: B_t depends on A_t-2, which is not available at t = 1; this
# exercises the unavailable-parent paths (skipped regressor for gaussian
# networks, CPT rows pooled over the missing parent for discrete ones)
structure_mo2 <- function() {
  dbn <- empty.dbn(dynamic_nodes = c("A", "B"), markov_order = 2)
  dbn <- add.arc.dbn(dbn, from = c("A", "t_0"), to = c("B", "t_0"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t-1"), to = c("A", "t"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t-2"), to = c("B", "t"))
  dbn <- add.arc.dbn(dbn, from = c("B", "t-1"), to = c("B", "t"))
  dbn
}

test_that("dbn.sampling matches dbn.sampling on a gaussian network", {
  truth <- structure_mo1()
  fit <- dbn.fit(DBN = truth, distribution = list(
    A_0 = cpd_for(truth, "A_0", 0.0, c(), 1.0),
    B_0 = cpd_for(truth, "B_0", 0.0, c(), 1.0),
    C_0 = cpd_for(truth, "C_0", 1.0, c(A_0 = 0.8, B_0 = -0.7), 0.5),
    A_t = cpd_for(truth, "A_t", 0.0, c(`A_t-1` = 0.6), 0.4),
    B_t = cpd_for(truth, "B_t", 0.2, c(`B_t-1` = 0.5), 0.4),
    C_t = cpd_for(truth, "C_t", -0.3, c(`C_t-1` = 0.4, A_t = 0.7, B_t = -0.6), 0.3)
  ))

  set.seed(81)
  reference <- dbn.sampling.R(fit, 60, 5)
  set.seed(81)
  sampled <- dbn.sampling(fit, 60, 5)
  expect_equal(sampled, reference, tolerance = 1e-12)
})

test_that("dbn.sampling matches dbn.sampling on a discrete network", {
  truth <- structure_mo1()
  fit <- dbn.fit(DBN = truth, distribution = list(
    A_0 = cpt_for(truth, "A_0", 0.45),
    B_0 = cpt_for(truth, "B_0", 0.55),
    C_0 = cpt_for(truth, "C_0", 0.10, c(A_0 = 0.35, B_0 = 0.45)),
    A_t = cpt_for(truth, "A_t", 0.12, c(`A_t-1` = 0.72)),
    B_t = cpt_for(truth, "B_t", 0.15, c(`B_t-1` = 0.65)),
    C_t = cpt_for(truth, "C_t", 0.05, c(`C_t-1` = 0.30, A_t = 0.30, B_t = 0.30))
  ))

  set.seed(82)
  reference <- dbn.sampling.R(fit, 80, 5)
  set.seed(82)
  sampled <- dbn.sampling(fit, 80, 5)
  expect_identical(sampled, reference)
})

test_that("dbn.sampling matches dbn.sampling with markov order 2 (gaussian)", {
  truth <- structure_mo2()
  fit <- dbn.fit(DBN = truth, distribution = list(
    A_0 = cpd_for(truth, "A_0", 0.0, c(), 1.0),
    B_0 = cpd_for(truth, "B_0", 0.5, c(A_0 = 0.9), 0.7),
    A_t = cpd_for(truth, "A_t", 0.1, c(`A_t-1` = 0.6), 0.4),
    B_t = cpd_for(truth, "B_t", -0.2, c(`A_t-2` = 0.8, `B_t-1` = 0.5), 0.3)
  ))

  set.seed(83)
  reference <- dbn.sampling.R(fit, 50, 6)
  set.seed(83)
  sampled <- dbn.sampling(fit, 50, 6)
  expect_equal(sampled, reference, tolerance = 1e-12)
})

test_that("dbn.sampling matches dbn.sampling with markov order 2 (discrete)", {
  truth <- structure_mo2()
  fit <- dbn.fit(DBN = truth, distribution = list(
    A_0 = cpt_for(truth, "A_0", 0.40),
    B_0 = cpt_for(truth, "B_0", 0.20, c(A_0 = 0.50)),
    A_t = cpt_for(truth, "A_t", 0.15, c(`A_t-1` = 0.70)),
    B_t = cpt_for(truth, "B_t", 0.10, c(`A_t-2` = 0.35, `B_t-1` = 0.45))
  ))

  set.seed(84)
  reference <- dbn.sampling.R(fit, 50, 6)
  set.seed(84)
  sampled <- dbn.sampling(fit, 50, 6)
  expect_identical(sampled, reference)
})

test_that("dbn.sampling handles CPTs whose levels are ordered differently", {
  # the levels of A are listed as no/yes in the G_0 CPTs but as yes/no in the
  # transition CPTs; values must still be matched by label, like filter_cpt()
  dbn <- empty.dbn(dynamic_nodes = c("A", "B"), markov_order = 1)
  dbn <- add.arc.dbn(dbn, from = c("A", "t_0"), to = c("B", "t_0"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t-1"), to = c("A", "t"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t"), to = c("B", "t"))
  dbn <- add.arc.dbn(dbn, from = c("B", "t-1"), to = c("B", "t"))

  A_lv_0 <- c("no", "yes")
  A_lv_t <- c("yes", "no")
  B_lv <- c("hi", "lo")

  dims_A_t <- list(A_t = A_lv_t)
  dims_A_t[["A_t-1"]] <- A_lv_t
  dims_B_t <- list(B_t = B_lv, A_t = A_lv_t)
  dims_B_t[["B_t-1"]] <- B_lv

  fit <- dbn.fit(DBN = dbn, distribution = list(
    A_0 = array(c(0.4, 0.6), dim = 2, dimnames = list(A_0 = A_lv_0)),
    B_0 = array(c(0.3, 0.7, 0.8, 0.2), dim = c(2, 2),
                dimnames = list(B_0 = B_lv, A_0 = A_lv_0)),
    A_t = array(c(0.75, 0.25, 0.2, 0.8), dim = c(2, 2), dimnames = dims_A_t),
    B_t = array(c(0.15, 0.85, 0.7, 0.3, 0.5, 0.5, 0.05, 0.95), dim = c(2, 2, 2),
                dimnames = dims_B_t)
  ))

  set.seed(85)
  reference <- dbn.sampling.R(fit, 60, 4)
  set.seed(85)
  sampled <- dbn.sampling(fit, 60, 4)
  expect_identical(sampled, reference)
})

test_that("dbn.sampling validates its inputs like dbn.sampling", {
  truth <- structure_mo1()
  fit <- dbn.fit(DBN = truth, distribution = list(
    A_0 = cpt_for(truth, "A_0", 0.45),
    B_0 = cpt_for(truth, "B_0", 0.55),
    C_0 = cpt_for(truth, "C_0", 0.10, c(A_0 = 0.35, B_0 = 0.45)),
    A_t = cpt_for(truth, "A_t", 0.12, c(`A_t-1` = 0.72)),
    B_t = cpt_for(truth, "B_t", 0.15, c(`B_t-1` = 0.65)),
    C_t = cpt_for(truth, "C_t", 0.05, c(`C_t-1` = 0.30, A_t = 0.30, B_t = 0.30))
  ))

  expect_error(dbn.sampling(fit, "5", 4), "N_samples must be an integer!")
  expect_error(dbn.sampling(fit, 5, "4"), "Time must be an integer!")
  expect_error(dbn.sampling(fit, 0, 4), "N_samples must be greater than 0!")
  expect_error(dbn.sampling(fit, 5, 0), "Time must be greater than 0!")
  expect_error(dbn.sampling(c("dbn"), 5, 4), "fitted_DBN must be a dbn.fit object")
})

# Extended G_0: when markov_order > 1 the prior network can span several initial
# slices (var_0, var_1, ...) linked by prior arcs running ALONG the markov
# dimension (e.g. A_0 -> A_1). dbn.sampling must sample each extended slice var_i
# and place it at Time == i (driven by its own G_0 regressor), only then handing
# over to the transition network. dbn.sampling.R cannot sample an extended G_0,
# so these tests check dbn.sampling's output directly rather than comparing.

test_that("dbn.sampling places extended-G_0 slices along the markov dimension", {
  # markov order 3, single variable, a prior chain A_0 -> A_1 -> A_2 inside G_0
  dbn <- empty.dbn(dynamic_nodes = c("A"), markov_order = 3)
  dbn <- add.arc.dbn(dbn, from = c("A", "t_0"), to = c("A", "t_1"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t_1"), to = c("A", "t_2"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t-1"), to = c("A", "t"))

  # clear G_0 regressors with (essentially) no residual noise, so the extended
  # slices are deterministic functions of A_0: A_1 = 4*A_0, A_2 = 3*A_1
  fit <- dbn.fit(DBN = dbn, distribution = list(
    A_0 = cpd_for(dbn, "A_0", 0.0, c(),           1.0),
    A_1 = cpd_for(dbn, "A_1", 0.0, c(A_0 = 4.0),  1e-6),
    A_2 = cpd_for(dbn, "A_2", 0.0, c(A_1 = 3.0),  1e-6),
    A_t = cpd_for(dbn, "A_t", 0.0, c(`A_t-1` = 0.6), 0.4)
  ))
  expect_equal(markov_order(fit), 3)

  set.seed(101)
  sampled <- dbn.sampling(fit, n_samples = 500, max_time = 5)

  # the extended slices must map to the single base column A (no duplicate
  # column) and the frame keeps the standard (max_time + 1) rows per sample
  expect_setequal(names(sampled), c("Time", "Sample_id", "A"))
  expect_equal(ncol(sampled), 3L)
  expect_equal(nrow(sampled), 500L * (5L + 1L))

  a0 <- sampled$A[sampled$Time == 0]
  a1 <- sampled$A[sampled$Time == 1]
  a2 <- sampled$A[sampled$Time == 2]
  a3 <- sampled$A[sampled$Time == 3]

  # G_0 regressors are honoured slice by slice, A_i placed at Time == i
  expect_equal(a1, 4 * a0, tolerance = 1e-3)
  expect_equal(a2, 3 * a1, tolerance = 1e-3)

  # from Time == 3 the transition network takes over (A_t = 0.6 * A_t-1)
  expect_equal(unname(coef(lm(a3 ~ a2))[2]), 0.6, tolerance = 0.1)
})

test_that("dbn.sampling handles cross-variable extended-G_0 arcs", {
  # two variables, a within-variable (A_0 -> A_1) and a cross-variable
  # (A_0 -> B_1) prior arc along the markov dimension
  dbn <- empty.dbn(dynamic_nodes = c("A", "B"), markov_order = 2)
  dbn <- add.arc.dbn(dbn, from = c("A", "t_0"), to = c("A", "t_1"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t_0"), to = c("B", "t_1"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t-1"), to = c("A", "t"))
  dbn <- add.arc.dbn(dbn, from = c("B", "t-1"), to = c("B", "t"))

  fit <- dbn.fit(DBN = dbn, distribution = list(
    A_0 = cpd_for(dbn, "A_0", 0.0, c(),            1.0),
    B_0 = cpd_for(dbn, "B_0", 0.0, c(),            1.0),
    A_1 = cpd_for(dbn, "A_1", 0.0, c(A_0 =  2.0),  1e-6),
    B_1 = cpd_for(dbn, "B_1", 1.0, c(A_0 = -3.0),  1e-6),
    A_t = cpd_for(dbn, "A_t", 0.0, c(`A_t-1` = 0.5), 0.3),
    B_t = cpd_for(dbn, "B_t", 0.0, c(`B_t-1` = 0.5), 0.3)
  ))

  set.seed(202)
  sampled <- dbn.sampling(fit, n_samples = 400, max_time = 4)

  expect_setequal(names(sampled), c("Time", "Sample_id", "A", "B"))
  expect_equal(ncol(sampled), 4L)

  a0 <- sampled$A[sampled$Time == 0]
  a1 <- sampled$A[sampled$Time == 1]
  b1 <- sampled$B[sampled$Time == 1]

  # A_1 reads A at Time 0; B_1 reads A (not B) at Time 0 -> the cross-variable,
  # cross-slice parent reference resolves to the right column and slice
  expect_equal(a1, 2 * a0, tolerance = 1e-3)
  expect_equal(b1, 1 - 3 * a0, tolerance = 1e-3)
})
