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


build_ground_truth_CPTs <- function() {
  lv <- c("yes", "no")

  # B_0: no parents
  B_0.prob <- array(c(0.4, 0.6),
                    dim = 2,
                    dimnames = list(B_0 = lv))
  # C_0: no parents
  C_0.prob <- array(c(0.5, 0.5),
                    dim = 2,
                    dimnames = list(C_0 = lv))
  # A_0: parents = [B_0, C_0]  (arc order: B first, then C)
  # 4 conditioning combos, each pair sums to 1
  A_0.prob <- array(
    c(
      0.7, 0.3,   # A_0 | B_0="yes", C_0="yes"
      0.4, 0.6,   # A_0 | B_0="no",  C_0="yes"
      0.6, 0.4,   # A_0 | B_0="yes", C_0="no"
      0.2, 0.8    # A_0 | B_0="no",  C_0="no"
    ),
    dim = c(2, 2, 2),
    dimnames = list(A_0 = lv, B_0 = lv, C_0 = lv)
  )

  # --- transition slice ---

  # B_t: parents = [B_t-1]
  dims_B_t <- list(B_t = lv)
  dims_B_t[["B_t-1"]] <- lv
  B_t.prob <- array(
    c(
      0.8, 0.2,   # B_t | B_t-1="yes"
      0.3, 0.7    # B_t | B_t-1="no"
    ),
    dim = c(2, 2),
    dimnames = dims_B_t
  )

  # C_t: parents = [C_t-1, B_t]  (arc order: C_t-1 first, then B_t)
  dims_C_t <- list(C_t = lv)
  dims_C_t[["C_t-1"]] <- lv
  dims_C_t[["B_t"]] <- lv
  C_t.prob <- array(
    c(
      0.9, 0.1,   # C_t | C_t-1="yes", B_t="yes"
      0.4, 0.6,   # C_t | C_t-1="no",  B_t="yes"
      0.6, 0.4,   # C_t | C_t-1="yes", B_t="no"
      0.1, 0.9    # C_t | C_t-1="no",  B_t="no"
    ),
    dim = c(2, 2, 2),
    dimnames = dims_C_t
  )

  # A_t: parents = [A_t-1, B_t]  (arc order: A_t-1 first, then B_t)
  dims_A_t <- list(A_t = lv)
  dims_A_t[["A_t-1"]] <- lv
  dims_A_t[["B_t"]] <- lv
  A_t.prob <- array(
    c(
      0.85, 0.15,  # A_t | A_t-1="yes", B_t="yes"
      0.3,  0.7,   # A_t | A_t-1="no",  B_t="yes"
      0.7,  0.3,   # A_t | A_t-1="yes", B_t="no"
      0.2,  0.8    # A_t | A_t-1="no",  B_t="no"
    ),
    dim = c(2, 2, 2),
    dimnames = dims_A_t
  )
  list(
    B_0 = B_0.prob,
    C_0 = C_0.prob,
    A_0 = A_0.prob,
    B_t = B_t.prob,
    C_t = C_t.prob,
    A_t = A_t.prob
  )
}



dbn_truth <- build_ground_truth_dbn()
cpts_truth <- build_ground_truth_CPTs()
fitted_truth <- dbn.fit(DBN = dbn_truth, CPTs = cpts_truth)


# sanity: the auto-detected dbn type must be gaussian
test_that("ground-truth fitted DBN is gaussian", {
  expect_equal(dbn_type(fitted_truth), "discrete")
})

# sample once and reuse across tests (keeps the suite fast).
# seed is set right before sampling so the dataset is pinned regardless of
# RNG state inherited from earlier test files.
set.seed(123)
sampled_data <- dbn.sampling(fitted_truth, n_samples = 1e4, max_time = 4)

# sanity: sampled data is numeric (so dataset_type should report gaussian) and
# has the expected columns
test_that("sampled dataset is character", {
  expect_equal(dataset_type(sampled_data), "discrete")
  expect_true(all(c("Sample_id", "Time", "A", "B", "C") %in% names(sampled_data)))
  expect_true(all(sapply(sampled_data[, c("A", "B", "C")], 
                          \(x){ is.character(x) || is.factor(x)})))
})


# score-based algorithms (hc, tabu) are the natural fit for continuous data.
# Passing score = NULL lets bnlearn pick its default gaussian score (bic-g),
# which is what works on numeric columns.
test_that("hc recovers the ground-truth discrete DBN structure", {
  set.seed(1)
  learned <- dbn.learn.structure(sampled_data, algorithm = "hc",
                                       score = NULL)
  expect_s3_class(learned, "dbn")
  tn_l = bnlearn::model2network(modelstring(learned)$g_t)
  tn_gt = bnlearn::model2network(modelstring(dbn_truth)$g_t)
  expect_equal(all.equal(tn_l, tn_gt), TRUE)

  g_0_l = bnlearn::model2network(modelstring(learned)$g_0)
  g_0_gt = bnlearn::model2network(modelstring(dbn_truth)$g_0)
  expect_equal(all.equal(g_0_l, g_0_gt), TRUE)
})


test_that("tabu recovers the ground-truth discrete DBN structure", {
  set.seed(2)
  learned <- dbn.learn.structure(sampled_data, algorithm = "tabu",
                                       score = NULL)
  expect_s3_class(learned, "dbn")
  tn_l = bnlearn::model2network(modelstring(learned)$g_t)
  tn_gt = bnlearn::model2network(modelstring(dbn_truth)$g_t)
  expect_equal(all.equal(tn_l, tn_gt), TRUE)

  g_0_l = bnlearn::model2network(modelstring(learned)$g_0)
  g_0_gt = bnlearn::model2network(modelstring(dbn_truth)$g_0)
  expect_equal(all.equal(g_0_l, g_0_gt), TRUE)
})
