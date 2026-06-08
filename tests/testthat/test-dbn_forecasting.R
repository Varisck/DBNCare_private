
# DBN definition (mirrors test-dbn_sampling.R setup)
DBN_example <- empty.dbn(dynamic_nodes = c("A"), markov_order = 1)

# Node Addition
DBN_example <-
  add.node.dbn(DBN = DBN_example, node = "B", type = 'Dynamic')

# Arcs Addition
DBN_example <-
  add.arc.dbn(DBN = DBN_example,
              from = c('A', 't'),
              to = c('B', 't'))
DBN_example <-
  add.arc.dbn(DBN = DBN_example,
              from = c('A', 't_0'),
              to = c('B', 't_0'))
DBN_example <-
  add.arc.dbn(DBN = DBN_example,
              from = c('A', 't-1'),
              to = c('A', 't'))


# CPTs external definition
A_lv <- c('yes', 'no')
B_lv <- c('high', 'low')
A_0.prob <- array(c(0.2, 0.8),
                  dim = length(A_lv),
                  dimnames = list(A_0 = A_lv))
B_0.prob <- array(
  c(0.25, 0.75, 0.6, 0.4),
  dim = c(length(B_lv), length(A_lv)),
  dimnames = list(B_0 = B_lv, A_0 = A_lv)
)

dims <- list(A_t = A_lv)
dims[['A_t-1']] <- A_lv
A_t.prob = array(c(0.8, 0.2, 0.05, 0.95),
                 dim = c(length(A_lv), length(A_lv)),
                 dimnames = dims)
B_t.prob <- array(
  c(0.2, 0.8, 0.6, 0.4),
  dim = c(length(B_lv), length(A_lv)),
  dimnames = list(B_t = B_lv, A_t = A_lv)
)
CPTs_toy <- list(
  A_0 = A_0.prob,
  B_0 = B_0.prob,
  A_t = A_t.prob,
  B_t = B_t.prob
)
fitted_DBN <- dbn.fit(DBN = DBN_example, CPTs = CPTs_toy)


# testing dbn.forecasting input validation
test_that("dbn.forecasting raises error in case of wrong inputs", {
  observations <- list(Time = c(0),
                       Sample_id = c("sample1"),
                       A = c("yes"),
                       B = c("high"))

  # timepoints must be numeric, not a string
  expect_error(dbn.forecasting(fitted_DBN, observations, "3"))

  # timepoints must be greater than 0
  expect_error(dbn.forecasting(fitted_DBN, observations, 0))
  expect_error(dbn.forecasting(fitted_DBN, observations, -1))

  # the dbn argument must be of class dbn.fit
  expect_error(dbn.forecasting(c("not", "a", "dbn"), observations, 3))
  expect_error(dbn.forecasting(DBN_example, observations, 3))
})

test_that("dbn.forecasting returns a list with the expected structure", {
  observations <- list(Time = c(0),
                       Sample_id = c("sample1"),
                       A = c("yes"),
                       B = c("high"))
  timepoints <- 3

  result <- dbn.forecasting(fitted_DBN, observations, timepoints)

  # the result should be a list containing entries for every variable
  expect_true(is.list(result))
  expect_true(all(c("Time", "Sample_id", "A", "B") %in% names(result)))

  # the forecasting horizon should extend the observations by `timepoints`
  expected_length <- length(observations[["Time"]]) + timepoints
  expect_equal(length(result[["A"]]), expected_length)
  expect_equal(length(result[["B"]]), expected_length)
  expect_equal(length(result[["Time"]]), expected_length)
})

test_that("dbn.forecasting preserves the original observations", {
  observations <- list(Time = c(0),
                       Sample_id = c("sample1"),
                       A = c("yes"),
                       B = c("high"))

  result <- dbn.forecasting(fitted_DBN, observations, 2)

  # the first values must match the provided observations
  expect_equal(result[["A"]][1], "yes")
  expect_equal(result[["B"]][1], "high")
  expect_equal(result[["Time"]][1], 0)
})

test_that("dbn.forecasting produces values from the variable's domain", {
  observations <- list(Time = c(0),
                       Sample_id = c("sample1"),
                       A = c("yes"),
                       B = c("high"))

  result <- dbn.forecasting(fitted_DBN, observations, 5)

  # all forecasted values must come from the declared levels
  expect_true(all(result[["A"]] %in% A_lv))
  expect_true(all(result[["B"]] %in% B_lv))
})

test_that(
  "dbn.forecasting forecasts converge to the transition CPT in distribution",
  {
    # When A_t-1 = 'yes', P(A_t = 'yes') = 0.8 according to A_t.prob.
    # Repeating one-step forecasts from the same observation should yield
    # an empirical frequency of 'yes' close to 0.8.
    set.seed(42)
    n_runs <- 1e+3
    tolerance <- 0.05

    observations <- list(Time = c(0),
                         Sample_id = c("sample1"),
                         A = c("yes"),
                         B = c("high"))

    forecasted_A <- character(n_runs)
    for (i in seq_len(n_runs)) {
      out <- dbn.forecasting(fitted_DBN, observations, 1)
      forecasted_A[i] <- out[["A"]][2]
    }

    empirical_yes <- mean(forecasted_A == "yes")
    expect_true(abs(empirical_yes - 0.8) <= tolerance)
  }
)

test_that(
  "dbn.forecasting works on a multi-node DBN and returns full trajectories",
  {
    # Build a 3-node DBN equivalent to the second sampling test
    my_dbn <-
      empty.dbn(dynamic_nodes = c("A", "B", "C"),
                markov_order = 1)

    # G_0
    my_dbn <- add.arc.dbn(my_dbn, from = c('A', 't_0'), to = c('B', 't_0'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('A', 't_0'), to = c('C', 't_0'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('B', 't_0'), to = c('C', 't_0'))

    # G_transition
    my_dbn <- add.arc.dbn(my_dbn, from = c('A', 't'), to = c('B', 't'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('A', 't'), to = c('C', 't'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('B', 't'), to = c('C', 't'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('A', 't-1'), to = c('A', 't'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('B', 't-1'), to = c('B', 't'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('C', 't-1'), to = c('C', 't'))

    # CPTs G_0
    A_lv <- c('yes', 'no')
    B_lv <- c('high', 'low')
    C_lv <- c('medium', 'high')
    A_0.prob <- array(c(0.2, 0.8),
                      dim = length(A_lv),
                      dimnames = list(A_0 = A_lv))
    B_0.prob <- array(
      c(0.25, 0.75, 0.6, 0.4),
      dim = c(length(B_lv), length(A_lv)),
      dimnames = list(B_0 = B_lv, A_0 = A_lv)
    )
    C_0.prob <- array(
      c(0.1, 0.9, 0.5, 0.5, 0.3, 0.7, 0.6, 0.4),
      dim = c(length(C_lv), length(B_lv), length(A_lv)),
      dimnames = list(C_0 = C_lv, B_0 = B_lv, A_0 = A_lv)
    )

    # CPTs G_transition
    dims_A_t <- list(A_t = A_lv)
    dims_A_t[['A_t-1']] = A_lv
    dims_B_t <- list(B_t = B_lv, A_t = A_lv)
    dims_B_t[['B_t-1']] = B_lv
    dims_C_t <- list(C_t = C_lv, A_t = A_lv, B_t = B_lv)
    dims_C_t[['C_t-1']] = C_lv

    A_t.prob <- array(c(0.8, 0.2, 0.05, 0.95),
                      dim = c(length(A_lv), length(A_lv)),
                      dimnames = dims_A_t)
    B_t.prob <- array(
      c(0.2, 0.8, 0.6, 0.4, 0.3, 0.7, 0.1, 0.9),
      dim = c(length(B_lv), length(A_lv), length(B_lv)),
      dimnames = dims_B_t
    )
    C_t.prob <- array(
      c(0.1, 0.9, 0.2, 0.8, 0.12, 0.88, 0.9, 0.1,
        0.3, 0.7, 0.4, 0.6, 0.25, 0.75, 0.45, 0.55),
      dim = c(length(C_lv), length(A_lv), length(B_lv), length(C_lv)),
      dimnames = dims_C_t
    )

    CPTs_mydbn <- list(
      A_0 = A_0.prob, B_0 = B_0.prob, C_0 = C_0.prob,
      A_t = A_t.prob, B_t = B_t.prob, C_t = C_t.prob
    )
    fitted_my_dbn <- dbn.fit(my_dbn, CPTs_mydbn)

    observations <- list(Time = c(0),
                         Sample_id = c("sample1"),
                         A = c("no"),
                         B = c("low"),
                         C = c("medium"))
    timepoints <- 3

    result <- dbn.forecasting(fitted_my_dbn, observations, timepoints)

    expected_length <- length(observations[["Time"]]) + timepoints
    expect_equal(length(result[["A"]]), expected_length)
    expect_equal(length(result[["B"]]), expected_length)
    expect_equal(length(result[["C"]]), expected_length)

    # initial observation preserved
    expect_equal(result[["A"]][1], "no")
    expect_equal(result[["B"]][1], "low")
    expect_equal(result[["C"]][1], "medium")

    # forecasted values respect the variables' domains
    expect_true(all(result[["A"]] %in% A_lv))
    expect_true(all(result[["B"]] %in% B_lv))
    expect_true(all(result[["C"]] %in% C_lv))
  }
)

test_that(
  "dbn.forecasting on a long horizon reproduces the transition CPTs",
  {
    # Forecast a single long trajectory and check that the empirical
    # conditional frequencies along it match the true transition CPTs
    # used to fit the DBN.
    set.seed(42)
    horizon <- 1e+4
    tolerance <- 0.02

    observations <- list(Time = c(0),
                         Sample_id = c("sample1"),
                         A = c("yes"),
                         B = c("high"))

    res <- dbn.forecasting(fitted_DBN, observations, horizon)

    compare_with_tolerance <-
      function(value1, value2, tolerance)
        abs(value1 - value2) <= tolerance

    A <- res[["A"]]
    B <- res[["B"]]
    A_prev <- A[1:(length(A) - 1)]
    A_curr <- A[2:length(A)]
    B_curr <- B[2:length(B)]

    # P(A_t | A_t-1) — A_t.prob
    expect_true(compare_with_tolerance(
      mean(A_curr[A_prev == "yes"] == "yes"),
      fitted_DBN$A_t$prob["yes", "yes"], tolerance))
    expect_true(compare_with_tolerance(
      mean(A_curr[A_prev == "yes"] == "no"),
      fitted_DBN$A_t$prob["no", "yes"], tolerance))
    expect_true(compare_with_tolerance(
      mean(A_curr[A_prev == "no"] == "yes"),
      fitted_DBN$A_t$prob["yes", "no"], tolerance))
    expect_true(compare_with_tolerance(
      mean(A_curr[A_prev == "no"] == "no"),
      fitted_DBN$A_t$prob["no", "no"], tolerance))

    # P(B_t | A_t) — B_t.prob
    expect_true(compare_with_tolerance(
      mean(B_curr[A_curr == "yes"] == "high"),
      fitted_DBN$B_t$prob["high", "yes"], tolerance))
    expect_true(compare_with_tolerance(
      mean(B_curr[A_curr == "yes"] == "low"),
      fitted_DBN$B_t$prob["low", "yes"], tolerance))
    expect_true(compare_with_tolerance(
      mean(B_curr[A_curr == "no"] == "high"),
      fitted_DBN$B_t$prob["high", "no"], tolerance))
    expect_true(compare_with_tolerance(
      mean(B_curr[A_curr == "no"] == "low"),
      fitted_DBN$B_t$prob["low", "no"], tolerance))
  }
)

test_that(
  "dbn.forecasting on a long horizon reproduces CPTs on a multi-node DBN",
  {
    # Build a 3-node DBN and check that a long forecast recovers the
    # transition CPTs of A_t (parent: A_t-1) and B_t (parents: A_t, B_t-1)
    # along the trajectory.
    my_dbn <-
      empty.dbn(dynamic_nodes = c("A", "B", "C"),
                markov_order = 1)

    my_dbn <- add.arc.dbn(my_dbn, from = c('A', 't_0'), to = c('B', 't_0'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('A', 't_0'), to = c('C', 't_0'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('B', 't_0'), to = c('C', 't_0'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('A', 't'), to = c('B', 't'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('A', 't'), to = c('C', 't'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('B', 't'), to = c('C', 't'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('A', 't-1'), to = c('A', 't'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('B', 't-1'), to = c('B', 't'))
    my_dbn <- add.arc.dbn(my_dbn, from = c('C', 't-1'), to = c('C', 't'))

    A_lv <- c('yes', 'no')
    B_lv <- c('high', 'low')
    C_lv <- c('medium', 'high')
    A_0.prob <- array(c(0.2, 0.8),
                      dim = length(A_lv),
                      dimnames = list(A_0 = A_lv))
    B_0.prob <- array(
      c(0.25, 0.75, 0.6, 0.4),
      dim = c(length(B_lv), length(A_lv)),
      dimnames = list(B_0 = B_lv, A_0 = A_lv)
    )
    C_0.prob <- array(
      c(0.1, 0.9, 0.5, 0.5, 0.3, 0.7, 0.6, 0.4),
      dim = c(length(C_lv), length(B_lv), length(A_lv)),
      dimnames = list(C_0 = C_lv, B_0 = B_lv, A_0 = A_lv)
    )

    dims_A_t <- list(A_t = A_lv); dims_A_t[['A_t-1']] = A_lv
    dims_B_t <- list(B_t = B_lv, A_t = A_lv); dims_B_t[['B_t-1']] = B_lv
    dims_C_t <- list(C_t = C_lv, A_t = A_lv, B_t = B_lv); dims_C_t[['C_t-1']] = C_lv

    A_t.prob <- array(c(0.8, 0.2, 0.05, 0.95),
                      dim = c(length(A_lv), length(A_lv)),
                      dimnames = dims_A_t)
    B_t.prob <- array(
      c(0.2, 0.8, 0.6, 0.4, 0.3, 0.7, 0.1, 0.9),
      dim = c(length(B_lv), length(A_lv), length(B_lv)),
      dimnames = dims_B_t
    )
    C_t.prob <- array(
      c(0.1, 0.9, 0.2, 0.8, 0.12, 0.88, 0.9, 0.1,
        0.3, 0.7, 0.4, 0.6, 0.25, 0.75, 0.45, 0.55),
      dim = c(length(C_lv), length(A_lv), length(B_lv), length(C_lv)),
      dimnames = dims_C_t
    )

    fitted_my_dbn <- dbn.fit(my_dbn, list(
      A_0 = A_0.prob, B_0 = B_0.prob, C_0 = C_0.prob,
      A_t = A_t.prob, B_t = B_t.prob, C_t = C_t.prob
    ))

    set.seed(7)
    horizon <- 2e+4
    tolerance <- 0.03

    observations <- list(Time = c(0),
                         Sample_id = c("sample1"),
                         A = c("yes"), B = c("high"), C = c("medium"))

    res <- dbn.forecasting(fitted_my_dbn, observations, horizon)

    compare_with_tolerance <-
      function(value1, value2, tolerance)
        abs(value1 - value2) <= tolerance

    A <- res[["A"]]; B <- res[["B"]]
    A_prev <- A[1:(length(A) - 1)]
    A_curr <- A[2:length(A)]
    B_prev <- B[1:(length(B) - 1)]
    B_curr <- B[2:length(B)]

    # P(A_t | A_t-1)
    expect_true(compare_with_tolerance(
      mean(A_curr[A_prev == "yes"] == "yes"),
      fitted_my_dbn$A_t$prob["yes", "yes"], tolerance))
    expect_true(compare_with_tolerance(
      mean(A_curr[A_prev == "no"] == "yes"),
      fitted_my_dbn$A_t$prob["yes", "no"], tolerance))

    # P(B_t | A_t, B_t-1) — a few representative cells
    mask <- A_curr == "yes" & B_prev == "high"
    expect_true(compare_with_tolerance(
      mean(B_curr[mask] == "high"),
      fitted_my_dbn$B_t$prob["high", "yes", "high"], tolerance))

    mask <- A_curr == "no" & B_prev == "high"
    expect_true(compare_with_tolerance(
      mean(B_curr[mask] == "high"),
      fitted_my_dbn$B_t$prob["high", "no", "high"], tolerance))

    mask <- A_curr == "yes" & B_prev == "low"
    expect_true(compare_with_tolerance(
      mean(B_curr[mask] == "high"),
      fitted_my_dbn$B_t$prob["high", "yes", "low"], tolerance))

    mask <- A_curr == "no" & B_prev == "low"
    expect_true(compare_with_tolerance(
      mean(B_curr[mask] == "high"),
      fitted_my_dbn$B_t$prob["high", "no", "low"], tolerance))
  }
)
