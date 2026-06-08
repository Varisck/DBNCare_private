DBN <- DynamicBayesianNetwork::random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
DBN_fitted <- DynamicBayesianNetwork::generate_dbn_nodes_distributions(DBN, TRUE, 2)
set.seed(123)
data <- DynamicBayesianNetwork::dbn.sampling(DBN_fitted, 2000, 5)

test_that("dbn.prep.data runs without error", {
  expect_no_error(DynamicBayesianNetwork::dbn.prep.data(data))
})

test_that("dbn.learn.structure runs with default settings", {
  expect_no_error(DynamicBayesianNetwork::dbn.learn.structure(data))
})

test_that("dbn.learn.structure accepts all score-based scores", {
  for (i in c("loglik", "aic", "bic", "ebic", "bde", "bds", "mbde", "k2", "fnml", "qnml", "nal", "pnal")) {
    expect_no_error(
      DynamicBayesianNetwork::dbn.learn.structure(data, score = i)
    )
  }
  expect_warning(
    DynamicBayesianNetwork::dbn.learn.structure(data, score = "bdla"),
    regexp = "bdla.*deprecated"
  )
})

test_that("dbn.learn.structure accepts all constraint-based tests via pc.stable", {
  for (i in c("mi", "mi-adf", "mc-mi", "smc-mi", "sp-mi", "mi-sh", "x2-adf", "mc-x2", "smc-x2", "sp-x2")) {
    expect_no_error(
      DynamicBayesianNetwork::dbn.learn.structure(data, algorithm = "pc.stable", test = i)
    )
  }
})

test_that("dbn.learn.structure accepts score-based algorithms", {
  for (i in c("hc", "tabu")) {
    expect_no_error(
      DynamicBayesianNetwork::dbn.learn.structure(data, algorithm = i)
    )
  }
})

test_that("dbn.learn.structure accepts constraint-based algorithms", {
  for (i in c("pc.stable", "gs", "iamb", "fast.iamb", "inter.iamb", "iamb.fdr")) {
    expect_no_error(
      DynamicBayesianNetwork::dbn.learn.structure(data, algorithm = i)
    )
  }
})

test_that("dbn.learn.structure accepts hybrid algorithms", {
  for (i in c("h2pc", "mmhc", "rsmax2")) {
    expect_no_error(
      DynamicBayesianNetwork::dbn.learn.structure(data, algorithm = i)
    )
  }
})
