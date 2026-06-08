library(testthat)

# ----- fixture ----------------------------------------------------------------

build_conv_dbn <- function() {
  dbn <- empty.dbn(dynamic_nodes = c("A", "B"), markov_order = 1)
  dbn <- add.arc.dbn(dbn, from = c("A", "t_0"), to = c("B", "t_0"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t"),   to = c("B", "t"))
  dbn <- add.arc.dbn(dbn, from = c("A", "t-1"), to = c("A", "t"))
  dbn
}

build_conv_dbn_fit <- function() {
  dbn  <- build_conv_dbn()
  A_lv <- c("yes", "no")
  B_lv <- c("high", "low")
  dims_At <- list(A_t = A_lv); dims_At[["A_t-1"]] <- A_lv
  dbn.fit(DBN = dbn, CPTs = list(
    A_0 = array(c(0.3, 0.7), dim = 2, dimnames = list(A_0 = A_lv)),
    B_0 = array(c(0.6, 0.4, 0.2, 0.8), dim = c(2, 2),
                dimnames = list(B_0 = B_lv, A_0 = A_lv)),
    A_t = array(c(0.9, 0.1, 0.4, 0.6), dim = c(2, 2), dimnames = dims_At),
    B_t = array(c(0.7, 0.3, 0.1, 0.9), dim = c(2, 2),
                dimnames = list(B_t = B_lv, A_t = A_lv))
  ))
}

# ----- get.g0.net -------------------------------------------------------------

test_that("get.g0.net(DBN) returns a bn with correct nodes and arc", {
  g0 <- get.g0.net(build_conv_dbn())
  expect_equal(class(g0), "bn")
  expect_setequal(bnlearn::node.ordering(g0), c("A_0", "B_0"))
  expect_true(any(g0$arcs[, "from"] == "A_0" & g0$arcs[, "to"] == "B_0"))
})

test_that("get.g0.net(dbn.fit) returns a bn.fit with correct nodes", {
  g0 <- get.g0.net(build_conv_dbn_fit())
  expect_equal(class(g0), "bn.fit")
  expect_setequal(names(g0), c("A_0", "B_0"))
})

# ----- get.transition.net -----------------------------------------------------

test_that("get.transition.net(DBN) returns a bn with correct nodes and arcs", {
  tn <- get.transition.net(build_conv_dbn())
  expect_equal(class(tn), "bn")
  nodes <- bnlearn::node.ordering(tn)
  expect_true(all(c("A_t", "B_t", "A_t-1", "B_t-1") %in% nodes))
  expect_true(any(tn$arcs[, "from"] == "A_t"   & tn$arcs[, "to"] == "B_t"))
  expect_true(any(tn$arcs[, "from"] == "A_t-1" & tn$arcs[, "to"] == "A_t"))
})

test_that("get.transition.net(dbn.fit) returns a bn.fit with correct nodes", {
  tn <- get.transition.net(build_conv_dbn_fit())
  expect_equal(class(tn), "bn.fit")
  expect_setequal(names(tn), c("A_t", "B_t"))
})
