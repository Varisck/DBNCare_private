test_that("split_variable_name returns the correct output", {
  expect_equal(split_variable_name("A_t-1"), list("name" = "A", "time" = "t-1"))
  expect_equal(split_variable_name("B_0"), list("name" = "B", "time" = "t_0"))
  expect_equal(split_variable_name("C_t"), list("name" = "C", "time" = "t"))
})
test_that("split_variable_name raises error when wrong input is given", {
  expect_error(split_variable_name(1))
})
test_that("generate_dbn_nodes_names returns the correct output", {
  expect_equal(generate_dbn_nodes_names(c("A", "B", "C")),
               list(
                 g_0 = c("A_0", "B_0", "C_0"),
                 g_t_1 = c("A_t-1", "B_t-1", "C_t-1"),
                 g_t = c("A_t", "B_t", "C_t")
               ))
})
test_that("generate_dbn_nodes_names raises error when wrong input is given",
          {
            #input not a vector
            expect_error(generate_dbn_nodes_names("A"))
            expect_error(generate_dbn_nodes_names(c(1, 2)))
          })
test_that("random.structure.dbn returns the correct output", {
  check_correspondence <- function(dbn_edges) {
    # Extract pairs ending with _0 and _t
    pairs_0 <-
      unique(dbn_edges[grepl("_0$", dbn_edges[, 1]) &
                         grepl("_0$", dbn_edges[, 2]), , drop = FALSE])
    pairs_t <-
      unique(dbn_edges[grepl("_t$", dbn_edges[, 1]) &
                         grepl("_t$", dbn_edges[, 2]), , drop = FALSE])

    # Function to create a set of pairs without suffix
    remove_suffix <- function(pairs) {
      unique(paste0(
        sub("_0$|_t$", "", pairs[, 1]),
        "_",
        sub("_0$|_t$", "", pairs[, 2])
      ))
    }

    # Create sets of pairs without suffix
    set_0 <- remove_suffix(pairs_0)
    set_t <- remove_suffix(pairs_t)

    # Check if all elements in set_0 are in set_t
    all(set_0 %in% set_t)
  }
  #g_0 does contain the same edges of g_transition when g_0_arcs = TRUE
  dbn_random_generated <-
    random.structure.dbn(c("A", "B", "C", "D"), 0.5, 0.5, markov_order = 1, g_0_arcs = TRUE)
  dbn_edges <- dbn_random_generated$arcs
  check_result <- check_correspondence(dbn_edges)
  expect_equal(check_result, TRUE)
  set.seed(2131)
  #g_0 does not contain the same edges of g_transition when g_0_arcs = FALSE
  dbn_random_generated <-
    random.structure.dbn(c("A", "B", "C", "D"), 0.2, 0.2, markov_order = 1,
                         g_0_arcs = FALSE)
  dbn_edges <- dbn_random_generated$arcs
  check_result <- check_correspondence(dbn_edges)
  expect_equal(check_result, FALSE)

  #the function returns a DBN object
  expect_equal(class(random.structure.dbn(c("A", "B"), 0.5, 0.5, markov_order = 1)), "dbn")

  #no error when the input parameters are correct
  expect_no_error(random.structure.dbn(c("A", "B", "C", "D"), 0.6, 0.6, markov_order = 1))
  expect_no_error(random.structure.dbn(c("A", "B", "C", "D"), 0.6, 0.6, markov_order = 1,
                                       g_0_arcs = FALSE, g_0_prob = 0.0))

  # markov_order > 1: per-lag interslice probability vector
  expect_no_error(random.structure.dbn(c("A", "B", "C"), 0.3, c(0.5, 0.2),
                                       markov_order = 2))
})
test_that("random.structure.dbn raises error when wrong input is given",
          {
            #input is not a char vector
            expect_error(random.structure.dbn(c(1, 2), 0.5, 0.5, markov_order = 1))
            #g_0_arcs is not a boolean
            expect_error(random.structure.dbn(c("A", "B"), 0.5, 0.5, markov_order = 1,
                                              g_0_arcs = "Not a boolean"))
            #intra-slice prob outside [0, 1]
            expect_error(random.structure.dbn(c("A", "B"), 10, 0.5, markov_order = 1))
            expect_error(random.structure.dbn(c("A", "B"), -0.1, 0.5, markov_order = 1))
            #inter-slice prob outside [0, 1]
            expect_error(random.structure.dbn(c("A", "B"), 0.5, 10, markov_order = 1))
            expect_error(random.structure.dbn(c("A", "B"), 0.5, -0.1, markov_order = 1))
            #markov_order must be a positive integer
            expect_error(random.structure.dbn(c("A", "B"), 0.5, 0.5, markov_order = 0))
            expect_error(random.structure.dbn(c("A", "B"), 0.5, 0.5, markov_order = 1.5))
            #length(prob_edge_interslice) must be 1 or markov_order
            expect_error(random.structure.dbn(c("A", "B"), 0.5, c(0.5, 0.3),
                                              markov_order = 1))
          })

test_that("generate_dbn_nodes_distributions raises error when wrong input is given",
          {
            #generated_dbn is not a DBN object
            expect_error(generate_dbn_nodes_distributions("A", TRUE, 2))
            generated_random_dbn <-
              random.structure.dbn(c("A", "B"), 0.5, 0.5, markov_order = 1)
            #fixed_cardinality not a boolean
            expect_error(generate_dbn_nodes_distributions(generated_random_dbn, "Not boolean", 2))
            #max_variables_cardinality <2
            expect_error(generate_dbn_nodes_distributions(generated_random_dbn, TRUE, 1))

          })

test_that("generate_dbn_nodes_distributions correctly istantiates nodes' distribution",
          {
            #no errors is produced when the input is correct
            generated_random_dbn <-
              random.structure.dbn(c("A", "B", "C"), 0.6, 0.6, markov_order = 1)
            expect_no_error(generate_dbn_nodes_distributions(generated_random_dbn, TRUE, 2))
            #variables levels are aligned with max_cardinality when fixed_cardinality = TRUE
            generated_random_dbn <-
              random.structure.dbn(c("A", "B", "C"), 0.6, 0.6, markov_order = 1)
            fitted_random_dbn <-
              generate_dbn_nodes_distributions(generated_random_dbn, TRUE, 2)
            all_variables <- names(fitted_random_dbn)
            correct_instantiation <- TRUE
            for (var in all_variables) {
              if (length(dimnames(fitted_random_dbn$B_t$prob)[[1]]) != 2) {
                correct_instantiation <- FALSE
              }
            }
            expect_equal(correct_instantiation, TRUE)
            #variables levels are aligned with max_cardinality when fixed_cardinality = FALSE
            fitted_random_dbn <-
              generate_dbn_nodes_distributions(generated_random_dbn, FALSE, 4)
            all_variables <- names(fitted_random_dbn)
            correct_instantiation <- TRUE
            for (var in all_variables) {
              if (!(length(dimnames(fitted_random_dbn[[var]][["prob"]])[[1]]) >= 2 &
                    length(dimnames(fitted_random_dbn[[var]][["prob"]])[[1]]) <= 4)) {
                correct_instantiation <- FALSE
              }
            }
            expect_equal(correct_instantiation, TRUE)
          })
