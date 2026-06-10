#' generate_dbn_nodes_names
#'
#' @description
#' the function returns the nodes names for g_0 and g_transition.
#'
#' @param nodes_names a vector containing the names of nodes.
#'
#' @return a list containing the nodes names for each time slice.
#' @export
#'
#' @examples generate_dbn_nodes_names(c("A", "B", "C"))
generate_dbn_nodes_names <- function(nodes_names) {
  if (!is.vector(nodes_names) ||
      !is.character(nodes_names) || length(nodes_names) <= 1) {
    stop("nodes_names must be a vector containing names of nodes!")
  }
  g_0_nodes_names <- c()
  g_t_1_nodes_names <- c()
  g_t_nodes_names <- c()
  
  for (n in nodes_names) {
    g_0_nodes_names <- c(g_0_nodes_names, paste(n, "_0", sep = ""))
    g_t_1_nodes_names <-
      c(g_t_1_nodes_names, paste(n, "_t-1", sep = ""))
    g_t_nodes_names <- c(g_t_nodes_names, paste(n, "_t", sep = ""))
    
  }
  list_of_nodes_names <-
    list(g_0 = g_0_nodes_names, g_t_1 = g_t_1_nodes_names, g_t = g_t_nodes_names)
  return(list_of_nodes_names)
}

#' generate_dbn_nodes_distributions
#'
#' @description
#' This function generates a set of multinomial distributions for each node in a network,
#' using the Dirichlet distribution to ensure uniform sampling of probability values.
#' The generated distributions can either have a fixed cardinality for all nodes or a
#' variable cardinality for each node, depending on the specified parameters.
#'
#'
#' @param generated_dbn an object of class DBN
#' @param fixed_cardinality A logical parameter (`TRUE` or `FALSE`). If `TRUE`,
#' all nodes will have the same cardinality as specified by
#' `max_variables_cardinality`. If `FALSE`, the cardinality for each node is randomly
#' determined, with a minimum of 2 and a maximum defined by `max_variables_cardinality`.
#' @param max_variables_cardinality An integer specifying the maximum cardinality to be used for the nodes.
#' If `fixed_cardinality` is `TRUE`, this parameter sets the exact cardinality for all nodes. If `fixed_cardinality`
#' is `FALSE`, this parameter defines the upper limit of the range from which the cardinality is randomly selected.
#'
#' @return a DBN.fit object
#' @export
#'
#' @examples generate_dbn_nodes_distributions(generated_dbn,c("A","B","C"), TRUE, 2)
generate_dbn_nodes_distributions <-
  function(generated_dbn,
           fixed_cardinality,
           max_variables_cardinality) {
    if (class(generated_dbn) != "dbn") {
      stop("generated_dbn must be a DBN object!")
    }
    if (!is.logical(fixed_cardinality)) {
      stop("fixed_cardinality must be TRUE or FALSE.")
    }
    if (max_variables_cardinality < 2) {
      stop("max_variables_cardinality must be greater than or equal to 2.")
    }
    nodes_names <- names(generated_dbn$nodes)
    list_of_nodes_names <- generate_dbn_nodes_names(nodes_names)
    g_0_nodes_names <- list_of_nodes_names[['g_0']]
    g_t_1_nodes_names <- list_of_nodes_names[['g_t_1']]
    g_t_nodes_names <- list_of_nodes_names[['g_t']]
    
    #generating levels for each variable
    list_of_levels <- list()
    #fixed_cardinality == TRUE all the nodes must have a fixed cardinality
    #the cardinality is equal to cardinality parameter
    if (fixed_cardinality == TRUE) {
      for (n in nodes_names) {
        list_of_levels[[n]] <- 0:(max_variables_cardinality - 1)
      }
    } else{
      #fixed_cardinality == FALSE for each node we will have a random cardinality
      # 2<= cardinality <=max_variables_cardinality
      for (n in nodes_names) {
        if (max_variables_cardinality == 2) {
          cardinality_n <- 1
        }
        else{
          cardinality_n <- sample(2:max_variables_cardinality, 1) - 1
        }
        list_of_levels[[n]] <- 0:cardinality_n
      }
    }
    #generate distribution for t_0
    list_cpt <- list()
    for (n_0 in g_0_nodes_names) {
      root_node_name <- split_variable_name(n_0)$name
      dim_names_n_0 <- list()
      cardinality_cpt_n_0 <- c()
      dim_names_n_0[[n_0]] <- list_of_levels[[root_node_name]]
      cardinality_cpt_n_0 <-
        c(cardinality_cpt_n_0, length(list_of_levels[[root_node_name]]))
      parents_n_0 <-
        generated_dbn[["nodes"]][[root_node_name]][["t_0"]][["parents"]]
      
      #adding each node's and parent's levels and creating the list of cardinality
      n_distributions <- 1
      
      for (p_n_0 in parents_n_0) {
        root_node_name_p <- split_variable_name(p_n_0)$name
        dim_names_n_0[[p_n_0]] <- list_of_levels[[root_node_name_p]]
        cardinality_cpt_n_0 <-
          c(cardinality_cpt_n_0, length(list_of_levels[[root_node_name_p]]))
        n_distributions <-
          n_distributions * length(list_of_levels[[root_node_name_p]])
      }
      k <- length(list_of_levels[[root_node_name]])
      alpha <-
        rep(1, k)  # Parameters for Dirichlet distribution, can be adjusted
      vector_random_distributions <- c()
      for (i in seq(1, n_distributions)) {
        random_probabilities <- MCMCpack::rdirichlet(1, alpha)
        distribution_n_0 <- as.vector(random_probabilities)
        vector_random_distributions <-
          c(vector_random_distributions, distribution_n_0)
      }
      list_cpt[[n_0]] <-
        array(vector_random_distributions,
              cardinality_cpt_n_0,
              dim_names_n_0)
    }
    #Generate distribution for t
    for (n_t in g_t_nodes_names) {
      root_node_name <- split_variable_name(n_t)$name
      dim_names_n_t <- list()
      cardinality_cpt_n_t <- c()
      dim_names_n_t[[n_t]] <- list_of_levels[[root_node_name]]
      cardinality_cpt_n_t <-
        c(cardinality_cpt_n_t, length(list_of_levels[[root_node_name]]))
      parents_n_t <-
        generated_dbn[["nodes"]][[root_node_name]][["t"]][["parents"]]
      
      #adding each node's and parent's levels and creating the list of cardinality
      n_distributions <- 1
      
      for (p_n_t in parents_n_t) {
        root_node_name_p <- split_variable_name(p_n_t)$name
        dim_names_n_t[[p_n_t]] <- list_of_levels[[root_node_name_p]]
        cardinality_cpt_n_t <-
          c(cardinality_cpt_n_t, length(list_of_levels[[root_node_name_p]]))
        n_distributions <-
          n_distributions * length(list_of_levels[[root_node_name_p]])
      }
      k <- length(list_of_levels[[root_node_name]])
      alpha <-
        rep(1, k)  # Parameters for Dirichlet distribution, can be adjusted
      vector_random_distributions <- c()
      for (i in seq(1, n_distributions)) {
        random_probabilities <- MCMCpack::rdirichlet(1, alpha)
        distribution_n_t <- as.vector(random_probabilities)
        vector_random_distributions <-
          c(vector_random_distributions, distribution_n_t)
      }
      list_cpt[[n_t]] <-
        array(vector_random_distributions,
              cardinality_cpt_n_t,
              dim_names_n_t)
    }
    CPTs <- list_cpt
    fitted_dbn <- dbn.fit(DBN = generated_dbn, CPTs = list_cpt)
    
    return(fitted_dbn)
  }

#' random.dbn
#'
#' @description
#' Thin convenience wrapper around [random.structure.dbn()]. Generates a random
#' DBN structure of arbitrary Markov order following the node ordering in
#' `nodes_names`.
#'
#' @param nodes_names a character vector of node names.
#' @param is_same TRUE if G_0 should mirror the intra-slice arcs of G_t (mapped to
#'   `g_0_arcs = TRUE` in [random.structure.dbn()]). FALSE samples G_0 independently
#'   with edge probability `prob_edges_g0`. Defaults to TRUE.
#' @param prob_edge_intraslice probability of an intra-slice edge (X_t -> Y_t).
#'   Defaults to `2 / (length(nodes_names) - 1)`, which targets ~2 expected
#'   intra-slice parents per node in a sparse random graph.
#' @param prob_edge_interslice probability of an inter-slice edge (X_t-k -> Y_t).
#'   Either a single number applied to every lag, or a length-`markov_order`
#'   vector of per-lag probabilities. Defaults to `2 / (length(nodes_names) - 1)`.
#' @param prob_edges_g0 probability of an arc in the prior network G_0; only used
#'   when `is_same = FALSE`. Passed through to [random.structure.dbn()] as
#'   `g_0_prob`. Defaults to `prob_edge_intraslice`.
#' @param markov_order positive integer Markov order of the generated DBN.
#'   Defaults to 1.
#'
#' @return a 'dbn' object.
#' @export
#'
#' @examples
#' random.dbn(c("A","B","C","D","E"))                    # uses all defaults
#' random.dbn(c("A","B","C"), TRUE, 0.6, 0.5)
#' random.dbn(c("A","B","C","D"), is_same = FALSE,
#'            prob_edge_intraslice = 0.4, prob_edge_interslice = 0.4,
#'            prob_edges_g0 = 0.2)
#' # markov_order > 1 with per-lag interslice probabilities
#' random.dbn(c("A","B","C"), prob_edge_intraslice = 0.3,
#'            prob_edge_interslice = c(0.5, 0.2), markov_order = 2)
random.dbn <- function(nodes_names,
                       is_same              = TRUE,
                       prob_edge_intraslice = 2 / (length(nodes_names) - 1),
                       prob_edge_interslice = 2 / (length(nodes_names) - 1),
                       prob_edges_g0        = prob_edge_intraslice,
                       markov_order         = 1) {
  random.structure.dbn(
    node_names           = nodes_names,
    prob_edge_intraslice = prob_edge_intraslice,
    prob_edge_interslice = prob_edge_interslice,
    markov_order         = markov_order,
    g_0_arcs             = is_same,
    g_0_prob             = prob_edges_g0
  )
}

#' fit_random_dbn
#'
#' @description
#' Fits random parameters to a DBN structure. Dispatches to the appropriate
#' subroutine based on `type`:
#' - `"discrete"`: random multinomial CPTs sampled from a Dirichlet prior
#'   (see [generate_dbn_nodes_distributions()]).
#' - `"continuous"`: random gaussian regression coefficients sampled uniformly
#'   on `[param.lower, param.upper]` with fixed residual sd.
#' - `"mixed"`: not yet implemented.
#'
#' @param generated_dbn a DBN class object
#' @param type one of `"discrete"`, `"continuous"`, or `"mixed"`.
#' @param fixed_cardinality (discrete only) A logical parameter (`TRUE` or `FALSE`). If `TRUE`,
#' all nodes will have the same cardinality as specified by
#' `max_variables_cardinality`. If `FALSE`, the cardinality for each node is randomly
#' determined, with a minimum of 2 and a maximum defined by `max_variables_cardinality`.
#' @param max_variables_cardinality (discrete only) An integer specifying the maximum cardinality to be used for the nodes.
#' If `fixed_cardinality` is `TRUE`, this parameter sets the exact cardinality for all nodes. If `fixed_cardinality`
#' is `FALSE`, this parameter defines the upper limit of the range from which the cardinality is randomly selected.
#' @param param.lower (continuous only) lower bound of the uniform from which regression coefficients are sampled.
#' @param param.upper (continuous only) upper bound of the uniform from which regression coefficients are sampled.
#' @param sd (continuous only) residual standard deviation assigned to every node.
#'
#' @return a dbn.fit object
#' @export
#'
#' @examples
#' \dontrun{
#' fit_random_dbn(dbn, type = "discrete", max_variables_cardinality = 2, fixed_cardinality = TRUE)
#' fit_random_dbn(dbn, type = "continuous", param.lower = .2, param.upper = .8, sd = 1)
#' }
fit_random_dbn <-
  function(generated_dbn,
           type = c("discrete", "continuous", "mixed"),
           max_variables_cardinality = 2,
           fixed_cardinality = FALSE,
           param.lower = .2,
           param.upper = .8,
           sd = 1) {
    type <- match.arg(type)

    if (type == "discrete") {
      return(generate_dbn_nodes_distributions(generated_dbn,
                                              fixed_cardinality,
                                              max_variables_cardinality))
    }
    if (type == "continuous") {
      return(fit_random_dbn_g(generated_dbn,
                              param.lower = param.lower,
                              param.upper = param.upper,
                              sd = sd))
    }
    if (type == "mixed") {
      stop("to be implemented")
    }
  }




fit_random_dbn_g = function(dbn, param.lower = .2, param.upper = .8,
                            sd = 1) {

  # fitted dbn
  node_info = list()
  variables = names(dbn$nodes)
  markov_order = dbn$markov_order

  # t_0 and t carry the full regression (intercept + parents).
  times = c("t_0", "t")

  for(variable in variables) {
    for(t in times) {
      parents = dbn$nodes[[variable]][[t]]$parents
      children = dbn$nodes[[variable]][[t]]$children
      p_len = length(parents)
      regs = runif(p_len + 1, param.lower, param.upper)
      names(regs) = c(intercept_name, parents)
      variable_name = gsub(" ", "", paste(variable, ifelse(t == "t_0", "_0", "_t")))

      node_info[[variable_name]] = list(
        node = variable_name,
        parents = parents,
        children = children,
        regs = regs,
        std = sd
      )
    }
  }

  # markov_order > 1 needs marginal priors for the initial slices
  # variable_1, ..., variable_(markov_order - 1) (no parents, just intercept).
  # Matches the shape produced by learn_param_g_data / get_var_priors.
  if (markov_order > 1) {
    priors = list()
    for (mo in seq_len(markov_order - 1)) {
      for (variable in variables) {
        coeff = runif(1, param.lower, param.upper)
        names(coeff) = intercept_name
        priors[[paste0(variable, "_", mo)]] = list(
          coeff = coeff,
          std = sd
        )
      }
    }
    node_info[["priors"]] = priors
  }

  class(node_info) = "dbn.fit"
  node_info
}

#' random.structure.dbn
#'
#' @description
#' Generates a random transition network (\eqn{G_t}) for a Dynamic Bayesian
#' Network of arbitrary Markov order. Two families of arcs are sampled
#' independently:
#' \itemize{
#'   \item \strong{Intra-slice arcs} \eqn{X_t \to Y_t}: arcs within the current
#'     time slice \eqn{t}. Each candidate arc is included with probability
#'     \code{prob_edge_intraslice}. Intra-slice arcs exist \emph{only} at slice
#'     \eqn{t}: past slices (\eqn{t-k}) carry no intra-slice arcs because they
#'     act as observed/frozen inputs to the transition model.
#'   \item \strong{Inter-slice (temporal) arcs} \eqn{X_{t-k} \to Y_t}: arcs from
#'     a past slice \eqn{t-k} (\eqn{1 \le k \le} \code{markov_order}) into the
#'     current slice. Each candidate arc is included with the probability
#'     associated to its lag (see \code{prob_edge_interslice}).
#' }
#'
#' @details
#' The generated graph is always a valid DAG. Inter-slice arcs always point
#' forward in time and no arc ever enters a past-slice node, so past-slice nodes
#' can never belong to a cycle; the only place a cycle could appear is among the
#' intra-slice arcs at \eqn{t}, which are sampled following the node ordering
#' implied by \code{node_names} (via \code{bnlearn::random.graph}), which
#' guarantees acyclicity. Self temporal arcs (\eqn{X_{t-k} \to X_t}) are allowed
#' and are typically desirable.
#'
#' The prior network (\eqn{G_0}) is populated according to \code{g_0_arcs}: when
#' \code{TRUE} (default) it mirrors the slice-\eqn{t} intra-slice structure
#' (every arc \eqn{X_t \to Y_t} also appears as \eqn{X_0 \to Y_0}); when
#' \code{FALSE} an independent random DAG is sampled over the \eqn{G_0} nodes
#' with edge probability \code{g_0_prob}. Both routes yield an acyclic \eqn{G_0},
#' as they reuse / follow the \code{node_names} ordering.
#'
#' @param node_names a character vector with the names of the (dynamic) nodes.
#' @param prob_edge_intraslice a single number in \code{[0, 1]}: probability of
#'   an intra-slice arc \eqn{X_t \to Y_t}.
#' @param prob_edge_interslice either a single number in \code{[0, 1]} (the same
#'   probability is used for every lag) or a numeric vector / list of length
#'   \code{markov_order}, where the \eqn{k}-th entry is the probability of a
#'   \eqn{t-k \to t} arc.
#' @param markov_order an integer (>= 1): the Markov order of the process.
#' @param g_0_arcs logical (default \code{TRUE}). If \code{TRUE} the prior
#'   network \eqn{G_0} is set to the same arcs as the intra-slice (every
#'   \eqn{X_t \to Y_t} arc is also added as \eqn{X_0 \to Y_0}). If \code{FALSE}
#'   an independent random \eqn{G_0} is generated using \code{g_0_prob}.
#' @param g_0_prob a single number in \code{[0, 1]} used only when
#'   \code{g_0_arcs = FALSE}: edge probability for the random \eqn{G_0} graph
#'   sampled via \code{bnlearn::random.graph}. Defaults to
#'   \code{prob_edge_intraslice}. If g_0_prob is set to 0 the g_0 network is empty.
#'
#' @return an object of class 'dbn'.
#' @export
#'
#' @examples
#' random.structure.dbn(c("A", "B", "C"), 0.5, 0.4, markov_order = 1)
#' random.structure.dbn(c("A", "B", "C"), 0.3, c(0.5, 0.2), markov_order = 2)
#' # independent prior network with its own edge density
#' random.structure.dbn(c("A", "B", "C"), 0.5, 0.4, markov_order = 1,
#'                       g_0_arcs = FALSE, g_0_prob = 0.2)
#' # empty g_0 network
#' random.structure.dbn(c("A", "B", "C"), 0.5, 0.4, markov_order = 1,
#'                       g_0_arcs = FALSE, g_0_prob = 0.0)

random.structure.dbn <- function(node_names,
                                 prob_edge_intraslice,
                                 prob_edge_interslice,
                                 markov_order,
                                 g_0_arcs = TRUE,
                                 g_0_prob = prob_edge_intraslice) {
  # ---- validation ----
  if (!is.character(node_names) || length(node_names) < 1) {
    stop("node_names must be a non-empty character vector of node names!")
  }
  if (length(markov_order) != 1 || markov_order < 1 ||
      markov_order != as.integer(markov_order)) {
    stop("markov_order must be a single integer >= 1.")
  }
  if (!is.numeric(prob_edge_intraslice) ||
      length(prob_edge_intraslice) != 1 ||
      prob_edge_intraslice < 0 || prob_edge_intraslice > 1) {
    stop("prob_edge_intraslice must be a single number in [0, 1].")
  }
  # prob_edge_interslice can be a float or a list/vector of per-lag probabilities
  prob_edge_interslice <- unlist(prob_edge_interslice)
  if (!is.numeric(prob_edge_interslice) ||
      any(prob_edge_interslice < 0) || any(prob_edge_interslice > 1)) {
    stop("prob_edge_interslice must be numeric with values in [0, 1].")
  }
  if (length(prob_edge_interslice) == 1) {
    prob_edge_interslice <- rep(prob_edge_interslice, markov_order)
  } else if (length(prob_edge_interslice) != markov_order) {
    stop("prob_edge_interslice must have length 1 or markov_order.")
  }
  if (!is.logical(g_0_arcs) || length(g_0_arcs) != 1 || is.na(g_0_arcs)) {
    stop("g_0_arcs must be a single logical value (TRUE/FALSE).")
  }
  if (!g_0_arcs &&
      (!is.numeric(g_0_prob) || length(g_0_prob) != 1 ||
       g_0_prob < 0 || g_0_prob > 1)) {
    stop("g_0_prob must be a single number in [0, 1].")
  }

  # ---- empty DBN scaffold ----
  generated_dbn <- empty.dbn(dynamic_nodes = node_names,
                             markov_order = markov_order)

  # ---- intra-slice arcs X_t -> Y_t (sampled as a DAG over the node ordering) ----
  intra_slice_edges <- character(0)
  if (prob_edge_intraslice > 0 && length(node_names) >= 2) {
    g_t_nodes_names <- paste0(node_names, "_t")
    intra_slice_edges <- as.vector(t(bnlearn::arcs(
      bnlearn::random.graph(g_t_nodes_names, prob = prob_edge_intraslice)
    )))
    for (i in seq_len(length(intra_slice_edges) / 2)) {
      node_from <- intra_slice_edges[2 * i - 1]
      node_to   <- intra_slice_edges[2 * i]
      generated_dbn <- add.arc.dbn(
        DBN  = generated_dbn,
        from = as.character(split_variable_name(node_from)),
        to   = as.character(split_variable_name(node_to))
      )
    }
  }

  # ---- prior network G_0 arcs ----
  if (g_0_arcs) {
    # mirror the slice-t structure into G_0: X_t -> Y_t  =>  X_0 -> Y_0
    for (i in seq_len(length(intra_slice_edges) / 2)) {
      from_name <- split_variable_name(intra_slice_edges[2 * i - 1])$name
      to_name   <- split_variable_name(intra_slice_edges[2 * i])$name
      generated_dbn <- add.arc.dbn(
        DBN  = generated_dbn,
        from = c(from_name, "t_0"),
        to   = c(to_name, "t_0")
      )
    }
  } else if (g_0_prob > 0 && length(node_names) >= 2) {
    # sample an independent random DAG for G_0
    g_0_nodes_names <- paste0(node_names, "_0")
    g_0_edges <- as.vector(t(bnlearn::arcs(
      bnlearn::random.graph(g_0_nodes_names, prob = g_0_prob)
    )))
    for (i in seq_len(length(g_0_edges) / 2)) {
      generated_dbn <- add.arc.dbn(
        DBN  = generated_dbn,
        from = as.character(split_variable_name(g_0_edges[2 * i - 1])),
        to   = as.character(split_variable_name(g_0_edges[2 * i]))
      )
    }
  }

  # ---- inter-slice (temporal) arcs X_t-k -> Y_t ----
  for (k in seq_len(markov_order)) {
    prob_k <- prob_edge_interslice[k]
    if (prob_k == 0) next
    for (from_var in node_names) {
      for (to_var in node_names) {
        if (rbinom(1, 1, prob_k) == 1) {
          generated_dbn <- add.arc.dbn(
            DBN  = generated_dbn,
            from = c(from_var, paste0("t-", k)),
            to   = c(to_var, "t")
          )
        }
      }
    }
  }

  return(generated_dbn)
}

