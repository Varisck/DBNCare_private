# Private helpers shared across DBN initialization, arc editing, and conversion.
# Not exported.

node_id <- function(name, time) {
  if (time == 't_0') paste0(name, '_0') else paste0(name, '_', time)
}

parse_node_id <- function(id) {
  parts <- strsplit(id, '_')[[1]]
  name <- paste(parts[1:(length(parts) - 1)], collapse = '_')
  time_token <- parts[length(parts)]
  time <- ifelse(time_token == '0', 't_0', time_token)
  list(name = name, time = time)
}

is_valid_parent_time <- function(s) {
  if (s == 't') return(TRUE)
  if (nchar(s) <= 2) return(FALSE)
  if (substring(s, 1, 2) != 't-') return(FALSE)
  !is.na(as.numeric(substring(s, 3)))
}

empty_node_details <- function() {
  list(
    mb = character(0),
    nbr = character(0),
    parents = character(0),
    children = character(0)
  )
}

empty_bn_shell <- function() {
  list(
    learning = list(
      whitelist = NULL,
      blacklist = NULL,
      test = 'none',
      ntests = 0,
      algo = 'empty',
      args = list()
    ),
    arcs = matrix(
      ncol = 2,
      nrow = 0,
      byrow = TRUE,
      dimnames = list(character(0), c("from", "to"))
    ),
    nodes = list()
  )
}

delete_arc_element <- function(x, element) {
  if (class(x)[1] == 'character') {
    x[!x == element]
  } else {
    if (dim(x)[1] != 2) {
      x[-prodlim::row.match(element, x), ]
    } else {
      matrix(
        x[-prodlim::row.match(element, x), ],
        nrow = 1,
        ncol = 2,
        dimnames = list(NULL, c('from', 'to'))
      )
    }
  }
}

#' create blacklist form list of nodes in dbn
#'
#' @param g_0_nodes all nodes in the g_0 network like: "A_0", ...
#' @param g_t_nodes all nodes ending in t like: "A_t", ...
#' @param markov_order markov_order to generate blacklist from
#' @param static_nodes character vector of time-independent node names (without suffix), e.g. c("S", "K")
#' @param allow_intraslice_edges Default TRUE weather or not to allow intraslice edges (A_t -> B_t)
#' @param allow_t_0_edges Default TRUE weather or not to allow edges at t = 0 (A_0 -> B_0)
#'
#' @return matrix of node pairs which edge is forbidden
#' @export
#'
#' @examples
#' blacklist_g0_gt(c("A_0", "B_0"), c("A_t", "B_t"), markov_order = 3)
#' blacklist_g0_gt(c("A_0", "B_0"), c("A_t", "B_t"), static_nodes = c("S"))
#'
blacklist_g0_gt = function(g_0_nodes, g_t_nodes, markov_order = 1,
                           static_nodes = c(),
                           allow_intraslice_edges = TRUE,
                           allow_t_0_edges = TRUE) {
  b = rbind(
    # A_0 -> A_t
    expand.grid(from = g_0_nodes, to = g_t_nodes),
    # A_t -> A_0
    expand.grid(from = g_t_nodes, to = g_0_nodes)
  ) |> as.matrix()

  for(mo in 1:markov_order) {
    g_t_i_nodes = sapply(g_t_nodes, concat_name_post, postfix = paste("-", mo))
    b = rbind(
      b,
      # A_t-i -> A_0
      expand.grid(from = g_t_i_nodes, to = g_0_nodes),
      # A_0 -> A_t-i
      expand.grid(from = g_0_nodes, to = g_t_i_nodes),
      # A_t -> A_t-i
      expand.grid(from = g_t_nodes, to = g_t_i_nodes)
    )

    for(mo2 in 1:markov_order) {
      g_t_j_nodes = sapply(g_t_nodes, concat_name_post, postfix = paste("-", mo2))
      b = rbind(
        b,
        # A_t-i -> A_t-j
        expand.grid(from = g_t_i_nodes, to = g_t_j_nodes)
      )
    }
  }

  if (length(static_nodes) > 0) {
    b = rbind(
      b,
      # dynamic_0 -> static_0 
      expand.grid(from = g_0_nodes, to = static_nodes),
      # dynamic_t -> static_0
      expand.grid(from = g_t_nodes, to = static_nodes),
      # static_0 -> dynamic_t
      expand.grid(from = static_nodes, to = g_t_nodes)
    ) |> as.matrix()
    for (mo in 1:markov_order) {
      g_t_i_nodes = sapply(g_t_nodes, concat_name_post, postfix = paste("-", mo))
      b = rbind(
        b,
        # dynamic_t-i -> static_0
        expand.grid(from = g_t_i_nodes, to = static_nodes),
        # static_0 -> dynamic_t-i
        expand.grid(from = static_nodes, to = g_t_i_nodes)
      )
    }
  }

  if(!allow_intraslice_edges)
    b = rbind(
      b,
      expand.grid(from = g_t_nodes, to = g_t_nodes)
    )
  if(!allow_t_0_edges)
    b = rbind(
      b,
      expand.grid(from = g_0_nodes, to = g_0_nodes)
    )
  return(b)
}



# DBN SUMMARY FUNCTION


#' Print a Dynamic Bayesian Network
#'
#' @param DBN object of class 'dbn'
#' @param ... further arguments passed to or from other methods
#'
#' @return the input \code{DBN} object, invisibly
#' @export
#'
#' @examples
#' DBN_example <- empty.dbn(static_nodes = c(), dynamic_nodes = c("A", "B", "C"), markov_order = 1)
#' print(DBN_example)
#' DBN_example
print.dbn <- function(DBN, ...) {
  summary(DBN)
  invisible(DBN)
}

#' Function returning a summary of input Dynamic Bayesian Network
#'
#' @param DBN object of class 'dbn'
#'
#' @details
#' The function returns from terminal:
#' - the model string corresponding to the DAG at G_0
#' - the model string corresponding to the DAG at G_transition
#' - the Markovian Order of the process
#' - the number of nodes, divided between Dynamic nodes and Static nodes
#' - the number of arcs, divided between arcs of G_0, inner-slice arcs and intra-slice arcs of G_transition
#' - the average Markov Blanket size of the nodes
#' - the average neighboorhood size of the nodes
#'
#' @export
#'
#' @examples
#' summary(DBN_example)
summary.dbn <- function(DBN, test = FALSE) {
  if (!class(DBN) == 'dbn')
    stop("ERROR: DBN argument is not of class 'dbn'")
  mo = DBN$markov_order
  nodes <- length(DBN$nodes)
  dynamic_nodes = 0
  mbs <- c()
  neighbs <- c()
  for (i in names(DBN$nodes)) {
    dynamic_nodes = dynamic_nodes + ifelse(DBN$nodes[[i]][['type']] == 'Dynamic', 1, 0)
    mbs <-
      c(mbs, length(DBN$nodes[[i]][['t_0']][['mb']]), length(DBN$nodes[[i]][['t']][['mb']]))
    neighbs <-
      c(neighbs, length(DBN$nodes[[i]][['t_0']][['nbr']]), length(DBN$nodes[[i]][['t']][['nbr']]))
  }
  static_nodes <- nodes - dynamic_nodes
  arcs = nrow(DBN$arcs)
  intra_slice_arcs = 0
  inner_slice_arcs = length(unlist(sapply(DBN$arcs[, 'from'], function(x) {
    quanteda::char_select(x, '*_t')
  })))
  prior_arcs = length(unlist(sapply(DBN$arcs[, 'to'], function(x) {
    quanteda::char_select(x, '*_0')
  })))
  for (j in 1:mo) {
    intra_slice_arcs <-
      intra_slice_arcs + length(unlist(sapply(DBN$arcs[, 'from'], function(x) {
        quanteda::char_select(x, paste0('*_t-', j))
      })))
  }
  if (test == FALSE) {
    cat('Dynamic Bayesian Network \n\n\n')
    if (bnlearn::directed(from_DBN_to_G_0(DBN))) {
      cat('Prior Network Model\n ',
          bnlearn::modelstring(from_DBN_to_G_0(DBN)),
          '\n\n')
    }
    else{
      cat('Prior Network Model\n  [partially directed graph]\n\n')
    }
    if (bnlearn::directed(from_DBN_to_G_transition(DBN))) {
      cat(
        'Trantision Network Model\n ',
        bnlearn::modelstring(from_DBN_to_G_transition(DBN)),
        '\n\n'
      )
    }
    else{
      cat('Trantision Network Model\n  [partially directed graph]\n\n')
    }
    cat("Markovian Order:           \t\t", mo, '\n\n')
    cat('Nodes:                     \t\t', nodes, '\n')
    cat('\tDynamic Nodes:           \t', dynamic_nodes, '\n')
    cat('\tStatic Nodes:            \t', static_nodes, '\n\n')
    cat('Arcs:                      \t\t', nrow(DBN$arcs), '\n')
    cat('\tPrior Arcs:              \t', prior_arcs, '\n')
    cat('\tInner-slice Arcs:        \t', inner_slice_arcs, '\n')
    cat('\tIntra-slice Arcs:        \t', intra_slice_arcs, '\n\n')
    cat('Avg Markov Blanket size:   \t\t', round(mean(mbs), digits = 4), '\n')
    cat('Avg Neighboorhood size:    \t\t', round(mean(neighbs), digits = 4), '\n')
  }
  result_list <-
    list(
      markov_order = mo,
      n_nodes = nodes,
      n_dynamic_nodes = dynamic_nodes,
      n_static_nodes = static_nodes,
      n_arcs = nrow(DBN$arcs),
      n_prior_arcs = prior_arcs,
      n_inner_slice_arcs = inner_slice_arcs,
      n_intra_slice_arcs = intra_slice_arcs,
      avg_mb = round(mean(mbs), digits = 4),
      avg_nbr = round(mean(neighbs), digits = 4)
    )
  invisible(result_list)
}


#' Print a fitted Dynamic Bayesian Network
#'
#' @param x object of class 'dbn.fit'
#' @param ... further arguments passed to or from other methods
#'
#' @return the input \code{x} object, invisibly
#' @method print dbn.fit
#' @export
#'
#' @examples
#' print(fitted_dbn)
#' fitted_dbn
print.dbn.fit <- function(x, ...) {
  summary(x)
  invisible(x)
}


#' Summary of a fitted Dynamic Bayesian Network
#'
#' @param object object of class 'dbn.fit'
#' @param ... further arguments passed to or from other methods
#'
#' @details
#' Prints a per-node overview in the style of \pkg{bnlearn}'s
#' \code{summary.bn.fit}:
#' \itemize{
#'   \item distribution type (discrete / gaussian / mixed)
#'   \item Markov order
#'   \item for each prior node (\code{_0}): parents and distribution info
#'   \item for each transition node (\code{_t}): parents and distribution info
#' }
#' For discrete nodes the number of levels is shown; for Gaussian nodes the
#' coefficient names and the residual standard deviation are shown.
#'
#' @return a list with elements \code{type}, \code{markov_order},
#'   \code{n_prior_nodes}, \code{n_transition_nodes}, and \code{nodes}
#'   (a named list of per-node summaries), invisibly.
#' @method summary dbn.fit
#' @export
#'
#' @examples
#' summary(fitted_dbn)
summary.dbn.fit <- function(object, ...) {
  if (!inherits(object, "dbn.fit"))
    stop("ERROR: object argument is not of class 'dbn.fit'")

  type <- dbn_type(object)
  mo   <- get_max_mo_dbn_fit(object)

  prior_nodes      <- sort(names(object)[grepl("_0$",  names(object))])
  transition_nodes <- sort(names(object)[grepl("_t$",  names(object))])

  .node_info <- function(n) {
    nd <- object[[n]]
    parents  <- if (length(nd$parents)  == 0) "none" else paste(nd$parents,  collapse = " ")
    children <- if (length(nd$children) == 0) "none" else paste(nd$children, collapse = " ")
    if (!is.null(nd$prob)) {
      lvls <- length(dimnames(nd$prob)[[1]])
      dist_str <- paste0("[", lvls, " levels]")
    } else {
      coef_str <- paste(names(nd$regs), collapse = " ")
      dist_str <- paste0("[coefficients: ", coef_str, " | sd: ", round(nd$std, 4), "]")
    }
    list(parents = parents, children = children, dist = dist_str)
  }

  cat("Fitted Dynamic Bayesian Network\n\n")
  cat("  Type:", type, "\n")
  cat("  Markov order:", mo, "\n\n")

  cat("  Prior Network (", length(prior_nodes), "nodes )\n", sep = "")
  for (n in prior_nodes) {
    info <- .node_info(n)
    cat("    -", n, info$dist, "\n")
    cat("       Parents:", info$parents, "\n")
  }

  cat("\n  Transition Network (", length(transition_nodes), "nodes )\n", sep = "")
  for (n in transition_nodes) {
    info <- .node_info(n)
    cat("    -", n, info$dist, "\n")
    cat("       Parents:", info$parents, "\n")
  }
  cat("\n")

  node_summaries <- setNames(
    lapply(c(prior_nodes, transition_nodes), .node_info),
    c(prior_nodes, transition_nodes)
  )
  invisible(list(
    type               = type,
    markov_order       = mo,
    n_prior_nodes      = length(prior_nodes),
    n_transition_nodes = length(transition_nodes),
    nodes              = node_summaries
  ))
}



#' Model string of a network
#'
#' @description
#' S3 generic returning the \pkg{bnlearn} model string of a network, dispatching
#' on the class of \code{x} just like \code{print}. For a \code{dbn} or
#' \code{dbn.fit} it returns a named list with the prior-network (\code{g_0}) and
#' transition-network (\code{g_t}) model strings, in the same format printed by
#' [summary.dbn()]. For any other object (e.g. \code{bn}, \code{bn.fit}) it
#' delegates to [bnlearn::modelstring()], returning a single string.
#'
#' @param x a network object (\code{dbn}, \code{dbn.fit}, \code{bn} or \code{bn.fit}).
#' @param ... further arguments passed to the method / [bnlearn::modelstring()].
#'
#' @return for \code{dbn} / \code{dbn.fit} a named list (\code{g_0}, \code{g_t});
#'   otherwise a single character string.
#' @export
modelstring <- function(x, ...) UseMethod("modelstring")

#' @rdname modelstring
#' @method modelstring default
#' @export
modelstring.default <- function(x, ...) bnlearn::modelstring(x, ...)


#' Model strings of a fitted Dynamic Bayesian Network
#'
#' @description
#' Returns the \pkg{bnlearn} model strings of the two networks of a fitted DBN:
#' the prior network \eqn{G_0} and the transition network \eqn{G_t}. The strings
#' use the same format printed by [summary.dbn()] (e.g.
#' \code{"[A_0][B_0][C_0|A_0:B_0]"}); the transition string also lists the lagged
#' root nodes \eqn{X_{t-k}}, which in a fitted DBN only appear as parent
#' references.
#'
#' @param dbn an object of class 'dbn.fit'
#'
#' @return a named list with two character strings: \code{g_0} (prior network
#'   model string) and \code{g_t} (transition network model string).
#' @method modelstring dbn.fit
#' @export
#' 
#'
#' @examples
#' \dontrun{
#' fit <- fit_random_dbn(random.structure.dbn(c("A","B","C"), .5, .5, 1), "continuous")
#' ms  <- modelstring.dbn.fit(fit)
#' ms$g_0
#' ms$g_t
#' }
modelstring.dbn.fit <- function(dbn, ...) {
  if (!inherits(dbn, "dbn.fit"))
    stop("ERROR: dbn argument is not of class 'dbn.fit'")

  # Build a bn over `all_nodes` whose arcs are parent -> child for every child,
  # then return its bnlearn model string (same format as summary.dbn).
  slice_modelstring <- function(child_nodes, all_nodes) {
    g <- bnlearn::empty.graph(all_nodes)
    arc_mat <- do.call(rbind, lapply(child_nodes, function(n) {
      parents <- dbn[[n]][["parents"]]
      if (length(parents)) cbind(from = parents, to = n) else NULL
    }))
    if (!is.null(arc_mat)) bnlearn::arcs(g) <- arc_mat
    bnlearn::modelstring(g)
  }

  # G_0: prior-network nodes (X_0); their parents are X_0 nodes too
  g_0_nodes <- grep("_0$", names(dbn), value = TRUE)

  # G_t: transition nodes (X_t) plus the lagged root nodes X_t-k (which only
  # appear as parent references in a fitted DBN) so the string matches the
  # transition network printed by summary.dbn
  g_t_nodes    <- grep("_t$", names(dbn), value = TRUE)
  variables    <- unique(vapply(g_t_nodes, get_variable_name, character(1)))
  markov_order <- get_max_mo_dbn_fit(dbn)
  lagged_nodes <- character(0)
  for (k in seq_len(markov_order))
    lagged_nodes <- c(lagged_nodes, paste0(variables, "_t-", k))

  list(
    g_0 = slice_modelstring(g_0_nodes, g_0_nodes),
    g_t = slice_modelstring(g_t_nodes, c(paste0(variables, "_t"), lagged_nodes))
  )
}

#' Model strings of a Dynamic Bayesian Network
#'
#'
#' @param dbn an object of class 'dbn'
#'
#' @return a named list with two character strings: \code{g_0} (prior network
#'   model string) and \code{g_t} (transition network model string).
#' @rdname modelstring
#' @method modelstring dbn
#' @export
#'
#' @examples
#' \dontrun{
#' dbn <- random.structure.dbn(c("A","B","C"), .5, .5, 1)
#' ms  <- modelstring(dbn)
#' ms$g_0
#' ms$g_t
#' }
modelstring.dbn <- function(dbn, ...) {
  g0 <- get.g0.net(dbn)
  gt <- get.transition.net(dbn)
  list(
    g_0 = if (bnlearn::directed(g0)) bnlearn::modelstring(g0) else "[partially directed graph]",
    g_t = if (bnlearn::directed(gt)) bnlearn::modelstring(gt) else "[partially directed graph]"
  )
}