# -------------------------------
#           Unrolling
# -------------------------------

#' rename_nodes_unroll
#' @description
#' Renames a list of nodes, making nodes names compliant to the dbn unrolled structure.
#' A node at time \code{t} (\code{X_t}) is mapped to \code{X_{time_slice}} and a lagged
#' node \code{X_t-k} (any Markov lag \eqn{k \ge 1}) is mapped to \code{X_{time_slice - k}}.
#'
#' @param time_slice an integer representing a time slice
#' @param nodes_names a list of nodes names
#'
rename_nodes_unroll <- function(time_slice, nodes_names) {
  renamed_nodes <- c()
  for (name_n in nodes_names) {
    # _t -> time_slice, _t-k -> time_slice - k (works for any Markov lag k)
    root_node_name <- get_variable_name(name_n)
    lag <- get_variable_time_index(name_n)
    new_time <- time_slice - lag
    new_node_name <-
      paste(root_node_name, "_", as.character(new_time), sep = "")
    renamed_nodes <- c(renamed_nodes, new_node_name)
  }
  return(renamed_nodes)
}

#' get_node_edges
#' @description
#' This function gets all the edges of a node
#'
#' @param dbn_transition a dbn transition
#' @param node name of the node in the dbn transition
#' @param time_slice an integer representing a time slice
#'
#' @return all the edges of a specific node in a time slice
#'
get_node_edges <- function(dbn_transition, node, time_slice) {
  if (!is.bn.fit(dbn_transition)) {
    stop("dbn must be a bn.fit object!")
  }
  if (is.character(node) == FALSE) {
    stop("node must be a character!")
  }
  if (is.numeric(time_slice)==FALSE){
    stop("time slice must be an integer!")
  }
  node_edges <- c()
  root_node_name_n_t <- get_variable_name(node)
  new_node_name_n_t <-
    paste(root_node_name_n_t, "_", as.character(time_slice), sep = "")
  node_parents <- get_parent_set(dbn_transition, node)

  if (length(node_parents) > 0) {
    renamed_parents <- rename_nodes_unroll(time_slice, node_parents)

    for (n_renamed_parent in renamed_parents) {
      #using "from" "to" standard
      single_edge <- c(n_renamed_parent, new_node_name_n_t)
      node_edges <- c(node_edges, single_edge)
    }
  }
  return(node_edges)
}

# Type of a fitted node, from its class (and the fields it carries).
node_dist_kind <- function(node) {
  if (is.dnode(node)) "discrete"
  else if (is.cgnode(node)) "mixed"
  else if (is.gnode(node)) "gaussian"
  else stop("get_unrolled_dbn: unrecognized node type (no prob/regs/coefficients)")
}

# Reorder the columns of a conditional-gaussian coefficient matrix (and its sd
# vector) from the stored expand.grid(dlevels) order to the configuration order
# bnlearn expects (discrete parents sorted alphabetically by name, the first one
# varying fastest). custom.fit does not accept an explicit `dlevels`, so the
# columns must already be in bnlearn's order for the discrete-parent
# configurations to line up.
cg_columns_bnlearn_order <- function(coef, sd, dlevels) {
  dpar <- names(dlevels)
  if (length(dpar) <= 1) return(list(coef = coef, sd = unname(sd)))
  combos <- expand.grid(dlevels, stringsAsFactors = FALSE)  # row i <-> column i
  ord_parents <- dpar[order(dpar)]                          # alphabetical
  ord_idx <- do.call(order, c(lapply(rev(ord_parents),
                                     function(p) combos[[p]])))
  list(coef = coef[, ord_idx, drop = FALSE], sd = unname(sd[ord_idx]))
}

# Distribution object (in bnlearn::custom.fit format) for a G_0 node, copied
# as-is - G_0 parents already carry their unrolled names (X_0, X_1, ...).
node_dist <- function(node) {
  switch(node_dist_kind(node),
    discrete = node$prob,
    gaussian = list(coef = node$regs, sd = node$std),
    mixed    = {
      r <- cg_columns_bnlearn_order(node$coefficients, node$sd, node$dlevels)
      list(coef = r$coef, sd = r$sd)
    }
  )
}

# ------------- transition-node reductions, one per node type -------------
# Each returns list(dist = <custom.fit distribution>, parents = <renamed
# parent>). Lagged parents whose resolved slice would be negative are dropped.

# Discrete: marginalize (normalize) the CPT over the dropped parents.
reduce_discrete_transition <- function(cpt, node_t, time_slice) {
  dn <- names(dimnames(cpt))
  lags <- vapply(dn, get_variable_time_index, numeric(1))
  keep <- (time_slice - lags) >= 0

  if (all(keep)) {
    names(dimnames(cpt)) <- rename_nodes_unroll(time_slice, dn)
    reduced <- cpt
  } else {
    keep_names <- dn[keep]
    # average over the dropped parents -> the conditional stays normalized
    reduced <- apply(cpt, keep_names, mean)
    if (is.null(dim(reduced))) {
      reduced <- array(reduced, dim = length(reduced),
                       dimnames = setNames(list(names(reduced)), keep_names))
    }
    names(dimnames(reduced)) <- rename_nodes_unroll(time_slice, keep_names)
  }
  parent_dims <- dn[keep & (dn != node_t)]
  list(dist = reduced, parents = rename_nodes_unroll(time_slice, parent_dims))
}

# Gaussian: drop the regressor of each dropped parent (keep intercept + std).
reduce_gaussian_transition <- function(node, time_slice) {
  regs <- node$regs
  parent_names <- setdiff(names(regs), intercept_name)
  lags <- vapply(parent_names, get_variable_time_index, numeric(1))
  keep <- (time_slice - lags) >= 0
  kept <- parent_names[keep]

  new_regs <- c(regs[intercept_name], regs[kept])
  names(new_regs) <- c(intercept_name, rename_nodes_unroll(time_slice, kept))
  list(dist = list(coef = new_regs, sd = node$std),
       parents = rename_nodes_unroll(time_slice, kept))
}

# Mixed (conditional gaussian): drop a dropped gaussian parent's regressor row;
# drop a dropped discrete parent's dimension (fix it to its first level). The
# surviving discrete-parent columns are re-ordered into bnlearn's config order.
reduce_mixed_transition <- function(node, time_slice) {
  coef <- node$coefficients          # [reg_names x combos], combos = expand.grid(dlevels)
  sdv  <- node$sd
  dlev <- node$dlevels               # named list: discrete parent -> levels
  reg_names <- rownames(coef)        # (Intercept) + gaussian parents
  gpar <- setdiff(reg_names, intercept_name)
  dpar <- names(dlev)

  g_keep <- (time_slice - vapply(gpar, get_variable_time_index, numeric(1))) >= 0
  d_keep <- (time_slice - vapply(dpar, get_variable_time_index, numeric(1))) >= 0
  kept_g <- gpar[g_keep]
  kept_d <- dpar[d_keep]

  # combos (one row per coef column). Drop a discrete parent by keeping only the
  # columns where it is fixed to its first level.
  combos <- expand.grid(dlev, stringsAsFactors = FALSE)
  keep_col <- rep(TRUE, nrow(combos))
  for (dp in dpar[!d_keep]) keep_col <- keep_col & (combos[[dp]] == dlev[[dp]][1])

  # drop dropped gaussian parents' rows + dropped columns
  row_keep <- reg_names %in% c(intercept_name, kept_g)
  coef_sub <- coef[row_keep, keep_col, drop = FALSE]
  sd_sub   <- sdv[keep_col]
  rownames(coef_sub) <- c(intercept_name, rename_nodes_unroll(time_slice, kept_g))

  renamed_parents <- rename_nodes_unroll(time_slice, c(kept_g, kept_d))

  if (length(kept_d) == 0) {
    # no discrete parents left -> the node is now purely gaussian
    coef_vec <- coef_sub[, 1]
    names(coef_vec) <- rownames(coef_sub)
    return(list(dist = list(coef = coef_vec, sd = unname(sd_sub[1])),
                parents = renamed_parents))
  }

  # the surviving columns are in expand.grid(dlev_kept) order (kept discrete
  # parents in their original relative order); reorder them into bnlearn's order.
  dlev_kept <- dlev[kept_d]
  names(dlev_kept) <- rename_nodes_unroll(time_slice, kept_d)
  ord <- cg_columns_bnlearn_order(coef_sub, sd_sub, dlev_kept)
  coef_mat <- ord$coef
  colnames(coef_mat) <- as.character(seq_len(ncol(coef_mat)))
  list(dist = list(coef = coef_mat, sd = ord$sd), parents = renamed_parents)
}

# Unroll one transition node at `time_slice`, dispatching on its type.
# Returns list(node, dist, edges) where `edges` is a flat from/to vector.
unroll_transition_node <- function(node, node_t, time_slice) {
  new_node <- paste(get_variable_name(node_t), "_",
                    as.character(time_slice), sep = "")
  reduced <- switch(node_dist_kind(node),
    discrete = reduce_discrete_transition(node$prob, node_t, time_slice),
    gaussian = reduce_gaussian_transition(node, time_slice),
    mixed    = reduce_mixed_transition(node, time_slice)
  )
  edges <- c()
  for (rp in reduced$parents) edges <- c(edges, rp, new_node)
  list(node = new_node, dist = reduced$dist, edges = edges)
}

#' get_unrolled_dbn
#'
#' @description this function generates an unrolled dynamic bayesian network.
#'
#' @details
#' The initial time slices are taken from the prior network \eqn{G_0}. When the
#' DBN was built with an extended \eqn{G_0} (\code{extend_g0 = TRUE}). 
#' When \eqn{G_0} only defines the slice 0 (\code{X_0}), a single initial 
#' slice is taken. The remaining slices, up to \code{slices}, are generated by 
#' unrolling the transition network.
#'
#' \code{extend_with_transition} controls how the initial Markov-order slices are
#' produced:
#' \itemize{
#'   \item \code{FALSE} (default): if \eqn{G_0} defines the extended slices
#'     they are used as-is; otherwise those slices are generated from the 
#'     transition network.
#'   \item \code{TRUE}: the extended \eqn{G_0} slices are ignored and every slice
#'     after slice 0 is generated from the transition network.
#' }
#' 
#' When a slice is generated from the transition network but a lagged parent
#' \code{X_t-k} would resolve to a negative slice that parent is dropped from 
#' the node's distribution according to its type:
#' \itemize{
#'   \item discrete node: the CPT is marginalized (normalized) over the parent;
#'   \item gaussian node: the parent's regressor is dropped;
#'   \item mixed node: a gaussian parent's regressor is dropped, a discrete
#'     parent's dimension is dropped (fixed to its first level).
#' }
#'
#' If \code{slices < markov_order} only the \eqn{G_0} slices up to \code{slices}
#' are returned and the transition network is not unrolled.
#'
#' @param dbn_fitted a dbn.fit class object
#' @param slices highest time index of the unrolled network (the network spans
#'   slices \code{0..slices})
#' @param extend_with_transition logical (default \code{FALSE}). If \code{TRUE},
#'   ignore the extended \eqn{G_0} slices and always rebuild slices
#'   \code{1..markov_order-1} (and beyond) from the transition network.
#'
#' @return a bn.fit object
#' @export
#'
#' @examples get_unrolled_dbn(my_fitted_dbn, 4)
get_unrolled_dbn <- function(dbn_fitted, slices, extend_with_transition = FALSE) {
  if (is.character(slices)) {
    stop("slices must be an integer!")
  }
  if (slices < 1){
    stop("slices must be greater than 0!")
  }
  if (!is.dbn.fit(dbn_fitted)) {
    stop("dbn_fitted must be a dbn.fit object")
  }
  if (!is.logical(extend_with_transition) ||
      length(extend_with_transition) != 1) {
    stop("extend_with_transition must be a single logical value!")
  }
  nodes_bn <- c()
  edges_bn <- c()
  dist_bn <- list()
  my_dbn_transition <- get.transition.net(dbn_fitted)
  # drop lagged (t-k) parent only the _t transition nodes are returned.
  nodes_in_transition <-
    get_nodes_t(remove_prev_time_from_bn_fit(my_dbn_transition))
  my_dbn_g_0 <- get.g0.net(dbn_fitted)
  nodes_in_0 <- bnlearn::node.ordering(my_dbn_g_0)

  # G_0 may cover several initial slices (X_0, X_1, ..., X_{g0_max}) when the DBN
  # was created with extend_g0 = TRUE. G_0 nodes carry a plain integer slice
  # suffix (X_0, X_1, ...).
  g0_slice_idx <- as.integer(sub(".*_([0-9]+)$", "\\1", nodes_in_0))
  g0_max <- if (length(g0_slice_idx) > 0) max(g0_slice_idx) else 0

  # Number of initial slices actually taken from G_0:
  #  - extend_with_transition = FALSE: use every slice G_0 provides (0..g0_max);
  #    when G_0 is not extended this is just slice 0.
  #  - extend_with_transition = TRUE : ignore the extended G_0 slices and take
  #    only slice 0, rebuilding the rest from the transition network.
  g0_use <- if (extend_with_transition) 0 else g0_max

  # --- initial slices: copy the G_0 nodes (X_0 .. X_{min(g0_use, slices)}) ---
  for (n_0 in nodes_in_0) {
    slice_n_0 <- as.integer(sub(".*_([0-9]+)$", "\\1", n_0))
    if (slice_n_0 > g0_use || slice_n_0 > slices) next
    #add the node in bn_nodes
    nodes_bn <- c(nodes_bn, n_0)
    dist_bn[[n_0]] <- node_dist(my_dbn_g_0[[n_0]])
    parents_n_0 <- my_dbn_g_0[[n_0]][["parents"]]
    if (length(parents_n_0) > 0) {
      for (p_n_0 in parents_n_0) {
        #using "from" "to" standard
        edges_bn <- c(edges_bn, p_n_0, n_0)
      }
    }
  }

  # --- remaining slices: unroll the transition net right after the G_0 slices ---
  if (slices > g0_use) {
    for (time_slice in (g0_use + 1):slices) {
      for (n_t in nodes_in_transition) {
        unrolled <- unroll_transition_node(
          my_dbn_transition[[n_t]], n_t, time_slice
        )
        dist_bn[[unrolled$node]] <- unrolled$dist
        #adding node
        nodes_bn <- c(nodes_bn, unrolled$node)
        #adding edges
        if (length(unrolled$edges) > 0) {
          edges_bn <- c(edges_bn, unrolled$edges)
        }
      }
    }
  }

  dag_unrolled <- bnlearn::empty.graph(nodes = nodes_bn)
  if (length(edges_bn) > 0) {
    arc_unrolled.set <- matrix(
      edges_bn,
      byrow = TRUE,
      ncol = 2,
      dimnames = list(NULL, c("from", "to"))
    )
    bnlearn::arcs(dag_unrolled) <- arc_unrolled.set
  }
  dbn_unrolled <- bnlearn::custom.fit(dag_unrolled, dist_bn)
  return(dbn_unrolled)
}
