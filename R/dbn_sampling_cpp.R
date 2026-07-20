# Rcpp-backed implementation of dbn.sampling().
#
# The sampling loops live in src/dbn_sampling.cpp; this file translates the
# dbn.fit object into the integer-indexed "plans" consumed by the C++ cores
# and re-assembles the resulting matrix into the same data.frame returned by
# dbn.sampling(). The C++ cores consume the R RNG stream with the same draws
# in the same order as the pure-R implementation, so for a fixed seed
# dbn.sampling() and dbn.sampling.cpp() return the same dataset.


# resolve a full variable name (e.g. "A_t-2") to the 0-based column index of
# its base variable and its time lag (0 for "A_0"/"A_t", k for "A_t-k")
plan_parent_ref <- function(parent, var_index) {
  name <- split_variable_name(parent)$name
  if (is.na(var_index[name])) {
    stop(paste("dbn.sampling.cpp: parent", parent,
               "does not match any variable of the network"))
  }
  list(var = var_index[[name]], lag = get_variable_time_index(parent))
}

# plan of a gaussian node: regression coefficients in parents order, indexed
# like sample_variable_gaussian() does (regs[1] = intercept, regs[i + 1] =
# coefficient of parents[i]).
#
# `g0 = TRUE` marks a prior-network node: `time` is the slice it writes to
# (0 for var_0, 1 for var_1, ...) and each parent lag is stored relative to
# that slice, so the C++ core reads the parent at block + (time - lag). For a
# transition node (g0 = FALSE) time = 0 and the lag stays the absolute t-k lag.
plan_gaussian_node <- function(net, variable, var_index, g0 = FALSE) {
  node <- net[[variable]]
  if (is.null(node$regs)) {
    stop(paste("dbn.sampling.cpp: node", variable,
               "has no regression coefficients (not a gaussian node)"))
  }
  parents <- node$parents
  regs <- as.numeric(node$regs)
  if (length(regs) != length(parents) + 1L) {
    stop(paste("dbn.sampling.cpp: node", variable, "has", length(regs),
               "coefficients for", length(parents), "parents"))
  }
  node_time <- if (g0) get_variable_time_index(variable) else 0L
  par_var <- integer(length(parents))
  par_lag <- integer(length(parents))
  for (i in seq_along(parents)) {
    ref <- plan_parent_ref(parents[i], var_index)
    par_var[i] <- ref$var
    par_lag[i] <- if (g0) node_time - ref$lag else ref$lag
  }

  list(var = var_index[[split_variable_name(variable)$name]],
       time = as.integer(node_time),
       intercept = regs[1],
       std = as.numeric(node$std),
       par_var = par_var,
       par_lag = par_lag,
       par_coef = regs[-1])
}

# canonical level labels of each base variable: the levels of its _t CPT,
# extended with any extra level of its _0 CPT. Values are stored as 0-based
# codes into these vectors while sampling.
discrete_levels <- function(bn_0, bn_transition, base_names) {
  levels_list <- list()
  for (b in base_names) {
    l_0 <- dimnames(bn_0[[paste0(b, "_0")]]$prob)[[1]]
    l_t <- dimnames(bn_transition[[paste0(b, "_t")]]$prob)[[1]]
    if (is.null(l_0) || is.null(l_t)) {
      stop(paste("dbn.sampling.cpp: missing CPT levels for variable", b))
    }
    levels_list[[b]] <- union(l_t, l_0)
  }
  levels_list
}

# plan of a discrete node: the CPT as its flat probability vector plus, per
# dimension, the translation of the level labels to canonical codes. The CPT
# layout (dim 1 = the node itself, then its parents) matches the row order of
# data.frame(as.table(prob)) used by sample_variable_discrete().
# `g0` behaves as in plan_gaussian_node(): prior nodes carry the slice they
# write to (`time`) and store parent lags relative to that slice.
plan_discrete_node <- function(net, variable, var_index, levels_list, g0 = FALSE) {
  node <- net[[variable]]
  prob <- node$prob
  if (is.null(prob)) {
    stop(paste("dbn.sampling.cpp: node", variable,
               "has no CPT (not a discrete node)"))
  }
  dn <- dimnames(prob)
  dnn <- names(dn)
  if (is.null(dnn) || dnn[1] != variable) {
    stop(paste("dbn.sampling.cpp: the first dimension of the CPT of",
               variable, "must be the node itself"))
  }
  if (!setequal(dnn[-1], node$parents)) {
    stop(paste("dbn.sampling.cpp: CPT dimensions of", variable,
               "do not match its parent set"))
  }
  base <- split_variable_name(variable)$name
  own_map <- match(dn[[1]], levels_list[[base]]) - 1L

  node_time <- if (g0) get_variable_time_index(variable) else 0L
  n_parents <- length(dnn) - 1L
  par_var <- integer(n_parents)
  par_lag <- integer(n_parents)
  par_map <- vector("list", n_parents)
  for (j in seq_len(n_parents)) {
    parent <- dnn[j + 1L]
    ref <- plan_parent_ref(parent, var_index)
    par_var[j] <- ref$var
    par_lag[j] <- if (g0) node_time - ref$lag else ref$lag
    codes <- match(dn[[j + 1L]], levels_list[[split_variable_name(parent)$name]]) - 1L
    codes[is.na(codes)] <- -1L
    par_map[[j]] <- as.integer(codes)
  }
  list(var = var_index[[base]],
       time = as.integer(node_time),
       dims = as.integer(dim(prob)),
       own_map = as.integer(own_map),
       freq = as.numeric(prob),
       par_var = par_var,
       par_lag = par_lag,
       par_map = par_map)
}


# `g0` behaves as in plan_gaussian_node(): prior nodes carry the slice they
# write to (`time`) and store both continuous and discrete parent lags relative
# to that slice.
plan_mixed_node = function(net, variable, var_index, levels_list, g0 = FALSE) {
  node <- net[[variable]]
  if (is.null(node$coefficients)) {
    stop(paste("dbn.sampling.cpp: node", variable,
               "has no conditional gaussian coefficients (not a cgnode)"))
  }

  coef_mat <- node$coefficients
  reg_names <- dimnames(coef_mat)[[1]]
  if (is.null(reg_names) || reg_names[1] != intercept_name) {
    stop(paste("dbn.sampling.cpp: the first regressor of node", variable,
               "must be", intercept_name))
  }
  n_combos <- ncol(coef_mat)
  if (length(node$sd) != n_combos) {
    stop(paste("dbn.sampling.cpp: node", variable, "has", length(node$sd),
               "residual sds for", n_combos, "discrete-parent combinations"))
  }

  node_time <- if (g0) get_variable_time_index(variable) else 0L

  # continuous part: like plan_gaussian_node, one (var, lag) per continuous
  # parent. The intercept (row 1) is kept separate, mirroring the gaussian plan.
  cont_names <- reg_names[-1]
  cont_par_var <- integer(length(cont_names))
  cont_par_lag <- integer(length(cont_names))
  for (i in seq_along(cont_names)) {
    ref <- plan_parent_ref(cont_names[i], var_index)
    cont_par_var[i] <- ref$var
    cont_par_lag[i] <- if (g0) node_time - ref$lag else ref$lag
  }

  # coefficients are one vector per discrete-parent combination
  par_coef <- lapply(seq_len(n_combos),
                     function(k) as.numeric(coef_mat[cont_names, k]))

  # discrete part: like the par_map of plan_discrete_node, one level->code
  # translation per discrete parent, kept in $dlevels (= expand.grid) order so
  # the combo index lines up with the columns of $coefficients / $sd.
  dlevels <- node$dlevels
  disc_names <- names(dlevels)
  disc_dims <- integer(length(dlevels))
  dis_par_var <- integer(length(dlevels))
  dis_par_lag <- integer(length(dlevels))
  par_code <- vector("list", length(dlevels))
  for (j in seq_along(dlevels)) {
    ref <- plan_parent_ref(disc_names[j], var_index)
    dis_par_var[j] <- ref$var
    dis_par_lag[j] <- if (g0) node_time - ref$lag else ref$lag
    disc_dims[j] <- length(dlevels[[j]])
    codes <- match(dlevels[[j]],
                   levels_list[[split_variable_name(disc_names[j])$name]]) - 1L
    codes[is.na(codes)] <- -1L
    par_code[[j]] <- as.integer(codes)
  }

  list(var = var_index[[split_variable_name(variable)$name]],
       time = as.integer(node_time),
       intercept = as.numeric(coef_mat[1, ]),
       std = as.numeric(node$sd),
       cont_par_var = cont_par_var,
       cont_par_lag = cont_par_lag,
       par_coef = par_coef,
       disc_dims = disc_dims,
       disc_par_var = dis_par_var,
       disc_par_lag = dis_par_lag,
       disc_par_map = par_code)
}

plan_all_mix = function(net, ordered_nodes, var_index, levels_list, g0 = FALSE) {
  plan = list()
  plan[["discrete"]] = list()
  plan[["gaussian"]] = list()
  plan[["mixed"]] = list()

  disc_idx <- 0L
  gauss_idx <- 0L
  mix_idx <- 0L

  type_vec <- integer(length(ordered_nodes))
  idx_vec <- integer(length(ordered_nodes))

  for (i in seq_along(ordered_nodes)) {
    nname = ordered_nodes[i]
    node = net[[nname]]
    if (inherits(node, "bn.fit.dnode")) {
      plan[["discrete"]][[nname]] = plan_discrete_node(net, nname, var_index, levels_list, g0 = g0)
      type_vec[i] = 0L
      idx_vec[i] = disc_idx
      disc_idx = disc_idx + 1
    } else if (inherits(node, "bn.fit.gnode")) {
      plan[["gaussian"]][[nname]] = plan_gaussian_node(net, nname, var_index, g0 = g0)
      type_vec[i] = 1L
      idx_vec[i] = gauss_idx
      gauss_idx = gauss_idx + 1
    } else if (inherits(node, "bn.fit.cgnode")) {
      plan[["mixed"]][[nname]] = plan_mixed_node(net, nname, var_index, levels_list, g0 = g0)
      type_vec[i] = 2L
      idx_vec[i] = mix_idx
      mix_idx = mix_idx + 1
    } else stop(paste("ERROR: class not recognized for node", nname, "after get.transition.net"))
  }
  list(plan = plan, 
       ordering = list(vector_type = type_vec, 
                       index_vec = idx_vec))
}


#' Generate a sampling dataset
#'
#'
#' @param fitted_dbn an object of class 'dbn.fit'
#' @param n_samples number of samples
#' @param max_time time series length
#' 
#' @returns the generated dataframe
#' @useDynLib DBNCare, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
#' 
#' @examples
#' dbn.sampling(DBN_example, N_samples, Time)
dbn.sampling <- function(fitted_dbn, n_samples, max_time) {
  if (is.character(max_time)) {
    stop("Time must be an integer!")
  }
  if (max_time < 1) {
    stop("Time must be greater than 0!")
  }
  if (is.character(n_samples)) {
    stop("N_samples must be an integer!")
  }
  if (n_samples < 1) {
    stop("N_samples must be greater than 0!")
  }
  if (!is.dbn.fit(fitted_dbn)) {
    stop("fitted_DBN must be a dbn.fit object")
  }

  dbn_type <- dbn_type(fitted_dbn)
  if (!dbn_type %in% c("discrete", "gaussian", "mixed")) {
    stop("Invalid dbn_type")
  }

  bn_0 <- get.g0.net(fitted_dbn)
  bn_transition <- get.transition.net(fitted_dbn)
  # sampling order: node ordering of G_0 at t = 0, then node ordering of the
  # transition network restricted to the _t nodes for t = 1..max_time
  nodes_0 <- bnlearn::node.ordering(bn_0)
  nodes_t <- get_nodes_t(remove_prev_time_from_bn_fit(bn_transition))

  # output columns follow the order in which the variables first enter the
  # time-series dictionary in dbn.sampling(), i.e. the t = 0 node ordering.
  # unique(): with an extended G_0 the same base variable appears once per
  # initial slice (var_0, var_1, ...) but maps to a single output column.
  base_names <- unique(vapply(nodes_0, function(v) split_variable_name(v)$name,
                              character(1), USE.NAMES = FALSE))
  for (v in nodes_t) {
    if (!split_variable_name(v)$name %in% base_names) {
      stop(paste("dbn.sampling.cpp: transition node", v,
                 "has no t = 0 counterpart"))
    }
  }
  var_index <- stats::setNames(seq_along(base_names) - 1L, base_names)

  n_samples <- as.integer(n_samples)
  max_time <- as.integer(max_time)

  if (dbn_type == "gaussian") {
    plan_0 <- lapply(nodes_0, function(v) plan_gaussian_node(bn_0, v, var_index, g0 = TRUE))
    plan_t <- lapply(nodes_t, function(v) plan_gaussian_node(bn_transition, v, var_index))
    values <- dbn_sample_gaussian_cpp(n_samples, max_time, length(base_names),
                                      plan_0, plan_t)
    columns <- lapply(seq_along(base_names), function(i) values[, i])
  } else if (dbn_type == "discrete") {
    levels_list <- discrete_levels(bn_0, bn_transition, base_names)
    plan_0 <- lapply(nodes_0, function(v) plan_discrete_node(bn_0, v, var_index, levels_list, g0 = TRUE))
    plan_t <- lapply(nodes_t, function(v) plan_discrete_node(bn_transition, v, var_index, levels_list))
    codes <- dbn_sample_discrete_cpp(n_samples, max_time, length(base_names),
                                     vapply(levels_list[base_names], length, integer(1)),
                                     plan_0, plan_t)
    columns <- lapply(seq_along(base_names),
                      function(i) levels_list[[base_names[i]]][codes[, i] + 1L])
  } else if (dbn_type == "mixed") {
    # split the base variables by node type (discrete / gaussian / cgnode)
    discrete_nodes <- base_names[vapply(
      base_names,
      \(n) inherits(fitted_dbn[[paste0(n, "_t")]], "dbn.fit.dnode"),
      logical(1)
    )]
    gaussian_nodes <- base_names[vapply(
      base_names,
      \(n) inherits(fitted_dbn[[paste0(n, "_t")]], "dbn.fit.gnode"),
      logical(1)
    )]
    levels_list <- discrete_levels(bn_0, bn_transition, discrete_nodes)

    parse_0 <- plan_all_mix(bn_0, nodes_0, var_index, levels_list, g0 = TRUE)
    parse_t <- plan_all_mix(bn_transition, nodes_t, var_index, levels_list)

    plan_0 <- parse_0$plan
    ordering_0 <- parse_0$ordering
    plan_t <- parse_t$plan
    ordering_t <- parse_t$ordering

    # the C++ core indexes n_levels by the global column, so it needs one entry
    # per base variable: the level count for discrete columns, 0 elsewhere
    nlevels <- integer(length(base_names))
    for (nm in discrete_nodes)
      nlevels[var_index[[nm]] + 1L] <- length(levels_list[[nm]])

    nvars = list("discrete" = length(discrete_nodes),
                 "gaussian" = length(gaussian_nodes),
                 "mixed" = length(base_names) - length(discrete_nodes)
                             - length(gaussian_nodes))

    values <- dbn_sample_mixed_cpp(n_samples, max_time, nvars,
                                 nlevels, plan_0, plan_t, ordering_0, ordering_t)

    # discrete columns: canonical codes -> labels; continuous columns: as-is
    columns <- lapply(seq_along(base_names), function(i) {
      nm <- base_names[i]
      if (nm %in% discrete_nodes)
        levels_list[[nm]][values[, i] + 1L]
      else
        values[, i]
    })
  }
  names(columns) <- base_names

  timeseries_dict <- c(
    list(
      Time = rep(c(0, seq_len(max_time)), times = n_samples),
      Sample_id = rep(paste("sample", seq_len(n_samples), sep = ""),
                      each = max_time + 1L)
    ),
    columns
  )
  data.frame(timeseries_dict, check.names = FALSE)
}
