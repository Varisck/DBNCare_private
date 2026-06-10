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
# coefficient of parents[i])
plan_gaussian_node <- function(net, variable, var_index) {
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
  par_var <- integer(length(parents))
  par_lag <- integer(length(parents))
  for (i in seq_along(parents)) {
    ref <- plan_parent_ref(parents[i], var_index)
    par_var[i] <- ref$var
    par_lag[i] <- ref$lag
  }
  list(var = var_index[[split_variable_name(variable)$name]],
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
plan_discrete_node <- function(net, variable, var_index, levels_list) {
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

  n_parents <- length(dnn) - 1L
  par_var <- integer(n_parents)
  par_lag <- integer(n_parents)
  par_map <- vector("list", n_parents)
  for (j in seq_len(n_parents)) {
    parent <- dnn[j + 1L]
    ref <- plan_parent_ref(parent, var_index)
    par_var[j] <- ref$var
    par_lag[j] <- ref$lag
    codes <- match(dn[[j + 1L]], levels_list[[split_variable_name(parent)$name]]) - 1L
    codes[is.na(codes)] <- -1L
    par_map[[j]] <- as.integer(codes)
  }
  list(var = var_index[[base]],
       dims = as.integer(dim(prob)),
       own_map = as.integer(own_map),
       freq = as.numeric(prob),
       par_var = par_var,
       par_lag = par_lag,
       par_map = par_map)
}


#' Generate a sampling dataset (Rcpp implementation)
#'
#' Drop-in replacement for [dbn.sampling()] with the sampling loops
#' implemented in C++ via \pkg{Rcpp}. The C++ cores draw from the R random
#' number generator in the same order as the pure-R implementation, so for a
#' fixed seed both functions return the same dataset.
#'
#' @param fitted_dbn an object of class 'dbn.fit'
#' @param n_samples number of samples
#' @param max_time time series length
#' @returns the generated dataframe
#' @useDynLib DynamicBayesianNetwork, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
#' @examples
#' dbn.sampling.cpp(DBN_example, N_samples, Time)
dbn.sampling.cpp <- function(fitted_dbn, n_samples, max_time) {
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
  if (class(fitted_dbn) != "dbn.fit") {
    stop("fitted_DBN must be a dbn.fit object")
  }

  dbn_type <- dbn_type(fitted_dbn)
  if (!dbn_type %in% c("discrete", "gaussian")) {
    stop("Invalid dbn_type")
  }

  bn_0 <- from_fitted_DBN_to_fitted_G_0(fitted_dbn)
  bn_transition <- from_fitted_DBN_to_fitted_G_transition(fitted_dbn)
  # sampling order: node ordering of G_0 at t = 0, then node ordering of the
  # transition network restricted to the _t nodes for t = 1..max_time
  nodes_0 <- bnlearn::node.ordering(bn_0)
  nodes_t <- get_nodes_t(remove_prev_time_from_bn_fit(bn_transition))

  # output columns follow the order in which the variables first enter the
  # time-series dictionary in dbn.sampling(), i.e. the t = 0 node ordering
  base_names <- vapply(nodes_0, function(v) split_variable_name(v)$name,
                       character(1), USE.NAMES = FALSE)
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
    plan_0 <- lapply(nodes_0, function(v) plan_gaussian_node(bn_0, v, var_index))
    plan_t <- lapply(nodes_t, function(v) plan_gaussian_node(bn_transition, v, var_index))
    values <- dbn_sample_gaussian_cpp(n_samples, max_time, length(base_names),
                                      plan_0, plan_t)
    columns <- lapply(seq_along(base_names), function(i) values[, i])
  } else {
    levels_list <- discrete_levels(bn_0, bn_transition, base_names)
    plan_0 <- lapply(nodes_0, function(v) plan_discrete_node(bn_0, v, var_index, levels_list))
    plan_t <- lapply(nodes_t, function(v) plan_discrete_node(bn_transition, v, var_index, levels_list))
    codes <- dbn_sample_discrete_cpp(n_samples, max_time, length(base_names),
                                     vapply(levels_list[base_names], length, integer(1)),
                                     plan_0, plan_t)
    columns <- lapply(seq_along(base_names),
                      function(i) levels_list[[base_names[i]]][codes[, i] + 1L])
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
  data.frame(timeseries_dict)
}
