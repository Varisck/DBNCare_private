# DEFINE FUNCTION FOR dbn.fit OBJECT CREATION -> NODES (with childrens, parents and CPTs)

get_static_nodes <- function(dbn) {
  bnlearn::node.ordering(get.g0.net(dbn))
}

get_dynamic_nodes <- function(dbn) {
  g_transition_graph <- get.transition.net(dbn)
  grep("_t$", bnlearn::node.ordering(g_transition_graph), value = TRUE)
}


#' Function for parameter learning in Dynamic Bayesian Networks
#'
#' Dispatcher for parameter learning. Routes to the discrete, gaussian, or
#' mixed subroutine based on the type of \code{distribution} provided:
#' \itemize{
#'   \item arrays (CPTs)                        -> discrete network
#'   \item named numeric vectors (CPDs)         -> gaussian network
#'   \item lists with \code{$coef}/\code{$sd} (CLG) -> mixed network (not yet supported)
#'   \item only \code{data}                     -> dataset_type(data) decides the subroutine
#' }
#' Exactly one of \code{distribution} or \code{data} must be non-NULL/non-empty.
#'
#' @param DBN object of class 'dbn'
#' @param distribution named list of per-node distributions. Each entry is one of:
#'   a multi-dimensional array (CPT, discrete), a named numeric vector
#'   (regression coefficients + std, gaussian), or a \code{list(coef, sd)} (CLG, mixed).
#'   The type is detected automatically by the function.
#' @param data data.frame object
#' @param replace.unidentifiable If TRUE conditional probabilities for unobserved
#'   parents combinations (unidentifiable parameters) are replaced by uniform
#'   conditional probabilities, if FALSE (default) they are set as NA, in the discrete case.
#'   In the conditional gaussian if replace.unidentifiable is TRUE set variable regressors to marginals
#'   if FALSE they are set to NA.
#' @param compute.priors (default to FALSE) if TRUE marginal distributions of nodes with < markov_order are included,
#'   if FALSE the distributions are conditional to the structure.
#'
#' @description If a data.frame is given, the function inspects the class of each column except
#'   \code{Time} and \code{Sample_id}: if all columns are factor, discrete parameter learning is
#'   used; if all columns are numeric, gaussian parameter learning is used. Mixed case is not
#'   supported yet.
#'
#' @return object of class 'dbn.fit'
#' @export
#'
#' @examples
#' learned_dbn <- dbn.fit(DBN = DBN_example, data = sampling_set)
dbn.fit <- function(DBN, distribution = NULL,
                    data = data.frame(),
                    replace.unidentifiable = FALSE,
                    compute.priors = FALSE) {
  if (!class(DBN) == "dbn") {
    stop("ERROR: DBN argument is not of class 'dbn'")
  }
  if (!class(data) == "data.frame") {
    stop("ERROR: data argument is not of class 'data.frame'")
  }
  if (any(is.na(DBN))) {
    stop("ERROR: missing data detected")
  }

  static_nodes <- get_static_nodes(DBN)
  dynamic_nodes <- get_dynamic_nodes(DBN)

  sources <- c(!is.null(distribution), nrow(data) > 0)
  if (sum(sources) == 0) {
    stop("ERROR: one of distribution or data must be provided
    to learn the parameters of the DBN")
  }
  if (sum(sources) > 1) {
    stop("ERROR: only one of distribution or data must be provided
     to learn the parameters of the DBN")
  }

  if (!is.null(distribution)) {
    # this trows error if distribution is not recognized
    distribution_type <- figure_out_distribution_type(distribution)

    if (distribution_type == "discrete") {
      return(learn_param_d_cpts(DBN,
        CPTs = distribution,
        static_nodes = static_nodes,
        dynamic_nodes = dynamic_nodes
      ))
    }

    if (distribution_type == "gaussian") {
      return(learn_param_g_cpds(DBN,
        CPDs = distribution,
        static_nodes = static_nodes,
        dynamic_nodes = dynamic_nodes
      ))
    }

    if (distribution_type == "mixed") {
      return(learn_param_mixed_dist(DBN,
        distribution = distribution,
        static_nodes = static_nodes,
        dynamic_nodes = dynamic_nodes
      ))
    }
  }

  # only data was provided: dispatch on dataset_type

  # checking data first dispatch later
  vars <- colnames(data)[!colnames(data) %in% c("Time", "Sample_id")]

  if (!setequal(
    names(DBN$nodes),
    setdiff(colnames(data), c("Sample_id", "Time"))
  )) {
    stop("ERROR: nodes in DBN and variables in dataframe do not match")
  }
  if (!(all(dplyr::count(data %>% dplyr::group_by(Sample_id, Time))$n == 1))) {
    stop("ERROR: One or mode combinations of IDs and time slices is repeated")
  }
  if (any(is.na(data %>% tidyr::complete(Sample_id, Time)))) {
    stop(paste(
      "ERROR: One or more sample/individual",
      "present an incomplete temporal sequences"
    ))
  }

  type <- dataset_type(data)
  if (type == "discrete") {
    return(learn_param_d_data(DBN,
      data = data,
      static_nodes = static_nodes,
      dynamic_nodes = dynamic_nodes,
      replace.unidentifiable = replace.unidentifiable,
      compute.priors = compute.priors
    ))
  }
  if (type == "gaussian") {
    return(learn_param_g_data(DBN,
      data = data,
      static_nodes = static_nodes,
      dynamic_nodes = dynamic_nodes,
      compute.priors = compute.priors
    ))
  }

  if (type == "mixed") {
    return(learn_param_data_mixed(DBN,
      data = data,
      static_nodes = static_nodes,
      dynamic_nodes = dynamic_nodes,
      replace.unidentifiable = replace.unidentifiable,
      compute.priors = compute.priors
    ))
  }
}

figure_out_distribution_type <- function(dist) {
  is_clg <- sapply(dist, \(x) is.list(x) && !is.null(x$coef) && !is.null(x$sd))
  is_discrete <- sapply(dist, is.array)
  is_gaussian <- sapply(dist, \(x) is.numeric(x) && is.vector(x))

  if (any(is_clg)) {
    return("mixed")
  }
  if (all(is_gaussian)) {
    return("gaussian")
  }
  if (all(is_discrete)) {
    return("discrete")
  }

  stop(paste(
    "ERROR: could not determine distribution type,",
    "entries must be arrays (CPT), numeric vectors (CPD),",
    "or list(coef, sd) (CLG)"
  ))
}

# ---- discrete subroutines -------------------------------------------------

compute_cpt <- function(data, variable, parents, lvs, replace.unidentifiable) {
  if (length(parents) > 0) {
    pr <-
      data[, c(rev(parents), variable)] %>%
      dplyr::group_by_all() %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      tidyr::complete(!!!rlang::syms(c(rev(parents), variable))) %>%
      replace(is.na(.), 0) %>%
      dplyr::group_by(dplyr::across(rev(parents))) %>%
      dplyr::reframe(!!dplyr::sym(variable), prob = n / sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange_all(.vars = c(rev(parents), variable))
    if (any(is.na(pr$prob))) {
      if (replace.unidentifiable) {
        pr$prob <- replace(pr$prob, is.na(pr$prob), 1 / length(lvs[[variable]]))
      } else {
        warning(paste(
          "WARNING: Probabilities of",
          "the conditioning set equal to 0: Relative frequency is NA"
        ))
      }
    }
    return(pr$prob)
  } else {
    return(unlist(lapply(
      lvs[[variable]],
      \(x) nrow(data[data[, variable] == x, ]) / nrow(data)
    )))
  }
}

fit_nodes_discrete <- function(nodes, data, dbn, lvs,
                               replace.unidentifiable) {
  fitted <- list()
  for (variable in nodes) {
    parents <- get_parent_set(dbn, variable)
    children <- get_children_set(dbn, variable)
    prob <- compute_cpt(data, variable, parents, lvs, replace.unidentifiable)
    fitted[[variable]] <- list(
      node = variable,
      parents = parents,
      children = children,
      prob = array(
        prob,
        dim      = unname(unlist(lapply(lvs[c(variable, parents)], length))),
        dimnames = lvs[c(variable, parents)]
      )
    )
    class(fitted[[variable]]) <- "dbn.fit.dnode"
  }
  fitted
}

learn_param_d_data <- function(dbn, data, static_nodes, dynamic_nodes,
                               replace.unidentifiable, compute.priors) {
  # df_0: time-0 slice with _0 suffix
  # get data.set for time 0 slice
  # substitue names of variables with _0 at the end
  # df_0 <- data[data$Time == 0, ]
  # names(df_0) <- lapply(names(df_0), concat_name_post, postfix = "_0")

  markov_order <- dbn$markov_order
  df_0 <- build_g0_df(data, markov_order = markov_order)

  # if markov order > 1 compute values of p(A_t=lvls), ...
  if (compute.priors) {
    if (markov_order > 1) {
      priors <- list()
      for (t in seq_len(markov_order - 1)) {
        df_t <- data[data$Time == t, ]
        for (var in vars) {
          key <- paste0(var, "_", t)
          counts <- table(df_0[[key]])
          priors[[key]] <- counts / sum(counts)
        }
      }
      dbn_fitted[["priors"]] <- priors
    }
  }

  # build transition dataframe A_t, ..., A_t-i
  df_transition <- build_shifted_df(data,
    markov_order = markov_order, separator = "-"
  )

  # build levels lookup for each node at _0, _t and _t-k
  lvs <- list()
  for (var in names(df_transition)) {
    lvs[[var]] <- sort(as.array(levels(factor(df_transition[[var]]))))
  }
  for (var in names(df_0)) {
    lvs[[var]] <- sort(as.array(levels(factor(df_0[[var]]))))
  }

  dbn_fitted <- c(
    fit_nodes_discrete(
      static_nodes, df_0, dbn,
      lvs, replace.unidentifiable
    ),
    fit_nodes_discrete(
      dynamic_nodes, df_transition, dbn,
      lvs, replace.unidentifiable
    )
  )

  class(dbn_fitted) <- "dbn.fit"
  dbn_fitted
}

# checks the necessary CPT errors and returns the node_info object
get_node_info_discrete <- function(CPT, variable, parents, children, defined_levels) {
  # check dimentions against parent set
  if (!setequal(setdiff(names(dimnames(CPT)), variable), parents)) {
    stop("ERROR: CPTs do not match parents set")
  }

  # checking parent levels metch stored levels of variable
  if (length(dimnames(CPT)) > 1) {
    dim_names <- dimnames(CPT)[2:length(dimnames(CPT))]
    for (parent in names(dim_names)) {
      parent_name <- get_variable_name(parent)
      if (!(setequal(dim_names[[parent]], defined_levels[[parent_name]]))) {
        stop(paste(
          "ERROR: Inconsistency in node's levels for variable",
          variable, "expected", toString(defined_levels[[parent_name]]),
          "for parent", parent_name, "got", toString(dim_names[[parent]])
        ))
      }
    }
  }

  list(
    node = variable,
    parents = parents,
    children = children,
    prob = aperm(CPT, c(variable, parents))
  )
}

learn_param_d_cpts <- function(dbn, CPTs, static_nodes, dynamic_nodes) {
  nodes <- c(static_nodes, dynamic_nodes)

  if (!setequal(names(CPTs), nodes)) {
    stop("ERROR: nodes in DBN and variables in CPTs do not match")
  }
  if (!(any(lapply(CPTs, class) %in% c("matrix", "array")))) {
    stop("ERROR: CPT must be of class 'matrix' or 'array'")
  }

  defined_levels <- list()
  nodes_info <- list()

  for (variable in nodes) {
    CPT <- CPTs[[variable]]
    # checking numeric and all probs sum to 1
    if (!(all(lapply(CPT, class) == "numeric"))) {
      stop("ERROR: Probabilities must be numeric")
    }
    if (length(dim(CPT)) > 1) {
      l <- length(dim(CPT))
      idx_target <- which(names(dimnames(CPT)) == variable)
      if (!(all(apply(CPT, setdiff(1:l, idx_target), sum) == as.character(1)))) {
        stop("ERROR: Probabilities for each conditioning set must sum to 1")
      }
    } else {
      if (sum(CPT) != as.character(1)) {
        stop("ERROR: Probabilities for each conditioning set must sum to 1")
      }
    }

    parents <- get_parent_set(dbn, variable)
    children <- get_children_set(dbn, variable)

    var_name <- get_variable_name(variable)
    defined_levels[[var_name]] <- dimnames(CPT)[[variable]]

    nodes_info[[variable]] <- get_node_info_discrete(
      CPT, variable, parents,
      children, defined_levels
    )
    class(nodes_info[[variable]]) <- "dbn.fit.dnode"
  }

  class(nodes_info) <- "dbn.fit"
  nodes_info
}


# ---- gaussian subroutines -------------------------------------------------

get_node_info_gaussian <- function(cpd, variable, parents, children) {
  if (is.list(cpd)) {
    stop(paste(
      "ERROR: CPD must be a numeric vector found",
      toString(class(cpd))
    ))
  }

  # check values in CPD are numeric
  if (!(all(lapply(cpd, class) == "numeric"))) {
    stop("ERROR: Probabilities must be numeric")
  }

  # check CPD is parents + intercept + std
  if (length(cpd) != (length(parents) + 2)) {
    stop(paste(
      "ERROR: Variable", variable, ".Expected CPD of length",
      length(parents) + 2, "got ", length(cpd)
    ))
  }

  # check that CPD names are ordered same as parents
  if (!identical(names(cpd)[seq_along(parents) + 1], parents)) {
    stop(paste(
      "ERROR: Variable", variable,
      "regressors have to be ordered according to parents!"
    ))
  }

  # check first and last elements are intercept and std
  if (names(cpd)[1] != intercept_name |
    names(cpd)[length(cpd)] != std_name) {
    stop(paste(
      "ERROR: Variable", variable,
      "first and last values of parameters must be",
      intercept_name, "and", std_name
    ))
  }

  return(list(
    node = variable,
    parents = parents,
    children = children,
    regs = cpd[1:length(cpd) - 1],
    std = cpd[length(cpd)]
  ))
}

learn_param_g_cpds <- function(dbn, CPDs, static_nodes, dynamic_nodes) {
  nodes <- c(static_nodes, dynamic_nodes)
  intercept_std <- c(intercept_name, std_name)

  nodes_info <- list()

  # check that CPDs contains distributions for each node in the net
  if (!setequal(names(CPDs), nodes)) {
    stop("ERROR: nodes in DBN and variables in CPDs do not match")
  }

  for (variable in names(CPDs)) {
    parents <- get_parent_set(dbn, variable)
    children <- get_children_set(dbn, variable)

    cpd <- CPDs[[variable]]

    nodes_info[[variable]] <- get_node_info_gaussian(cpd, variable, parents, children)
    class(nodes_info[[variable]]) <- "dbn.fit.gnode"
  }

  class(nodes_info) <- "dbn.fit"
  nodes_info
}

replace_minus_unerscore <- function(name) {
  gsub("-([0-9]+)$", "_\\1", name)
}

replace_underscore_minus <- function(name) {
  gsub("t_([0-9]+)$", "t-\\1", name)
}

# given the data variable and parent set find the parameters of the linear reg.
get_variable_model <- function(data, variable, parents) {
  bt <- function(x) paste0("`", x, "`")
  if (length(parents) == 0) {
    formula <- as.formula(paste(bt(variable), "~ 1"))
  } else {
    formula <- reformulate(sapply(parents, \(p) bt(p)),
      response = bt(variable)
    )
  }

  # run linear regression model
  model <- lm(formula = formula, data = data)

  coeff <- model$coefficients
  # removing the `` used to escape the - character
  names(coeff) <- sapply(names(coeff), \(name) gsub("`", "", name))
  std <- sigma(model)
  return(
    list(
      coeff = coeff,
      std = sigma(model)
    )
  )
}

# find priors of 1:markov_order for all columns
get_var_priors <- function(markov_order, data) {
  priors <- list()

  columns <- setdiff(names(data), c("Sample_id", "Time"))

  for (mo in seq(markov_order - 1)) {
    d <- data[data$Time == mo, ]
    for (col in columns) {
      priors[[paste0(col, "_", mo)]] <-
        get_variable_model(d, col, c())
    }
  }
  priors
}


fit_nodes_gaussian <- function(nodes, data, dbn) {
  nodes_info <- list()
  for (variable in nodes) {
    parents <- get_parent_set(dbn, variable)
    children <- get_children_set(dbn, variable)
    res <- get_variable_model(data, variable, parents)

    nodes_info[[variable]] <- list(
      node = variable,
      parents = parents,
      children = children,
      regs = res$coeff,
      std = res$std
    )
    class(nodes_info[[variable]]) <- "dbn.fit.gnode"
  }
  nodes_info
}


learn_param_g_data <- function(dbn, data, static_nodes, dynamic_nodes, compute.priors) {
  nodes_info <- list()

  nodes <- c(static_nodes, dynamic_nodes)
  markov_order <- dbn$markov_order

  # build the G_0 frame. When markov_order > 1 this spans the initial slices
  # (var_0, var_1, ..., var_{markov_order-1}) so the extended G_0 nodes can be
  # fitted, mirroring the discrete path (learn_param_d_data).
  df_0 <- build_g0_df(data, markov_order = markov_order)

  # if markov order > 1 compute values of p(A_t=1), ...
  if (markov_order > 1) {
    df_priors <- data[data$Time %in% seq(markov_order), ]
    priors <- get_var_priors(markov_order, df_priors)
    nodes_info[["priors"]] <- priors
  }

  # build transition dataframe A_t, ..., A_t-i
  df_transition <- build_shifted_df(data, markov_order = markov_order, separator = "-")

  nodes_info <- c(
    fit_nodes_gaussian(static_nodes, df_0, dbn),
    fit_nodes_gaussian(dynamic_nodes, df_transition, dbn)
  )
  class(nodes_info) <- "dbn.fit"
  nodes_info
}


# ---- Mixed subroutines -------------------------------------------------

get_node_info_mixed <- function(dist, variable, parents, children, defined_levels) {
  if (!is.array(dist$coef)) {
    stop(paste(
      "ERROR: distribution of variable", variable,
      "attribute coef must be of class array"
    ))
  }
  if (!is.numeric(dist$sd)) {
    stop(paste(
      "ERROR: distribution of variable", variable,
      "attribute sd must be of of class numeric"
    ))
  }

  coef <- dist$coef

  if (!identical(unname(sapply(dimnames(coef), length)), dim(coef))) {
    stop(paste(
      "ERROR: distribution of variable", variable,
      "discrete parents levels size is mismatched"
    ))
  }

  add_parents <- c()
  if (length(dimnames(coef)) > 1) {
    add_parents <- setdiff(dimnames(coef)[[1]], intercept_name)
    # checking levels for discrete parents match the defined ones
    dim_names <- dimnames(coef)[2:length(dimnames(coef))]
    for (parent in names(dim_names)) {
      parent_name <- get_variable_name(parent)
      if (!(setequal(dim_names[[parent]], defined_levels[[parent_name]]))) {
        stop(paste(
          "ERROR: Inconsistency in node's levels for variable",
          variable, "expected", toString(defined_levels[[parent_name]]),
          "for parent", parent_name, "got", toString(dim_names[[parent]])
        ))
      }
    }
  }

  if (!setequal(setdiff(union(names(dimnames(coef)), add_parents), variable), parents)) {
    stop(paste("ERROR: distribution of variable", variable, "do not match parent set"))
  }

  if (length(dim(coef)) > 1) {
    num_regs <- prod(dim(coef)[2:length(dim(coef))])
    if (length(dist$sd) != num_regs) {
      stop(paste(
        "ERROR: distribution of variable", variable, "number of regressions and",
        "number of sd do not match:", num_regs, "against", length(dist$sd)
      ))
    }
  }

  d_levels <- dimnames(coef)[2:length(dimnames(coef))]
  comb <- (1:nrow(expand.grid(d_levels)))
  # comb = list(as.character(1:nrow(expand.grid(d_levels, stringsAsFactors = FALSE))))

  prob <- array(c(coef), dim = c(dim(coef)[1], num_regs), dimnames = list(
    variable = dimnames(coef)[[1]],
    comb
  ))
  # this sets dimnames(coef)[[1]] above with the name of the variable
  names(dimnames(prob))[1] <- variable

  d_parents <- match(names(d_levels), parents)
  g_parents <- match(setdiff(dimnames(coef)[[1]], intercept_name), parents)

  return(list(
    node = variable,
    parents = parents,
    children = children,
    coefficients = prob,
    sd = dist$sd,
    dlevels = d_levels,
    dparents = d_parents,
    gparents = g_parents
  ))
}

learn_param_mixed_dist <- function(dbn, distribution, static_nodes, dynamic_nodes) {
  is_clg <- \(x) is.list(x) && !is.null(x$coef) && !is.null(x$sd)
  is_discrete <- is.array
  is_gaussian <- \(x) is.numeric(x) && is.vector(x)

  nodes <- c(static_nodes, dynamic_nodes)
  intercept_std <- c(intercept_name, std_name)

  defined_levels <- list()
  nodes_info <- list()

  # check that distribution contains distributions for each node in the net
  if (!setequal(names(distribution), nodes)) {
    stop("ERROR: nodes in DBN and variables in distribution do not match")
  }

  # parents may be lagged "D_t-1", map each base variable name to its
  # distribution to look up parent types regardless of slice.
  dist_by_var <- list()
  for (n in names(distribution)) {
    dist_by_var[[get_variable_name(n)]] <- distribution[[n]]
  }

  for (variable in names(distribution)) {
    parents <- get_parent_set(dbn, variable)
    children <- get_children_set(dbn, variable)

    dist <- distribution[[variable]]

    if (is_discrete(dist)) {
      # discrete can only have discrete parents
      parent_dists <- lapply(parents, \(p) dist_by_var[[get_variable_name(p)]])
      if (!all(sapply(parent_dists, is_discrete))) {
        stop(paste(
          "ERROR: found discrete variable", variable,
          "with non discrete parents"
        ))
      }
      # update levels map and get node_info
      var_name <- get_variable_name(variable)
      defined_levels[[var_name]] <- dimnames(dist)[[variable]]
      nodes_info[[variable]] <- get_node_info_discrete(
        dist, variable, parents,
        children, defined_levels
      )
      class(nodes_info[[variable]]) <- "dbn.fit.dnode"
    } else if (is_gaussian(dist)) {
      nodes_info[[variable]] <- get_node_info_gaussian(
        dist, variable,
        parents, children
      )
      class(nodes_info[[variable]]) <- "dbn.fit.gnode"
    } else if (is_clg(dist)) {
      nodes_info[[variable]] <- get_node_info_mixed(
        dist, variable, parents,
        children, defined_levels
      )
      class(nodes_info[[variable]]) <- "dbn.fit.cgnode"
    } else {
      stop(paste(
        "ERROR: found unrecognized node type in distribution for variable", variable, "please use:",
        "\n1) a multi-dimensional array (CPT, discrete)",
        "\n2) a named numeric vector (regression coefficients + std, gaussian)",
        "\n3) a list(coef, sd) (CLG, mixed)"
      ))
    }
  }

  class(nodes_info) <- "dbn.fit"
  nodes_info
}

discrete_columns <- function(df) {
  names(df)[sapply(df, \(x) is.factor(x) || is.character(x))]
}

get_node_distribution_mixed <- function(data, variable, parents, children, lvs,
                                        replace.unidentifiable) {
  # get names of discrete parents
  discrete_parents <- intersect(parents, names(lvs))

  d_parents <- match(discrete_parents, parents)
  g_parents <- match(setdiff(parents, discrete_parents), parents)
  d_levels <- lvs[discrete_parents]

  combos <- expand.grid(d_levels, stringsAsFactors = FALSE)

  coef_list <- vector("list", nrow(combos))
  sd_vec <- numeric(nrow(combos))

  # this is so that each res$coef from regression is always in the same order
  regression_names <- c(intercept_name, parents[g_parents])

  for (i in seq_len(nrow(combos))) {
    mask <- Reduce("&", lapply(
      names(combos),
      \(col) data[[col]] == combos[i, col]
    ))
    # no combination of discrete parents exist in df
    if (sum(mask) == 0) {
      if (replace.unidentifiable) {
        warning(paste(
          "WARNING: no data for discrete parent combination", i,
          "for variable", variable, "- using marginal model"
        ))
        coef_list[[i]] <- setNames(
          c(mean(data[[variable]], na.rm = TRUE), rep(0, length(parents[g_parents]))),
          regression_names
        )
        sd_vec[[i]] <- sd(data[[variable]], na.rm = TRUE)
      } else {
        warning(paste(
          "WARNING: variable", variable, "comination of discrete parents",
          "not in dataframe, replace.unidentifiable set to FALSE - setting to NA"
        ))
        coef_list[[i]] <- setNames(rep(NA_real_, length(regression_names)), regression_names)
        sd_vec[[i]] <- NA_real_
      }
    } else {
      res <- get_variable_model(data[mask, ], variable, parents[g_parents])
      coef_list[[i]] <- res$coeff[regression_names]
      sd_vec[[i]] <- res$std
    }
  }

  coef_array <- array(
    unlist(coef_list),
    dim = c(length(regression_names), prod(sapply(d_levels, length))),
    dimnames = c(list(regression_names), list(as.character(1:nrow(combos))))
  )
  names(dimnames(coef_array))[1] <- variable

  return(list(
    node = variable,
    parents = parents,
    children = children,
    coefficients = coef_array,
    sd = sd_vec,
    dlevels = d_levels,
    dparents = d_parents,
    gparents = g_parents
  ))
}


fit_nodes_mixed <- function(nodes, data, dbn, lvs, replace.unidentifiable) {
  nodes_info <- list()
  for (variable in nodes) {
    parents <- get_parent_set(dbn, variable)
    children <- get_children_set(dbn, variable)

    # discrete case
    if (is.factor(data[[variable]]) || is.character(data[[variable]])) {
      prob <- compute_cpt(data, variable, parents, lvs, replace.unidentifiable)
      nodes_info[[variable]] <- list(
        node = variable,
        parents = parents,
        children = children,
        prob = array(
          prob,
          dim      = unname(unlist(lapply(lvs[c(variable, parents)], length))),
          dimnames = lvs[c(variable, parents)]
        )
      )
      class(nodes_info[[variable]]) <- "dbn.fit.dnode"
    } else if (is.numeric(data[[variable]])) {
      # here need to check if variable is gaussian or conditional gaussian
      # if all parents are gaussian -> variable is gaussian
      if (all(sapply(parents, \(p) is.numeric(data[[p]])))) {
        # gaussian
        res <- get_variable_model(data, variable, parents)
        nodes_info[[variable]] <- list(
          node = variable,
          parents = parents,
          children = children,
          regs = res$coeff,
          std = res$std
        )
        class(nodes_info[[variable]]) <- "dbn.fit.gnode"
        # still check if some parents are discrete
      } else if (any(sapply(parents, \(p) is.factor(data[[p]]) ||
        is.character(data[[p]])))) {
        # conditional gaussian
        nodes_info[[variable]] <- get_node_distribution_mixed(
          data, variable,
          parents, children,
          lvs,
          replace.unidentifiable
        )
        class(nodes_info[[variable]]) <- "dbn.fit.cgnode"
      } else { # some error in parent type
        stop(paste(
          "ERROR: dbn.fit unrecognize data type in parents of variable", variable,
          "variable is not gaussian nor conditional gaussian"
        ))
      }
    }
  }
  nodes_info
}


learn_param_data_mixed <- function(dbn, data, static_nodes, dynamic_nodes,
                                   replace.unidentifiable, compute.priors) {
  #
  df_0 <- data[data$Time == 0, ]
  names(df_0) <- lapply(names(df_0), concat_name_post, postfix = "_0")

  # if markov order > 1 compute values of p(A_t=lvls), ...
  markov_order <- dbn$markov_order
  # if (markov_order > 1) {
  #   priors <- list()
  #   for (t in seq_len(markov_order - 1)) {
  #     df_t <- data[data$Time == t, ]
  #     for (var in vars) {
  #       key <- paste0(var, "_", t)
  #       counts   <- table(df_t[[var]])
  #       priors[[key]] <- counts / sum(counts)
  #     }
  #   }
  #   dbn_fitted[["priors"]] <- priors
  # }

  # build tr ansition dataframe A_t, ..., A_t-i
  df_transition <- build_shifted_df(data, markov_order = markov_order, separator = "-")

  # build levels lookup for each node at _0, _t and _t-k

  lvs <- list()
  for (var in discrete_columns(df_transition)) {
    lvs[[var]] <- sort(as.array(levels(factor(df_transition[[var]]))))
  }
  for (var in discrete_columns(df_0)) {
    lvs[[var]] <- sort(as.array(levels(factor(df_0[[var]]))))
  }

  nodes_info <- c(
    fit_nodes_mixed(static_nodes, df_0, dbn, lvs, replace.unidentifiable),
    fit_nodes_mixed(
      dynamic_nodes,
      df_transition, dbn, lvs, replace.unidentifiable
    )
  )
  class(nodes_info) <- "dbn.fit"
  nodes_info
}
