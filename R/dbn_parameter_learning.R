# DEFINE FUNCTION FOR dbn.fit OBJECT CREATION -> NODES (with childrens, parents and CPTs)

get_static_nodes = function(dbn) {
  bnlearn::node.ordering(get.g0.net(dbn))
}

get_dynamic_nodes = function(dbn) {
  g_transition_graph = get.transition.net(dbn)
  quanteda::char_select(bnlearn::node.ordering(g_transition_graph),
                        "*t",
                        valuetype = "glob")
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
#'   conditional probabilities, if FALSE (default) they are set as NA. Discrete only.
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
                    replace.unidentifiable = FALSE) {
  if (!class(DBN) == 'dbn')
    stop("ERROR: DBN argument is not of class 'dbn'")
  if (!class(data) == 'data.frame')
    stop("ERROR: data argument is not of class 'data.frame'")
  if (any(is.na(DBN)))
    stop("ERROR: missing data detected")

  static_nodes = get_static_nodes(DBN)
  dynamic_nodes = get_dynamic_nodes(DBN)

  sources = c(!is.null(distribution), nrow(data) > 0)
  if (sum(sources) == 0)
    stop("ERROR: one of distribution or data must be provided to learn the parameters of the DBN")
  if (sum(sources) > 1)
    stop("ERROR: only one of distribution or data must be provided to learn the parameters of the DBN")

  if (!is.null(distribution)) {
    # this trows error if distribution is not recognized
    distribution_type = figure_out_distribution_type(distribution)

    if (distribution_type == 'discrete')
      return(learn_param_d_cpts(DBN, CPTs = distribution,
                                static_nodes = static_nodes,
                                dynamic_nodes = dynamic_nodes))

    if (distribution_type == 'gaussian')
      return(learn_param_g_cpds(DBN, CPDs = distribution,
                                static_nodes = static_nodes,
                                dynamic_nodes = dynamic_nodes))

    if (distribution_type == 'mixed')
      return(learn_param_mixed(DBN, distribution = distribution,
                               static_nodes = static_nodes,
                               dynamic_nodes = dynamic_nodes))
  }

  # only data was provided: dispatch on dataset_type

  # checking data first dispatch later
  vars = colnames(data)[!colnames(data) %in% c('Time', 'Sample_id')]

  if (!setequal(names(DBN$nodes), setdiff(colnames(data), c('Sample_id', 'Time'))))
    stop("ERROR: nodes in DBN and variables in dataframe do not match")
  if (!(all(dplyr::count(data %>% dplyr::group_by(Sample_id, Time))$n == 1)))
    stop("ERROR: One or mode combinations of IDs and time slices is repeated")
  if (any(is.na(data %>% tidyr::complete(Sample_id, Time))))
    stop("ERROR: One or more sample/individual present an incomplete temporal sequences")

  type = dataset_type(data)
  if (type == "discrete")
    return(learn_param_d_data(DBN, data = data,
                              static_nodes = static_nodes,
                              dynamic_nodes = dynamic_nodes,
                              replace.unidentifiable = replace.unidentifiable))
  if (type == "gaussian")
    return(learn_param_g_data(DBN, data = data,
                              static_nodes = static_nodes,
                              dynamic_nodes = dynamic_nodes))

  stop("ERROR: mixed datasets are not supported yet")
}

figure_out_distribution_type = function(dist) {
  is_clg      <- sapply(dist, \(x) is.list(x) && !is.null(x$coef) && !is.null(x$sd))
  is_discrete <- sapply(dist, is.array)
  is_gaussian <- sapply(dist, \(x) is.numeric(x) && is.vector(x))

  if (any(is_clg))       return("mixed")
  if (all(is_gaussian))  return("gaussian")
  if (all(is_discrete))  return("discrete")

  stop("ERROR: could not determine distribution type — entries must be arrays (CPT), numeric vectors (CPD), or list(coef, sd) (CLG)")
}

# ---- discrete subroutines -------------------------------------------------

compute_cpt = function(data, variable, parents, lvs, replace.unidentifiable) {
  if(length(parents) > 0){
    pr <-
    data[, c(rev(parents), variable)] %>% 
    dplyr::group_by_all() %>% 
    dplyr::count() %>% 
    dplyr::ungroup() %>% 
    tidyr::complete(!!! rlang::syms(c(rev(parents), variable))) %>% 
    replace(is.na(.), 0) %>% 
    dplyr::group_by(dplyr::across(rev(parents))) %>% 
    dplyr::reframe(!!dplyr::sym(variable), prob = n / sum(n)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange_all(.vars = c(rev(parents), variable))
    if (any(is.na(pr$prob))) {
      if (replace.unidentifiable) {
        pr$prob <- replace(pr$prob, is.na(pr$prob), 1 / length(lvs[[variable]]))
      } else {
        warning("WARNING: Probabilities of the conditioning set equal to 0: Relative frequency is NA")
      }
    }
    return(pr$prob)
  } else {
    return(unlist(lapply(lvs[[variable]], \(x)
        nrow(data[data[, variable] == x, ]) / nrow(data))))
  }  
}

fit_nodes_discrete <- function(nodes, data, dbn, lvs, replace.unidentifiable) {
  fitted <- list()
  for (variable in nodes) {
    parents  <- get_parent_set(dbn, variable)
    children <- get_children_set(dbn, variable)
    prob <- compute_cpt(data, variable, parents, lvs, replace.unidentifiable)
    fitted[[variable]] <- list(
      node     = variable,
      parents  = parents,
      children = children,
      prob     = array(
        prob,
        dim      = unname(unlist(lapply(lvs[c(variable, parents)], length))),
        dimnames = lvs[c(variable, parents)]
      )
    )
    class(fitted[[variable]]) <- "dbn.fit.dnode"
  }
  fitted
}

learn_param_d_data = function(dbn, data,
                              static_nodes,
                              dynamic_nodes,
                              replace.unidentifiable) {
  # df_0: time-0 slice with _0 suffix
  # get data.set for time 0 slice
  # substitue names of variables with _0 at the end
  df_0 = data[data$Time == 0,]
  names(df_0) = lapply(names(df_0), concat_name_post, postfix = "_0")

  # if markov order > 1 compute values of p(A_t=lvls), ...
  markov_order = dbn$markov_order
  if (markov_order > 1) {
    priors <- list()
    for (t in seq_len(markov_order - 1)) {
      df_t <- data[data$Time == t, ]
      for (var in vars) {
        key <- paste0(var, "_", t)
        counts   <- table(df_t[[var]])
        priors[[key]] <- counts / sum(counts)
      }
    }
    dbn_fitted[["priors"]] <- priors
  }

  # build transition dataframe A_t, ..., A_t-i
  df_transition = build_shifted_df(data, markov_order = markov_order, separator = "-")

  # build levels lookup for each node at _0, _t and _t-k
  lvs = list()
  for (var in names(df_transition)) {
    lvs[[var]] <- sort(as.array(levels(factor(df_transition[[var]]))))
  }
  for (var in names(df_0)) {
    lvs[[var]] <- sort(as.array(levels(factor(df_0[[var]]))))
  }

  dbn_fitted <- c(
    fit_nodes_discrete(static_nodes,  df_0,          dbn, lvs, replace.unidentifiable),
    fit_nodes_discrete(dynamic_nodes, df_transition, dbn, lvs, replace.unidentifiable)
  )

  class(dbn_fitted) <- "dbn.fit"
  dbn_fitted
}

# checks the necessary CPT errors and returns the node_info object
get_node_info_discrete = function(CPT, variable, parents, children, defined_levels) {
  # check dimentions against parent set
  if (!setequal(setdiff(names(dimnames(CPT)), variable), parents))
    stop("ERROR: CPTs do not match parents set")

  # checking parent levels metch stored levels of variable
  if(length(dimnames(CPT)) > 1){ 
    dim_names = dimnames(CPT)[2:length(dimnames(CPT))]
    for (parent in names(dim_names)) {
      parent_name = get_variable_name(parent)
      if (!(setequal(dim_names[[parent]], defined_levels[[parent_name]])))
        stop(paste("ERROR: Inconsistency in node's levels for variable",
                variable, "expected", toString(defined_levels[[parent_name]]),
                "for parent", parent_name, "got", toString(dim_names[[parent]])))
    }
  }

  return(list(
      node = variable,
      parents = parents,
      children = children,
      prob = aperm(CPT, c(variable, parents))
    ))
}

learn_param_d_cpts = function(dbn, CPTs,
                              static_nodes,
                              dynamic_nodes) {
  nodes = c(static_nodes, dynamic_nodes)

  if (!setequal(names(CPTs), nodes))
    stop("ERROR: nodes in DBN and variables in CPTs do not match")
  if (!(any(lapply(CPTs, class) %in% c('matrix', 'array'))))
    stop("ERROR: CPT must be of class 'matrix' or 'array'")

  defined_levels <- list()
  nodes_info <- list()

  for (variable in nodes) {
    CPT <- CPTs[[variable]]
    # checking numeric and all probs sum to 1
    if (!(all(lapply(CPT, class) == 'numeric')))
      stop("ERROR: Probabilities must be numeric")
    if (length(dim(CPT)) > 1) {
      l <- length(dim(CPT))
      idx_target <- which(names(dimnames(CPT)) == variable)
      if (!(all(apply(CPT, setdiff(1:l, idx_target), sum) == as.character(1))))
        stop("ERROR: Probabilities for each conditioning set must sum to 1")
    } else {
      if (sum(CPT) != as.character(1))
        stop("ERROR: Probabilities for each conditioning set must sum to 1")
    }

    parents = get_parent_set(dbn, variable)
    children = get_children_set(dbn, variable)

    var_name = get_variable_name(variable)
    defined_levels[[var_name]] = dimnames(CPT)[[variable]]

    nodes_info[[variable]] = get_node_info_discrete(CPT, variable, parents,
                                                    children, defined_levels)
    class(nodes_info[[variable]]) = "dbn.fit.dnode"
  }

  class(nodes_info) <- "dbn.fit"
  nodes_info
}


# ---- gaussian subroutines -------------------------------------------------

get_node_info_gaussian = function(cpd, variable, parents, children) {
    if (is.list(cpd)) 
      stop(paste("ERROR: CPD must be a numeric vector found", toString(class(cpd))))

    # check values in CPD are numeric
    if(!(all(lapply(cpd, class) == 'numeric')))
      stop("ERROR: Probabilities must be numeric")

    # check CPD is parents + intercept + std
    if(length(cpd) != (length(parents) + 2))
      stop(paste("ERROR: Variable", variable, ".Expected CPD of length",
                 length(parents) + 2, "got ", length(cpd)))

    # check that CPD names are ordered same as parents
    if (!identical(names(cpd)[seq_along(parents) + 1], parents)) {
      stop(paste("ERROR: Variable", variable,
                 "regressors have to be ordered according to parents!"))
    }

    # check first and last elements are intercept and std
    if(names(cpd)[1] != intercept_name |
       names(cpd)[length(cpd)] != std_name)
      stop(paste("ERROR: Variable", variable,
                 "first and last values of parameters must be",
                 intercept_name, "and", std_name))

    return(list(
      node = variable,
      parents = parents,
      children = children,
      regs = cpd[1:length(cpd) - 1],
      std = cpd[length(cpd)]
    ))
}

learn_param_g_cpds = function(dbn, CPDs = list(),
                              static_nodes = static_nodes,
                              dynamic_nodes = dynamic_nodes) {

  nodes = c(static_nodes, dynamic_nodes)
  intercept_std = c(intercept_name, std_name)

  nodes_info = list()

  # check that CPDs contains distributions for each node in the net
  if (!setequal(names(CPDs), nodes))
    stop("ERROR: nodes in DBN and variables in CPDs do not match")

  for(variable in names(CPDs)) {

    parents = get_parent_set(dbn, variable)
    children = get_children_set(dbn, variable)

    cpd = CPDs[[variable]]

    nodes_info[[variable]] = get_node_info_gaussian(cpd, variable, parents, children)
    class(nodes_info[[variable]]) = "dbn.fit.gnode"
  }

  class(nodes_info) = "dbn.fit"
  nodes_info
}

replace_minus_unerscore = function(name) {
  gsub("-([0-9]+)$", "_\\1", name)
}

replace_underscore_minus = function(name) {
  gsub("t_([0-9]+)$", "t-\\1", name)
}

# given the data variable and parent set find the parameters of the linear reg.
get_variable_model = function(data, variable, parents) {
  if (length(parents) == 0) {
    formula <- as.formula(paste(variable, "~ 1"))
  } else {
    formula <- reformulate(sapply(parents, replace_minus_unerscore),
                           variable)
  }

  # run linear regression model
  model = lm(formula = formula, data = data)

  coeff = model$coefficients
  names(coeff) = sapply(names(coeff), replace_underscore_minus)
  std = sigma(model)
  return(
    list(
      coeff = coeff,
      std = sigma(model)
    )
  )
}

# find priors of 1:markov_order for all columns
get_var_priors = function(markov_order, data) {
  priors = list()

  columns = setdiff(names(data), c("Sample_id", "Time"))

  for(mo in seq(markov_order - 1)) {
    d = data[data$Time == mo, ]
    for(col in columns) {
      priors[[paste0(col, "_", mo)]] =
        get_variable_model(d, col, c())
    }
  }
  priors
}


learn_param_g_data = function(dbn, data = data.frame(),
                              static_nodes = static_nodes,
                              dynamic_nodes = dynamic_nodes) {

  nodes_info = list()

  nodes = c(static_nodes, dynamic_nodes)
  # get data.set for time 0 slice
  # substitue names of variables with _0 at the end
  df_0 = data[data$Time == 0,]
  names(df_0) = lapply(names(df_0), concat_name_post, postfix = "_0")

  # if markov order > 1 compute values of p(A_t=1), ...
  markov_order = dbn$markov_order
  if(markov_order > 1) {
    df_priors = data[data$Time %in% seq(markov_order), ]
    priors = get_var_priors(markov_order, df_priors)
    nodes_info[["priors"]] = priors
  }

  # build transition dataframe A_t, ..., A_t-i
  df_transition = build_shifted_df(data, markov_order = markov_order, separator = "_")

  for(variable in static_nodes) {
    parents = get_parent_set(dbn, variable)
    children = get_children_set(dbn, variable)
    res = get_variable_model(df_0, variable, parents)

    nodes_info[[variable]] = list(
      node = variable,
      parents = parents,
      children = children,
      regs = res$coeff,
      std = res$std
    )
    class(nodes_info[[variable]]) = "dbn.fit.gnode"
  }

  for(variable in dynamic_nodes) {
    parents = get_parent_set(dbn, variable)
    children = get_children_set(dbn, variable)
    res = get_variable_model(df_transition, variable, parents)

    nodes_info[[variable]] = list(
      node = variable,
      parents = parents,
      children = children,
      regs = res$coeff,
      std = res$std
    )
    class(nodes_info[[variable]]) = "dbn.fit.gnode"
  }

  class(nodes_info) = "dbn.fit"
  nodes_info
}


# ---- Mixed subroutines -------------------------------------------------

get_node_info_mixed = function(dist, variable, parents, 
                               children, defined_levels) {
  if(!is.array(dist$coef))
    stop(paste("ERROR: distribution of variable", variable,
               "attribute coef must be of class array"))
  if(!is.numeric(dist$sd))
    stop(paste("ERROR: distribution of variable", variable,
               "attribute sd must be of of class numeric"))

  coef = dist$coef

  if (!identical(unname(sapply(dimnames(coef), length)), dim(coef)))
    stop(paste("ERROR: distribution of variable", variable,
                "discrete parents levels size is mismatched"))

  add_parents = c()
  if(length(dimnames(coef)) > 1) {
    add_parents = setdiff(dimnames(coef)[[1]], intercept_name)
    # checking levels for discrete parents match the defined ones
    dim_names = dimnames(coef)[2:length(dimnames(coef))]
    for (parent in names(dim_names)) {
      parent_name = get_variable_name(parent)
      if (!(setequal(dim_names[[parent]], defined_levels[[parent_name]])))
        stop(paste("ERROR: Inconsistency in node's levels for variable",
                variable, "expected", toString(defined_levels[[parent_name]]),
                "for parent", parent_name, "got", toString(dim_names[[parent]])))
    }
  }
  
  if(!setequal(setdiff(union(names(dimnames(coef)), add_parents), variable), parents))
    stop(paste("ERROR: distribution of variable", variable, "do not match parent set"))

  if(length(dim(coef)) > 1) {
    num_regs = prod(dim(coef)[2:length(dim(coef))])
    if(length(dist$sd) != num_regs)
      stop(paste("ERROR: distribution of variable", variable, "number of regressions and",
                 "number of sd do not match:", num_regs, "against", length(dist$sd)))
  }

  d_levels = dimnames(coef)[2:length(dimnames(coef))]
  comb = (1:nrow(expand.grid(d_levels)))  

  prob = array(c(coef), dim = c(dim(coef)[1], num_regs), dimnames = list(
    variable = dimnames(coef)[[1]],
    comb
  ))
  d_parents = match(names(d_levels), parents)
  g_parents = match(setdiff(dimnames(coef)[[1]], intercept_name), parents)

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

learn_param_mixed = function(dbn, distribution = list(),
                             static_nodes = static_nodes,
                             dynamic_nodes = dynamic_nodes) {
  
  is_clg      <- \(x) is.list(x) && !is.null(x$coef) && !is.null(x$sd)
  is_discrete <- is.array
  is_gaussian <- \(x) is.numeric(x) && is.vector(x)

  nodes = c(static_nodes, dynamic_nodes)
  intercept_std = c(intercept_name, std_name)

  defined_levels = list()
  nodes_info = list()

  # check that distribution contains distributions for each node in the net
  if (!setequal(names(distribution), nodes))
    stop("ERROR: nodes in DBN and variables in distribution do not match")

  # parents may be lagged "D_t-1", map each base variable name to its
  # distribution to look up parent types regardless of slice.
  dist_by_var = list()
  for (n in names(distribution))
    dist_by_var[[get_variable_name(n)]] = distribution[[n]]

  for(variable in names(distribution)) {

    parents = get_parent_set(dbn, variable)
    children = get_children_set(dbn, variable)

    dist = distribution[[variable]]

    if (is_discrete(dist)) {
      # discrete can only have discrete parents
      parent_dists = lapply(parents, \(p) dist_by_var[[get_variable_name(p)]])
      if (!all(sapply(parent_dists, is_discrete)))
        stop(paste("ERROR: found discrete variable", variable,
                   "with non discrete parents"))
      # update levels map and get node_info
      var_name = get_variable_name(variable)
      defined_levels[[var_name]] = dimnames(dist)[[variable]]
      nodes_info[[variable]] = get_node_info_discrete(dist, variable, parents,
                                                      children, defined_levels)
      class(nodes_info[[variable]]) = "dbn.fit.dnode"
    } else if(is_gaussian(dist)) {
      nodes_info[[variable]] = get_node_info_gaussian(dist, variable, 
                                                      parents, children)
      class(nodes_info[[variable]]) = "dbn.fit.gnode"

    } else if(is_clg(dist)) {
      nodes_info[[variable]] = get_node_info_mixed(dist, variable, parents,
                                                   children, defined_levels)
      class(nodes_info[[variable]]) = "dbn.fit.cgnode"
    } else {
      stop(paste("ERROR: found unrecognized node type in distribution for variable", variable, "please use:",
                  "\n1) a multi-dimensional array (CPT, discrete)",
                  "\n2) a named numeric vector (regression coefficients + std, gaussian)",
                  "\n3) a list(coef, sd) (CLG, mixed)"))
    }
  }

  class(nodes_info) = "dbn.fit"
  nodes_info
}


stuff = function() {
  A_t.prob = array(c(1, 0.2, 0.3,
                     1, 0.4, 0.2,
                     1, 0.1, 0.4,
                     1, 0.2, 0.6),
                     dim = c(3, 2, 2),
                     dimnames = list(
                      A_t = c(intercept_name, "A_t-1", "C_t"),
                      D_t = c("yes", "no"),
                      B_t = c("yes", "no")
                     ))

  test = array(c(1, 0.2, 0.3), dim = 3, dimnames = list(
                      A_t = c(intercept_name, "A_t-1", "C_t")
                     ))

  expand.grid(dimnames(A_t.prob)[2:3])

  fitted.prob = array(c(1, 0.2, 0.3,
                        1, 0.4, 0.2,
                        1, 0.1, 0.4,
                        1, 0.2, 0.6),
                        dim = c(3, 4),
                        dimnames = list(
                          A_t = c(intercept_name, "A_t-1", "C_t"),
                          comb = (1:nrow(expand.grid(dimnames(A_t.prob)[2:3])))
                        ))

  A_t.dist = list(coef = A_t.prob, sd = c(.1, .1, .1, .1))
}





