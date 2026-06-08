# DEFINE FUNCTION FOR dbn.fit OBJECT CREATION -> NODES (with childrens, parents and CPTs)

get_static_nodes = function(dbn) {
  bnlearn::node.ordering(from_DBN_to_G_0(dbn))
}

get_dynamic_nodes = function(dbn) {
  g_transition_graph = from_DBN_to_G_transition(dbn)
  quanteda::char_select(bnlearn::node.ordering(g_transition_graph),
                        "*t",
                        valuetype = "glob")
}


#' Function for parameter learning in Dynamic Bayesian Networks
#'
#' Dispatcher for parameter learning. Routes to the discrete or gaussian
#' subroutine based on what the user provides:
#' \itemize{
#'   \item \code{CPTs} non-empty -> discrete network learned from CPTs
#'   \item \code{CPDs} non-empty -> gaussian network learned from CPDs
#'   \item only \code{data}      -> dataset_type(data) decides the subroutine
#' }
#' Exactly one of \code{CPTs}, \code{CPDs} or \code{data} must be non-empty.
#'
#' @param DBN object of class 'dbn'
#' @param CPTs list of multi-dimensional arrays (CPTs for each DBN node, discrete)
#' @param CPDs list of named numeric vectors (regression coefficients + std, gaussian)
#' @param data data.frame object
#' @param replace.unidentifiable If TRUE conditional probabilities for unobserved
#'   parents combinations (unidentifiable parameters) are replaced by uniform
#'   conditional probabilities, if FALSE (default) they are set as NA. Discrete only.
#' 
#' @description If a data.frame is given, the function inspects the class of each column exept
#'   \code{Time} and \code{Sample_id}, if all columns are factor discrete parameter learning is used
#'   if they are all numeric gaussian parameter learning. Mixed case is not supported at the time.
#'
#' @return object of class 'dbn.fit'
#' @export
#'
#' @examples
#' learned_dbn <- dbn.fit(DBN = DBN_example, data = sampling_set)
dbn.fit <- function(DBN, CPTs = list(), CPDs = list(), 
                    data = data.frame(), 
                    replace.unidentifiable = FALSE) {
  if (!class(DBN) == 'dbn')
    stop("ERROR: DBN argument is not of class 'dbn'")
  if (!class(CPTs) == 'list')
    stop("ERROR: CPTs argument is not of class 'list'")
  if (!class(CPDs) == 'list')
    stop("ERROR: CPDs argument is not of class 'list'")
  if (!class(data) == 'data.frame')
    stop("ERROR: data argument is not of class 'data.frame'")
  if (any(is.na(DBN)))
    stop("ERROR: missing data detected")

  # exactly one of CPTs / CPDs / data must be non-empty
  sources = c(length(CPTs) > 0, length(CPDs) > 0, nrow(data) > 0)
  if (sum(sources) == 0)
    stop("ERROR: one of CPTs, CPDs or data must be provided to learn the parameters of the DBN")
  if (sum(sources) > 1)
    stop("ERROR: only one between CPTs, CPDs or data must be defined to learn the parameters of the DBN")

  static_nodes = get_static_nodes(DBN)
  dynamic_nodes = get_dynamic_nodes(DBN)

  # discrete from CPTs
  if (length(CPTs) > 0)
    return(learn_param_d_cpts(DBN, CPTs = CPTs,
                              static_nodes = static_nodes,
                              dynamic_nodes = dynamic_nodes))

  # gaussian from CPDs
  if (length(CPDs) > 0)
    return(learn_param_g_cpds(DBN, CPDs = CPDs,
                              static_nodes = static_nodes,
                              dynamic_nodes = dynamic_nodes))

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


# ---- discrete subroutines -------------------------------------------------

learn_param_d_data = function(dbn, data,
                              static_nodes,
                              dynamic_nodes,
                              replace.unidentifiable = FALSE) {
  vars = colnames(data)[!colnames(data) %in% c('Time', 'Sample_id')]

  # build levels lookup for each node at _0, _t and _t-1
  lvs = list()
  for (i in vars) {
    var   <- gsub(" ", "", paste(i, '_0'))
    var_1 <- gsub(" ", "", paste(i, '_t'))
    var_2 <- gsub(" ", "", paste(i, '_t-1'))
    lvs[[var]]   <- sort(as.array(levels(factor(data[[i]]))))
    lvs[[var_1]] <- sort(as.array(levels(factor(data[[i]]))))
    lvs[[var_2]] <- sort(as.array(levels(factor(data[[i]]))))
  }

  # df_0: time-0 slice with _0 suffix
  df_0 <- data[data$Time == 0,]
  names(df_0) <-
    lapply(names(data), function(x)
      ifelse(x %in% c('Sample_id', 'Time'), x, gsub(" ", "", paste(x, '_0'))))

  # df_transition: full data renamed with _t-1, then leaded into _t
  df_transition <- data
  names(df_transition) <-
    lapply(names(data), function(x)
      ifelse(x %in% c('Sample_id', 'Time'), x, gsub(" ", "", paste(x, '_t-1'))))
  for (k in dynamic_nodes) {
    df_transition <-
      df_transition %>% dplyr::group_by(Sample_id) %>% dplyr::mutate(!!k := dplyr::lead(get(gsub(
        " ", "", paste(k, '-1')
      )), n = 1, default = NA))
  }
  df_transition <-
    df_transition %>% dplyr::ungroup() %>% stats::na.omit()

  dbn_fitted = list()

  for (f in static_nodes) {
    f_star <- strsplit(f, '_')[[1]]
    f_1 <- paste(f_star[1:(length(f_star) - 1)], collapse = '_')
    temp_f <-
      ifelse(f_star[length(f_star)] == '0', 't_0', f_star[length(f_star)])
    parents <- dbn[['nodes']][[f_1]][[temp_f]][['parents']]
    children <- dbn[['nodes']][[f_1]][[temp_f]][['children']]
    if (length(parents) > 0) {
      pr <-
        df_0[, c(rev(parents), f)] %>% dplyr::group_by_all() %>% dplyr::count() %>% dplyr::ungroup() %>% tidyr::complete(!!! rlang::syms(c(rev(parents), f))) %>% replace(is.na(.), 0) %>% dplyr::group_by(dplyr::across(rev(parents))) %>% dplyr::reframe(!!dplyr::sym(f), prob = n/sum(n)) %>% dplyr::ungroup() %>% dplyr::arrange_all(.vars = c(rev(parents), f))
      if (any(is.na(pr$prob))) {
        if (replace.unidentifiable) {
          pr$prob <- replace(pr$prob, is.na(pr$prob), 1/length(lvs[[f]]))
        } else {
          warning("WARNING: Probabilities of the conditioning set equal to 0: Relative frequency is NULL")
        }
      }
      dbn_fitted[[f]] <-
        list(
          node = f,
          parents = parents,
          children = children,
          prob = array(
            pr$prob,
            dim = unname(unlist(lapply(lvs[c(f, parents)], length))),
            dimnames = lvs[c(f, parents)]
          )
        )
    } else {
      prob_vec <-
        unlist(lapply(lvs[[f]], function(x)
          nrow(df_0[df_0[, f] == x,]) / nrow(df_0)))
      dbn_fitted[[f]] <-
        list(
          node = f,
          parents = parents,
          children = children,
          prob = array(
            prob_vec,
            dim = unname(unlist(lapply(lvs[c(f, parents)], length))),
            dimnames = lvs[c(f, parents)]
          )
        )
    }
  }

  for (g in dynamic_nodes) {
    g_star <- strsplit(g, '_')[[1]]
    g_1 <- paste(g_star[1:(length(g_star) - 1)], collapse = '_')
    temp_g <-
      ifelse(g_star[length(g_star)] == '0', 't_0', g_star[length(g_star)])
    parents <- dbn[['nodes']][[g_1]][[temp_g]][['parents']]
    children <- dbn[['nodes']][[g_1]][[temp_g]][['children']]
    if (length(parents) > 0) {
      pr <-
        df_transition[, c(rev(parents), g)] %>% dplyr::group_by_all() %>% dplyr::count() %>% dplyr::ungroup() %>% tidyr::complete(!!! rlang::syms(c(rev(parents), g))) %>% replace(is.na(.), 0) %>% dplyr::group_by(dplyr::across(rev(parents))) %>% dplyr::reframe(!!dplyr::sym(g), prob = n/sum(n)) %>% dplyr::ungroup() %>% dplyr::arrange_all(.vars = c(rev(parents), g))
      if (any(is.na(pr$prob))) {
        if (replace.unidentifiable) {
          pr$prob <- replace(pr$prob, is.na(pr$prob), 1/length(lvs[[g]]))
        } else {
          warning("WARNING: Probabilities of the conditioning set equal to 0: Relative frequency is NULL")
        }
      }
      dbn_fitted[[g]] <-
        list(
          node = g,
          parents = parents,
          children = children,
          prob = array(
            pr$prob,
            dim = unname(unlist(lapply(lvs[c(g, parents)], length))),
            dimnames = lvs[c(g, parents)]
          )
        )
    } else {
      prob_vec <-
        unlist(lapply(lvs[[g]], function(x)
          nrow(df_transition[df_transition[, g] == x,]) / nrow(df_transition)))
      dbn_fitted[[g]] <-
        list(
          node = g,
          parents = parents,
          children = children,
          prob = array(
            prob_vec,
            dim = unname(unlist(lapply(lvs[c(g, parents)], length))),
            dimnames = lvs[c(g, parents)]
          )
        )
    }
  }

  class(dbn_fitted) <- "dbn.fit"
  dbn_fitted
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

  for (i in nodes) {
    CPT <- CPTs[[i]]
    if (!(all(lapply(CPT, class) == 'numeric')))
      stop("ERROR: Probabilities must be numeric")
    if (length(dim(CPT)) > 1) {
      l <- length(dim(CPT))
      idx_target <- which(names(dimnames(CPT)) == i)
      if (!(all(apply(CPT, setdiff(1:l, idx_target), sum) == as.character(1))))
        stop("ERROR: Probabilities for each conditioning set must sum to 1")
    } else {
      if (sum(CPT) != as.character(1))
        stop("ERROR: Probabilities for each conditioning set must sum to 1")
    }

    i_star <- strsplit(i, '_')[[1]]
    i_1 <- paste(i_star[1:(length(i_star) - 1)], collapse = '_')
    temp_i <-
      ifelse(i_star[length(i_star)] == '0', 't_0', i_star[length(i_star)])
    parents <- dbn[['nodes']][[i_1]][[temp_i]][['parents']]
    children <- dbn[['nodes']][[i_1]][[temp_i]][['children']]

    if (!setequal(setdiff(names(dimnames(CPT)), i), parents))
      stop("ERROR: CPTs do not match parents set")

    def_lev <-
      c(dimnames(CPT),
        setNames(dimnames(CPT), array(unlist(
          sapply(names(dimnames(CPT)), function(x) {
            gsub(" ", "", paste(substring(x, 1, (nchar(x) - 2)), '_t'))
          })
        ))),
        setNames(dimnames(CPT), array(unlist(
          sapply(names(dimnames(CPT)), function(x) {
            gsub(" ", "", paste(substring(x, 1, (nchar(x) - 2)), '_t-1'))
          })
        ))))
    int_nodes <- intersect(names(def_lev), names(defined_levels))
    for (c in int_nodes) {
      if (!(setequal(def_lev[[c]], defined_levels[[c]])))
        stop("ERROR: Inconsistency in node's levels")
    }
    defined_levels <-
      c(defined_levels, def_lev[setdiff(names(def_lev), names(defined_levels))])
    nodes_info[[i]] <-
      list(
        node = i,
        parents = parents,
        children = children,
        prob = aperm(CPT, c(i, parents))
      )
  }

  class(nodes_info) <- "dbn.fit"
  nodes_info
}


# ---- gaussian subroutines -------------------------------------------------

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

    # check values in CPD are numeric
    if(!(all(lapply(cpd, class) == 'numeric')))
      stop("ERROR: Probabilities must be numeric")

    # need do check cpd is not list!

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

    nodes_info[[variable]] = list(
      node = variable,
      parents = parents,
      children = children,
      regs = cpd[1:length(cpd) - 1],
      std = cpd[length(cpd)]
    )

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
  }

  class(nodes_info) = "dbn.fit"
  nodes_info
}
