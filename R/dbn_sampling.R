#' Get the time point to retrieve the values of variables in the time series
#'
#' @param t time point
#' @param index_n time index of node n
#' @returns time-point in the series used to retrieve the value of the node n
#' @export
#' @examples
#' get_time_point(3, 1)
get_time_point <- function(t, index_n) {
  if (is.character(t) | is.character(index_n)) {
    stop("t and index_n must be numeric")
  }
  # compute t required to get the value of the variable
  # one unit is added cause the series starts from zero
  time_index <- as.integer(t) - as.integer(index_n) + 1
  return(time_index)
}


#' Filter the cpt of node n given the value of its parents
#'
#' @param df_cpt_n a dataframe representing the cpt
#' @param parents_values a list with the parents of n and their values
#' @returns the filtered cpt
#' @export
#' @examples
#' filter_cpt(df_cp_n, parents_values)
filter_cpt<- function(df_cpt_n, parents_values) {
  # Check if df_cpt_n is a dataframe
  if (!is.data.frame(df_cpt_n)) {
    stop("df_cpt_n must be a dataframe!")
  }
  
  #Check if parents_values is a list
  if (!is.list(parents_values)) {
    stop("parents_values must be a list!")
  }
  
  # Create a logical index vector for all conditions at once
  filter_index <- rep(TRUE, nrow(df_cpt_n))
  
  for (key in names(parents_values)) {
    filter_index <- filter_index & (df_cpt_n[[key]] == parents_values[[key]])
  }
  
  # Apply the combined logical index to filter the dataframe in one go
  filtered_df <- df_cpt_n[filter_index, ]
  
  return(filtered_df)
}


# builds a cache of variables information:
# given variable A_t-1 in a network stores:
# - name A
# - time t-1
# - time_index 1
# - parent set
build_cache <- function(bn_0, bn_transition, dbn_type = "gaussian") {
  cache_names <- new.env(parent = emptyenv())

  add_value <- function(variable, dbn, cache_names) {
    if (!exists(variable, envir = cache_names)) {
      name_time <- split_variable_name(variable)
      parents <- get_parent_set(dbn, variable)
      time_index <- get_variable_time_index(variable)
      result <- list(
        name = name_time$name,
        time = name_time$time,
        parents = parents,
        time_index = time_index
      )
      if (dbn_type == "discrete" && time_index == 0) {
        if(variable %in% names(bn_0)){
          df_cp_n <- data.frame(as.table(bn_0[[variable]][["prob"]]))
        } else {
          df_cp_n <- data.frame(as.table(bn_transition[[variable]][["prob"]]))
        }
        colnames(df_cp_n) <- gsub("\\.", "-", colnames(df_cp_n))
        result[["cp"]] = df_cp_n
      }
      assign(variable, result, envir = cache_names)
    }
  }

  for (variable in bnlearn::node.ordering(bn_0)) {
    add_value(variable, bn_0, cache_names)
    for (parent in get(variable, envir = cache_names)$parents) {
      add_value(parent, bn_0, cache_names)
    }
  }
  for (variable in bnlearn::node.ordering(bn_transition)) {
    add_value(variable, bn_transition, cache_names)
    for (parent in get(variable, envir = cache_names)$parents) {
      add_value(parent, bn_transition, cache_names)
    }
  }
  cache_names
}

sample_variable_gaussian <- function(ts_dict, bn_0t, variable, obs, t, time,
                                      cache_names, sd_sampled) {
  # get variable name (no time indx) and parent set
  parents <- get(variable, envir = cache_names)$parents
  # first element is 1 for product of intercept
  parents_values <- double(length = length(parents) + 1)
  parents_values[1] <- 1
  parents_idx <- c(1)

  # extract parents values
  for (i in seq_along(parents)) {
    parent_infos <- get(parents[i], envir = cache_names)
    parent_name <- parent_infos$name
    # get the time point for the parent in order to retrieve the
    # value from the time series
    time_index <- parent_infos$time_index
    time_point <- get_time_point(t, time_index)
    if (time_point > 0) {
      parents_values[i + 1] <- ts_dict[[parent_name]][
        (obs - 1) * (time + 1) + time_point
      ]
      parents_idx <- c(parents_idx, i + 1)
    } else {
      parents_values[i + 1] <- 0
    }
  }

  sum(parents_values[parents_idx] * bn_0t[[variable]]$regs[parents_idx]) + sd_sampled
}

sample_variable_discrete <- function(ts_dict, bn_0t, variable, obs, t, time,
                                      cache_names, sd_sampled) {
  # get variable name (no time indx) and parent set
  parents <- get(variable, envir = cache_names)$parents
  # first element is 1 for product of intercept
  parents_values <- list()
  parents_idx <- c()

  # extract parents values
  for (i in seq_along(parents)) {
    parent_infos <- get(parents[i], envir = cache_names)
    parent_name <- parent_infos$name
    # get the time point for the parent in order to retrieve the
    # value from the time series
    time_index <- parent_infos$time_index
    time_point <- get_time_point(t, time_index)
    if (time_point > 0) {
      # get correct name for the names list
      name = node_id(parent_name, parent_infos$time)
      parents_values[name] <- ts_dict[[parent_name]][
        (obs - 1) * (time + 1) + time_point
      ]
      parents_idx <- c(parents_idx, name)
    }
  }

  parents_values = parents_values[parents_idx]

  df_cp_n = get(variable, envir = cache_names)$cp
  filtered_cpt_n <- filter_cpt(df_cp_n, parents_values)
  #sampling n using frequencies
  as.character(sample(filtered_cpt_n[[variable]], size = 1, 
                      prob = filtered_cpt_n[["Freq"]]))
}

remove_prev_time_from_bn_fit = function(bn) {
  class(bn) = "list"
  for(variable in names(bn)) {
    parents = bn[[variable]]$parents
    bn[[variable]]$parents = parents[!grepl("t-[0-9]+$", parents)]
    
    children = bn[[variable]]$children
    bn[[variable]]$children = children[!grepl("t-[0-9]+$", children)]
  }
  class(bn) = "bn.fit"
  bn
}


#' Generate a sampling dataset
#'
#' @param fitted_DBN an object of class 'dbn'
#' @param n_samples number of samples
#' @param max_time time series length
#' @returns the generated dataframe
#' @export
#' @examples
#' dbn.sampling(DBN_example, N_samples, Time)
dbn.sampling.R <- function(fitted_dbn, n_samples, max_time) {
 if (is.character(max_time)) {
    stop("Time must be an integer!")
  }
  if (max_time < 1){
    stop("Time must be greater than 0!")
  }
  if (is.character(n_samples)) {
    stop("N_samples must be an integer!")
  }
  if (n_samples < 1){
    stop("N_samples must be greater than 0!")
  }
  if (class(fitted_dbn) != "dbn.fit") {
    stop("fitted_DBN must be a dbn.fit object")
  }

  dbn_type = dbn_type(fitted_dbn)
  sampling_fun <- switch(
    dbn_type,
    "discrete" = sample_variable_discrete,
    "gaussian" = sample_variable_gaussian,
    stop("Invalid dbn_type")
  )

  bn_0 <- from_fitted_DBN_to_fitted_G_0(fitted_dbn)
  bn_transition <- from_fitted_DBN_to_fitted_G_transition(fitted_dbn)
  # first remove t-1 parents in bn_tranistion
  # need to to this to get the correct node.ordering
  bn_transition_2 <- remove_prev_time_from_bn_fit(bn_transition)
  # get node ordering (only _t)
  nodes_t <- get_nodes_t(bn_transition_2)

  n_samples <- as.integer(n_samples)
  max_time <- as.integer(max_time)
  timeseries_dict <- list(Time = c())

  # build cache of node info
  cache_names <- build_cache(bn_0, bn_transition, dbn_type = dbn_type)
  sd = 0
  for (observation in seq(n_samples)) {
    # start with t = 0
    timeseries_dict[["Time"]] <- c(timeseries_dict[["Time"]], 0)
    timeseries_dict[["Sample_id"]] <- c(
      timeseries_dict[["Sample_id"]],
      paste("sample", observation, sep = "")
    )

    for (variable in bnlearn::node.ordering(bn_0)) {
      # current time is 0
      var_name <- get(variable, envir = cache_names)$name

      # discrete don't need std
      if(dbn_type == "gaussian")
        sd <- rnorm(1, 0, bn_0[[variable]]$std)

      v <- sampling_fun(
        timeseries_dict, bn_0,
        variable, observation, 0, max_time, cache_names,
        sd
      )

      timeseries_dict[[var_name]] <- append(
        timeseries_dict[[var_name]], v
      )
    }

    if(dbn_type == "gaussian") {
      # sampling normals with sds for the entire traj for every variables
      # saves a lot of time
      sds <- list()

      for (variable in nodes_t) {
        sds[[variable]] <- rnorm(max_time, 0, bn_transition[[variable]]$std)
        var_name <- get(variable, envir = cache_names)$name
        timeseries_dict[[var_name]] <- append(
          timeseries_dict[[var_name]], double(length = max_time)
        )
      }
    }
    
    # now iterate over time
    for (t in 1:max_time) {
      for (variable in nodes_t) {
        var_name <- get(variable, envir = cache_names)$name
        v <- sampling_fun(
          timeseries_dict, bn_transition,
          variable, observation, t, max_time,
          cache_names, sds[[variable]][t]
        )

        timeseries_dict[[var_name]][(observation - 1) * (max_time + 1) + t + 1] <- v
      }
    }

    # add time
    timeseries_dict[["Time"]] <-
      c(timeseries_dict[["Time"]], seq(max_time))
    timeseries_dict[["Sample_id"]] <-
      c(
        timeseries_dict[["Sample_id"]],
        rep(paste("sample", observation, sep = ""), max_time)
      )
  }
  df_timeseries <- data.frame(timeseries_dict)
  return(df_timeseries)
}



#' Forecast next values from observations
#'
#' @param fitted_DBN an object of class 'dbn'
#' @param observations a list (or data.frame) containing the observations for each variable.
#'   data.frames are coerced to a list internally; the forecast trajectory is grown
#'   in-place via append(), which is a list operation.
#' @param timepoints number of timepoints to be forecasted
#' @returns a dataframe containing the forecasting
#' @export
#' @examples
#' dbn.forecasting(fitted_DBN, observations, timepoints)
dbn.forecasting <- function(dbn, observations, timepoints) {
   if (is.character(timepoints)) {
    stop("timepoints must be an integer!")
  }
  if (timepoints < 1){
    stop("timepoints must be greater than 0!")
  }
  #check that observation is a list
  if (class(dbn) != "dbn.fit") {
    stop("fitted_DBN must be a dbn.fit object")
  }

  # accept a data.frame for convenience (dbn.sampling returns one): internally
  # we grow per-variable vectors with append(), which fights the nrow constraint
  # of data.frames.
  if (is.data.frame(observations)) {
    observations <- as.list(observations)
  }

  if (!is.list(observations)) {
    stop("Parameter type not valid observations must be a list or data.frame")
  }

  markov_order = get_max_mo_dbn_fit(dbn)
  
  # maybe here use priors instead of erroring out
  if(length(observations[["Time"]]) < markov_order) {
    stop(paste("Dbn_forecasting, provided", length(observations[["Time"]]),
               "observations for a dbn of markov order:", markov_order))
  }
  
  dbn_type = dbn_type(dbn)
  sampling_fun <- switch(
    dbn_type,
    "discrete" = sample_variable_discrete,
    "gaussian" = sample_variable_gaussian,
    stop("Invalid dbn_type")
  )

  bn_0 <- from_fitted_DBN_to_fitted_G_0(dbn)
  bn_transition <- from_fitted_DBN_to_fitted_G_transition(dbn)
  # first remove t-1 parents in bn_tranistion
  # need to to this to get the correct node.ordering
  # get node ordering (only _t)
  nodes_t <- get_nodes_t(remove_prev_time_from_bn_fit(bn_transition))
  # build cache of node info
  cache_names <- build_cache(bn_0, bn_transition, dbn_type = dbn_type)


  # make observations as long as the forecasting horizon
  observations$Sample_id = append(observations$Sample_id, 
                                  rep(observations$Sample_id[1], timepoints))
  observations$Time = append(observations$Time, markov_order:(markov_order + timepoints - 1))
  for(node in nodes_t) {
    var_name = get(node, envir = cache_names)$name
    observations[[var_name]] = append(observations[[var_name]], 
                                      double(length = timepoints))
  }
  
  # predict next values
  for(t in 1:timepoints) {
    for(variable in nodes_t) {
      var_name = get(variable, envir = cache_names)$name
      v = sampling_fun(observations, bn_transition,
                            variable, 1, markov_order + t - 1, 0,
                            cache_names, 0)

      observations[[var_name]][markov_order + t] = v
    }
  }
  
  observations
}