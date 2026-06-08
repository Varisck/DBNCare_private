# -------------------------------
#           Variables
# -------------------------------

intercept_name = "(Intercept)"
std_name = "Std (res)"


# -------------------------------
#           Functions
# -------------------------------


#' rename_nodes_unroll
#' @description
#' Renames a list of nodes, making nodes names compliant to the dbn unrolled structure
#' 
#' @param time_slice an integer representing a time slice
#' @param nodes_names a list of nodes names 
#'
#' @export
#'
rename_nodes_unroll <- function(time_slice, nodes_names) {
  renamed_nodes <- c()
  for (name_n in nodes_names) {
    #does the node end with t-1
    end_t_1 <- grepl("_t-1$", name_n)
    end_with_t <- grepl("_t$", name_n)
    if (end_t_1) {
      previous_time_slice <- time_slice - 1
      root_node_name <- get_generic_node_name_rex(name_n)
      new_node_name <-
        paste(root_node_name,
              "_",
              as.character(previous_time_slice),
              sep = "")
    }
    if (end_with_t) {
      root_node_name <- get_generic_node_name_rex(name_n)
      new_node_name <-
        paste(root_node_name, "_", as.character(time_slice), sep = "")
    }
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
#' @export
#'
get_node_edges <- function(dbn_transition, node, time_slice) {
  if (class(dbn_transition) != "bn.fit") {
    stop("dbn must be a bn.fit object!")
  }
  if (is.character(node) == FALSE) {
    stop("node must be a character!")
  }
  if (is.numeric(time_slice)==FALSE){
    stop("time slice must be an integer!")
  }
  node_edges <- c()
  root_node_name_n_t <- get_generic_node_name_rex(node)
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

#' get_unrolled_dbn
#' 
#' @description 
#' this function generates an unrolled dynamic bayesian network
#' @param dbn_fitted a dbn.fit class object
#' @param slices number of time-slice of the unrolled network
#'
#' @return a bn.fit object 
#' @export
#'
#' @examples get_unrolled_dbn(my_fitted_dbn, 4)
get_unrolled_dbn <- function(dbn_fitted, slices) {
  if (is.character(slices)) {
    stop("slices must be an integer!")
  }
  if (slices < 1){
    stop("slices must be greater than 0!")
  }
  if (class(dbn_fitted) != "dbn.fit") {
    stop("dbn_fitted must be a dbn.fit object")
  }
  nodes_bn <- c()
  edges_bn <- c()
  cpt_bn <- list()
  my_dbn_transition <-
    from_fitted_DBN_to_fitted_G_transition(dbn_fitted)
  nodes_in_transition <- get_nodes_t(my_dbn_transition)
  my_dbn_g_0 <- from_fitted_DBN_to_fitted_G_0(dbn_fitted)
  nodes_in_0 <- bnlearn::node.ordering(my_dbn_g_0)
  #iterate on nodes and get,for each one, edges and probability tables
  for (n_0 in nodes_in_0) {
    #add the node in bn_nodes
    nodes_bn <- c(nodes_bn, n_0)
    cpt_bn[[n_0]] <- my_dbn_g_0[[n_0]][["prob"]]
    parents_n_0 <- my_dbn_g_0[[n_0]][["parents"]]
    if (length(parents_n_0) > 0) {
      for (p_n_0 in parents_n_0) {
        #using "from" "to" standard
        edge_n_0 <- c(p_n_0, n_0)
        edges_bn <- c(edges_bn, edge_n_0)
      }
    }
  }
  for (time_slice in seq(slices)) {
    for (n_t in nodes_in_transition) {
      cpt_n_t <- my_dbn_transition[[n_t]][["prob"]]
      names_cpt_n_t <- names(dimnames(cpt_n_t))
      #rename each single node if node has t rename t with time_slice
      #if name has t-1 subtract 1 to time_slice.
      renamed_nodes <- rename_nodes_unroll(time_slice, names_cpt_n_t)
      root_node_name_n_t <- get_generic_node_name_rex(n_t)
      new_node_name_n_t <-
        paste(root_node_name_n_t, "_", as.character(time_slice), sep = "")
      names(dimnames(cpt_n_t)) <- renamed_nodes
      cpt_bn[[new_node_name_n_t]] <- cpt_n_t
      #adding node
      nodes_bn <- c(nodes_bn, new_node_name_n_t)
      #adding edges
      edges_time_slice <-
        get_node_edges(my_dbn_transition, n_t, time_slice)
      if (length(edges_time_slice) > 0) {
        edges_bn <- c(edges_bn, edges_time_slice)
        
      }
    }
  }
  dag_unrolled <- bnlearn::empty.graph(nodes = nodes_bn)
  arc_unrolled.set <- matrix(
    edges_bn,
    byrow = TRUE,
    ncol = 2,
    dimnames = list(NULL, c("from", "to"))
  )
  bnlearn::arcs(dag_unrolled) <- arc_unrolled.set
  dbn_unrolled <- bnlearn::custom.fit(dag_unrolled, cpt_bn)
  return(dbn_unrolled)
}


# takes variables col names and add postfix at the end
concat_name_post = function(name, postfix, 
                            exclude_cols = c("Sample_id", "Time")) {
  name = ifelse(name %in% exclude_cols, 
                name, 
                gsub(" ", "", paste(name, postfix)))
}

#' Get the time index of a node
#'
#' @param n name of node 
#' 
#' @returns time index of a node n
#' @export
#' 
#' @examples
#' time_lag <- get_index_regular_expression("A_t-1")
#' time_lag == 1
get_variable_time_index = function(variable) {
  var_split = strsplit(variable, "_")[[1]]
  if(var_split[length(var_split)] == '0' ||
     var_split[length(var_split)] == 't')
    return(as.numeric(0))
  if(!grepl("^t-\\d+$", var_split[length(var_split)]))
    stop("Error get_variable_time_index: invalid time format")
  time = strsplit(var_split[length(var_split)], "t-")[[1]]
  as.numeric(time[length(time)])
}


#' Split the name of a variable in a list with name and time
#'
#' @param n name of node 
#' 
#' @returns list containing name and time
#' @export
#' 
#' @examples
#' var_split <- split_variable_name("A_t-1")
#' var_split == ("A", "t-1")
split_variable_name = function(variable) {
  var_split = strsplit(variable, "_")[[1]]
  if(length(var_split) == 1) stop("Error split_variable_name input is not a recognized variable") 
  var_name = paste(var_split[1:(length(var_split) - 1)], collapse = '_')
  if(!grepl("^(t-\\d+|0|t)$", var_split[length(var_split)])) stop("Error split_variable_name time format not recognized")
  var_time = ifelse(
    var_split[length(var_split)] == '0', 't_0', var_split[length(var_split)]
  )
  list("name" = var_name, "time" = var_time)
}

#' Get the name of a variable
#'
#' @param n variable name like (A_t-1)
#' 
#' @returns characther the name of the variable (A)
#' @export
#' 
#' @details equivalent to \code{split_variable_name(n)$name}
#' 
#' @examples
#' name <- split_variable_name("A_t-1")
#' name == "A"
get_variable_name <- function(n) {
  split_variable_name(n)$name
}

#' Get the time of a variable
#'
#' @param n variable name like (A_0)
#' 
#' @returns characther the time of the variable (t_0)
#' @export
#' 
#' @details equivalent to \code{split_variable_name(n)$time}
#' 
#' @examples
#' time <- split_variable_time("A_0")
#' time == "t_0"
get_variable_time <- function(n) {
  split_variable_name(n)$name
}

#' Get the parents of a node
#'
#' @param G A DBN/dbn.fit/bn.fit object
#' @param n name of node n
#' 
#' @returns parents of the node n.
#' @export
#' 
#' @examples
#' get_parent_set(G, "A_t")
get_parent_set = function(dbn, variable) {
  if(class(dbn) == "dbn") {
    split = split_variable_name(variable)
    name = split$name
    time = split$time
    
    dbn$nodes[[name]][[time]]$parents
  } else if(class(dbn) == "dbn.fit" | 
            class(dbn) == "bn.fit")
    dbn[[variable]]$parents
  else
    stop(paste("Get_parent_set dbn class not recognized, got",
               class(dbn)))
}

#' Get the children of a node
#'
#' @param G A DBN/dbn.fit/bn.fit object
#' @param n name of node n
#' 
#' @returns childrens of the node n.
#' @export
#' 
#' @examples
#' get_children_set(G, "A_t")
get_children_set = function(dbn, variable) {
  if(class(dbn) == "dbn") {
    split = split_variable_name(variable)
    name = split$name
    time = split$time
    
    dbn$nodes[[name]][[time]]$children
  } else if(class(dbn) == "dbn.fit" | 
            class(dbn) == "bn.fit")
    dbn[[variable]]$children
  else
    stop(paste("Get_parent_set dbn class not recognized, got",
               class(dbn)))
}

#' Return the variables indexes by t in G_transition
#'
#' @param G_transition a G_transition graph.
#' 
#' @returns nodes indexed by t only e.g. A_t, B_t.
#' @export
#' 
#' @examples
#' get_nodes_t(G_transition)
get_nodes_t <- function(G_transition) {
  if (class(G_transition) != "bn.fit") {
    stop("G_transition must be a bn.fit object!")
  }
  nodes_dbn <- bnlearn::node.ordering(G_transition)
  nodes_t <- c()
  for (n in nodes_dbn) {
    ends_with_t <- substr(n, nchar(n), nchar(n)) == "t"
    if (ends_with_t == TRUE) {
      nodes_t <- c(nodes_t, n)
    }
  }
  return(nodes_t)
}


#' Return the maximum markov order of a dbn.fit object
#'
#' @param dbn an object of class dbn.fit
#' 
#' @returns markov order of the dbn
#' @export
#' 
#' @examples
#' mo <- get_max_mo_dbn_fit(fit)
get_max_mo_dbn_fit = function(dbn) {
  if(class(dbn) != "dbn.fit")
    stop("get_max_mo_dbn_fit input object is not a dbn.fit")
  nodes_t = get_nodes_t(remove_prev_time_from_bn_fit(dbn))
  mx = 0
  for(n in nodes_t) {
    if("prob" %in% names(dbn[[n]])) {
      regnames = names(dimnames(dbn[[n]]$prob))
    } else {
      regnames = names(dbn[[n]]$regs)
      if(length(regnames) > 1) regnames = regnames[2:length(regnames)]
      else regnames = character(0)
    }
    lagged = regnames[grepl("t-[0-9]+$", regnames)]
    if(length(lagged) > 0) {
      tl = as.double(sub(".*t-([0-9]+)$", "\\1", lagged))
      mx = max(mx, max(tl))
    }
  }
  mx
}



#' Return type of dbn
#'
#' @param dbn object of type dbn.fit
#' @returns string "discrete", "gaussian" or "mixed"
#'
#' @details if the type is not recognized raise error
#' 
#' @export
#' @examples
#' type <- dbn_type(dbn)
dbn_type = function(dbn) {
  if(class(dbn) != "dbn.fit") {
    stop("dbn_type input dbn must be object of class dbn.fit")
  }

  nodes = names(dbn)[grepl("_t$", names(dbn))]

  discrete = all(sapply(nodes, \(n) {
    !is.null(dbn[[n]]$prob)
  }))

  if(all(discrete)) return("discrete")

  gaussian = all(sapply(nodes, \(n) {
    !is.null(dbn[[n]]$regs)
  }))

  if(all(gaussian)) return("gaussian")

  if(any(discrete) && any(gaussian)) return("mixed")

  stop("dbn_type dbn not recognized as neither discrete, gaussian nor mixed")
}


#' Return type of a dataset
#'
#' Inspects the variable columns of a dataset (every column except
#' \code{Sample_id} and \code{Time}) and reports whether the data is
#' discrete, gaussian (continuous numeric) or mixed.
#'
#' @param data a data.frame
#' @returns string "discrete", "gaussian" or "mixed"
#'
#' @details factor and character columns count as discrete; numeric
#'   columns count as gaussian. Raises an error if no variable columns
#'   are found or if any column is of an unsupported type.
#'
#' @export
#' @examples
#' type <- data_type(data)
dataset_type = function(data) {
  if(!is.data.frame(data)) {
    stop("data_type input data must be a data.frame")
  }

  vars = setdiff(colnames(data), c("Sample_id", "Time"))
  if(length(vars) == 0) {
    stop("data_type no variable columns found (only Sample_id/Time present)")
  }

  discrete = sapply(vars, \(v) is.factor(data[[v]]) || is.character(data[[v]]))
  gaussian = sapply(vars, \(v) is.numeric(data[[v]]))

  unsupported = !(discrete | gaussian)
  if(any(unsupported)) {
    stop(paste("data_type unsupported column type(s):",
               paste(vars[unsupported], collapse = ", ")))
  }

  if(all(discrete)) return("discrete")
  if(all(gaussian)) return("gaussian")
  "mixed"
}

#' Given a data.frame and a markov_order builds a transition data.frame
#'
#' @param data a data.frame
#' @param markov_order = 1 markov order of the data frame
#' @param separator = "-" separator to use when creating variable X_t-k
#' @returns data.frame shifted df
#'
#' @details for each variable excluding Sample_id and Time creates 
#'  variables X_t, ..., X_t-k where k is the markov_order
#'
#' @export
#' @examples
#' shifted_df <- build_shifted_df(data, 4)
build_shifted_df = function(data, markov_order = 1, separator = "-"){
  names(data) = lapply(names(data), concat_name_post, postfix = "_t")
  
  vars = setdiff(names(data), c("Sample_id", "Time"))
  
  for(shift in 1:markov_order) {
    for(variable in vars) {
      data = data %>% dplyr::group_by(Sample_id) %>%
        dplyr::mutate(!!paste0(variable, separator, shift) :=
                 dplyr::lag(.data[[variable]], n = shift, default = NA))
    }
  }
  data = data %>% dplyr::ungroup() %>% stats::na.omit()
  as.data.frame(data)
}





