# -------------------------------
#           Variables
# -------------------------------

intercept_name = "(Intercept)"
std_name = "Std (res)"


# -------------------------------
#        Class predicates
# -------------------------------

is.dbn     <- function(x) inherits(x, "dbn")
is.dbn.fit <- function(x) inherits(x, "dbn.fit")
is.bn.fit  <- function(x) inherits(x, "bn.fit")
is.bn      <- function(x) inherits(x, "bn")


is.dnode  <- function(x) inherits(x, c("dbn.fit.dnode",  "bn.fit.dnode"))  && !is.null(x$prob)
is.gnode  <- function(x) inherits(x, c("dbn.fit.gnode",  "bn.fit.gnode"))  && !is.null(x$regs)
is.cgnode <- function(x) inherits(x, c("dbn.fit.cgnode", "bn.fit.cgnode")) && !is.null(x$coefficients)

# -------------------------------
#           Functions
# -------------------------------


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
#' time_lag <- get_variable_time_index("A_t-1")
#' time_lag == 1
get_variable_time_index = function(variable) {
  var_split = strsplit(variable, "_")[[1]]
  tok = var_split[length(var_split)]
  if(tok == 't') return(as.numeric(0))                          # current slice X_t
  if(grepl("^[0-9]+$", tok)) return(as.numeric(tok))            # initial slices X_0, X_1, ...
  if(grepl("^t-[0-9]+$", tok)) return(as.numeric(sub("^t-", "", tok)))  # lagged parents X_t-k
  stop("Error get_variable_time_index: invalid time format")
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

  if(!grepl("^(t-\\d+|[0-9]+|t)$", var_split[length(var_split)])) stop("Error split_variable_name time format not recognized")
  
  # if time is of format like _2 -> t_2
  var_time = ifelse(
    grepl("^[0-9]+$", var_split[length(var_split)]),
        paste0('t_', var_split[length(var_split)]),
        var_split[length(var_split)]
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
  split_variable_name(n)$time
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
  if(is.dbn(dbn)) {
    split = split_variable_name(variable)
    name = split$name
    time = split$time

    dbn$nodes[[name]][[time]]$parents
  } else if(is.dbn.fit(dbn) | is.bn.fit(dbn))
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
  if(is.dbn(dbn)) {
    split = split_variable_name(variable)
    name = split$name
    time = split$time

    dbn$nodes[[name]][[time]]$children
  } else if(is.dbn.fit(dbn) | is.bn.fit(dbn))
    dbn[[variable]]$children
  else
    stop(paste("Get_children_set dbn class not recognized, got",
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
  if (!is.bn.fit(G_transition)) {
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

