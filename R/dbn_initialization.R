# CREATING AN EMPTY DBN GIVEN THE SETS OF STATIC NODES AND DYNAMIC NODES AND THE MARKOVIAN ORDER

#' Function to create a Dynamic Bayesian Network
#'
#' @param static_nodes character vector listing the time-independent nodes
#' @param dynamic_nodes character vector listing the time-dependent nodes
#' @param markov_order integer value (> 1): markovian order of the process
#'
#' @return object of class 'dbn'
#' @export
#'
#' @examples
#' DBN_example <- empty.dbn(static_nodes = c("K"), dynamic_nodes = c("A", "S", "E", "O", "R", "T", "B"), markov_order = 1)
#' class(DBN_example) # 'dbn'
empty.dbn <- function(static_nodes = c(), dynamic_nodes, markov_order) {
  if (is.null(static_nodes)) {
    if (length(intersect(static_nodes, dynamic_nodes)) != 0) {
      stop('A node could not be both Static and Dynamic!!!')
    }
  }
  if (markov_order < 1) {
    stop('Markov order of the process must be 1 or higher!!!')
  }

  blacklist = blacklist_g0_gt(sapply(dynamic_nodes, concat_name_post, postfix = "_0"),
                              sapply(dynamic_nodes, concat_name_post, postfix = "_t"),
                              markov_order = markov_order,
                              static_nodes = sapply(static_nodes, concat_name_post, postfix = "_0"))

  DBN <- list(
    learning = list(
      whitelist = NULL,
      blacklist = blacklist,
      test = 'none',
      ntests = 0,
      algo = 'empty',
      args = list()
    ),
    markov_order = markov_order,
    arcs = matrix(
      ncol = 2,
      nrow = 0,
      byrow = TRUE,
      dimnames = list(character(0), c("from", "to"))
    ),
    nodes = list()
  )

  dyn_details <- list(type = 'Dynamic', 't_0' = empty_node_details())
  for (i in 0:markov_order) {
    key <- if (i == 0) 't' else paste0('t-', i)
    dyn_details[[key]] <- empty_node_details()
  }
  for (i in dynamic_nodes) {
    DBN$nodes[[i]] <- dyn_details
  }
  for (i in static_nodes) {
    DBN$nodes[[i]][['type']] <- 'Static'
    DBN$nodes[[i]][['t_0']] <- empty_node_details()
  }
  class(DBN) <- "dbn"
  DBN
}

# ADDING NODE TO DBN

#' Function for node addiction in Dynamic Bayesian Network
#'
#' @param DBN object of class 'dbn'
#' @param node object of class 'character': name of the node to add
#' @param type 'Dynamic' (default) or 'Static' node
#'
#' @return object of class 'dbn'
#' @export
#'
#' @examples
#' DBN_example <- add.node.dbn(DBN=DBN_example, node="F", type='Dynamic')
#' DBN_example <- add.node.dbn(DBN=DBN_example, node="K", type='Static')
add.node.dbn <- function(DBN, node, type = 'Dynamic') {
  if (!class(DBN) == 'dbn')
    stop("ERROR: DBN argument is not of class 'dbn'")
  if (!is.character(node))
    stop("ERROR: node is not a character")
  if (node %in% names(DBN$nodes)) {
    stop("ERROR: node name already exists")
  }
  if (type == 'Dynamic') {
    dyn_details <- list(type = 'Dynamic', 't_0' = empty_node_details())
    for (i in 0:DBN$markov_order) {
      key <- if (i == 0) 't' else paste0('t-', i)
      dyn_details[[key]] <- empty_node_details()
    }
    DBN$nodes[[node]] <- dyn_details
  }
  else if (type == 'Static') {
    DBN$nodes[[node]][['type']] <- 'Static'
    DBN$nodes[[node]][['t_0']] <- empty_node_details()
  }
  else{
    stop("ERROR: type must be 'Dynamic' or 'Static")
  }
  DBN
}


# ADDING ARC FUNCTION

# Validates an add.arc.dbn call and returns the canonical from/to ids.
# Stops with the appropriate error otherwise.
validate_arc_DBN <- function(DBN, from, to) {
  if (!class(DBN) == 'dbn')
    stop("ERROR: DBN argument is not of class 'dbn'")
  if (!is.character(from))
    stop("ERROR: from is not a character")
  if (!is.character(to))
    stop("ERROR: to is not a character")
  if (!(length(from) == 2 & length(to) == 2)) {
    stop("ERROR: from and to must be vectors of type (variable, time)")
  }
  if (!(from[1] %in% names(DBN$nodes) &
        to[1] %in% names(DBN$nodes))) {
    stop("ERROR: Defined node not in DBN!!!")
  }
  if ((from[1] == to[1] & from[2] == to[2])) {
    stop("ERROR: The arc defines a loop!!!")
  }

  if (to[2] == 't_0') {
    if (from[2] != 't_0') {
      stop("ERROR: Prior arcs must be between nodes at 't_0'")
    }
    from_ <- node_id(from[1], 't_0')
    to_ <- node_id(to[1], 't_0')
  }
  else if (to[2] == 't') {
    if (!is_valid_parent_time(from[2])) {
      stop("ERROR: Parent nodes in DBN must be at time 't_0', 't' or 't-k' (with k <= Markov Order)")
    }
    if (from[2] != 't' &
        as.numeric(substring(from[2], 3)) > DBN$markov_order) {
      stop("ERROR: Transition arcs could not be of order higher than DBN's Markov Order")
    }
    from_ <- node_id(from[1], from[2])
    to_ <- node_id(to[1], to[2])
  }
  else{
    stop("ERROR: Children nodes in DBN must be at time 't' or 't_0'")
  }

  if ((DBN[['nodes']][[from[1]]][['type']] == 'Dynamic' &
       DBN[['nodes']][[to[1]]][['type']] == 'Static')) {
    stop("ERROR: A Dynamic Node could not be parent of a Static Node!!!")
  }

  if ((from[2] != 't_0' &
       DBN[['nodes']][[from[1]]][['type']] == 'Static')) {
    stop(
      "ERROR: Static Node have not order higher than 0, thus they are not included in transition arcs!!!"
    )
  }

  if ((to[2] != 't_0' &
       DBN[['nodes']][[to[1]]][['type']] == 'Static')) {
    stop(
      "ERROR: Static Node have not order higher than 0, thus they are not included in transition arcs!!!"
    )
  }

  list(from_id = from_, to_id = to_)
}

# Walks the existing arcs to detect whether adding from -> to would create a
# cycle. Stops with the cycle error if so. Mirrors the behavior of the original
# inline check (which has no visited-set; preserved unchanged).
check_arc_cycle_DBN <- function(DBN, from, to) {
  if (identical(DBN[['nodes']][[from[1]]][[from[2]]][['parents']], character(0)) ||
      identical(DBN[['nodes']][[to[1]]][[to[2]]][['children']], character(0))) {
    return(invisible(NULL))
  }
  chi <- DBN[['nodes']][[to[1]]][[to[2]]][['children']]
  s <- TRUE
  while (s) {
    parents_from <- DBN[['nodes']][[from[1]]][[from[2]]][['parents']]
    overlap <- intersect(parents_from, chi)
    if (!identical(overlap, character(0)) & !is.null(overlap)) {
      stop("ERROR: The arc create a cycle!!!")
    }
    chi1 <- c()
    for (j in chi) {
      pj <- parse_node_id(j)
      chi1 <- c(chi1, DBN[['nodes']][[pj$name]][[pj$time]][['children']])
    }
    chi <- chi1
    if (identical(chi, character(0)) | is.null(chi)) {
      s <- FALSE
    }
  }
  invisible(NULL)
}

# Records the arc in the DBN: updates children/parents/nbr/mb on the involved
# nodes, propagates Markov-blanket entries from existing parents of `to`, and
# appends the row to DBN$arcs.
record_arc_DBN <- function(DBN, from, to, from_id, to_id) {
  if (to_id %in% DBN[['nodes']][[from[1]]][[from[2]]][['children']]) {
    return(DBN)
  }
  DBN[['nodes']][[from[1]]][[from[2]]][['children']] <-
    c(DBN[['nodes']][[from[1]]][[from[2]]][['children']], to_id)
  DBN[['nodes']][[to[1]]][[to[2]]][['parents']] <-
    c(DBN[['nodes']][[to[1]]][[to[2]]][['parents']], from_id)
  if (!to_id %in% DBN[['nodes']][[from[1]]][[from[2]]][['nbr']]) {
    DBN[['nodes']][[from[1]]][[from[2]]][['nbr']] <-
      c(DBN[['nodes']][[from[1]]][[from[2]]][['nbr']], to_id)
    DBN[['nodes']][[to[1]]][[to[2]]][['nbr']] <-
      c(DBN[['nodes']][[to[1]]][[to[2]]][['nbr']], from_id)
  }
  if (!to_id %in% DBN[['nodes']][[from[1]]][[from[2]]][['mb']]) {
    DBN[['nodes']][[from[1]]][[from[2]]][['mb']] <-
      c(DBN[['nodes']][[from[1]]][[from[2]]][['mb']], to_id)
    DBN[['nodes']][[to[1]]][[to[2]]][['mb']] <-
      c(DBN[['nodes']][[to[1]]][[to[2]]][['mb']], from_id)
  }
  if (!is.na(DBN$arcs[DBN$arcs[, 'to'] == to_id, ]['from'])) {
    for (i in DBN$arcs[DBN$arcs[, 'to'] == to_id, ]['from']) {
      pi <- parse_node_id(i)
      if (!i %in% DBN[['nodes']][[from[1]]][[pi$time]][['mb']]) {
        DBN[['nodes']][[from[1]]][[from[2]]][['mb']] <-
          c(DBN[['nodes']][[from[1]]][[from[2]]][['mb']], i)
        DBN[['nodes']][[pi$name]][[pi$time]][['mb']] <-
          c(DBN[['nodes']][[pi$name]][[pi$time]][['mb']], from_id)
      }
    }
  }
  DBN[['arcs']] <- rbind(DBN$arcs, c(from_id, to_id))
  DBN
}

#' Function for arc addiction in Dynamic Bayesian Network
#'
#' @param DBN object of class 'dbn'
#' @param from two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param to two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param cycle_OK boolean value, TRUE if cycle are admitted (default)
#'
#' @return object of class 'dbn'
#' @export
#'
#' @examples
#' DBN_example <- add.arc.dbn(DBN=DBN_example,from=c('A','t_0'),to=c('R','t_0'))
#' DBN_example <- add.arc.dbn(DBN=DBN_example,from=c('S','t'),to=c('O','t'))
#' DBN_example <- add.arc.dbn(DBN=DBN_example,from=c('S','t-1'),to=c('S','t'))
add.arc.dbn <- function(DBN, from, to, cycle_OK = TRUE) {
  ids <- validate_arc_DBN(DBN, from, to)
  if (cycle_OK == FALSE) {
    check_arc_cycle_DBN(DBN, from, to)
  }
  record_arc_DBN(DBN, from, to, ids$from_id, ids$to_id)
}

# REMOVING ARC FUNCTION

#' Function for arc deletion in Dynamic Bayesian Network
#'
#' @param DBN object of class 'dbn'
#' @param from two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param to two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param cycle_OK boolean value, TRUE if cycle are admitted (default)
#'
#' @return object of class 'dbn'
#' @export
#'
#' @examples
#' DBN_example <- delete.arc.dbn(DBN=DBN_example,from=c('A','t_0'),to=c('R','t_0'))
#' DBN_example <- delete.arc.dbn(DBN=DBN_example,from=c('S','t'),to=c('O','t'))
#' DBN_example <- delete.arc.dbn(DBN=DBN_example,from=c('S','t-1'),to=c('S','t'))
delete.arc.dbn <- function(DBN, from, to) {
  if (!class(DBN) == 'dbn')
    stop("Error: DBN argument is not of class 'dbn'")
  if (!is.character(from))
    stop("Error: from is not a character")
  if (!is.character(to))
    stop("Error: to is not a character")
  from_ <- node_id(from[1], from[2])
  to_ <- node_id(to[1], to[2])
  if (to_ %in% DBN[['nodes']][[from[1]]][[from[2]]][['children']]) {
    DBN[['arcs']] <- delete_arc_element(DBN$arcs, c(from_, to_))
    DBN[['nodes']][[from[1]]][[from[2]]][['children']] <-
      delete_arc_element(DBN[['nodes']][[from[1]]][[from[2]]][['children']], to_)
    DBN[['nodes']][[to[1]]][[to[2]]][['parents']] <-
      delete_arc_element(DBN[['nodes']][[to[1]]][[to[2]]][['parents']], from_)
    if (!to_ %in% DBN[['nodes']][[from[1]]][[from[2]]][['parents']]) {
      DBN[['nodes']][[from[1]]][[from[2]]][['nbr']] <-
        delete_arc_element(DBN[['nodes']][[from[1]]][[from[2]]][['nbr']], to_)
      DBN[['nodes']][[to[1]]][[to[2]]][['nbr']] <-
        delete_arc_element(DBN[['nodes']][[to[1]]][[to[2]]][['nbr']], from_)
      DBN[['nodes']][[from[1]]][[from[2]]][['mb']] <-
        delete_arc_element(DBN[['nodes']][[from[1]]][[from[2]]][['mb']], to_)
      DBN[['nodes']][[to[1]]][[to[2]]][['mb']] <-
        delete_arc_element(DBN[['nodes']][[to[1]]][[to[2]]][['mb']], from_)
    }
    else if (identical(intersect(DBN[['nodes']][[from[1]]][[from[2]]][['children']], DBN[['nodes']][[to[1]]][[to[2]]][['children']]), character(0))) {
      DBN[['nodes']][[from[1]]][[from[2]]][['mb']] <-
        delete_arc_element(DBN[['nodes']][[from[1]]][[from[2]]][['mb']], to_)
      DBN[['nodes']][[to[1]]][[to[2]]][['mb']] <-
        delete_arc_element(DBN[['nodes']][[to[1]]][[to[2]]][['mb']], from_)
    }
    if (!identical(DBN[['nodes']][[to[1]]][[to[2]]][['parents']], character(0))) {
      for (i in DBN[['nodes']][[to[1]]][[to[2]]][['parents']]) {
        pi <- parse_node_id(i)
        if (identical(intersect(DBN[['nodes']][[from[1]]][[from[2]]][['children']], DBN[['nodes']][[pi$name]][[pi$time]][['children']]), character(0))) {
          DBN[['nodes']][[from[1]]][[from[2]]][['mb']] <-
            delete_arc_element(DBN[['nodes']][[from[1]]][[from[2]]][['mb']], i)
          DBN[['nodes']][[pi$name]][[pi$time]][['mb']] <-
            delete_arc_element(DBN[['nodes']][[pi$name]][[pi$time]][['mb']], from_)
        }
      }
    }
  }
  else{
    stop("ERROR: Arc do not exists!!!")
  }
  DBN
}

# REVERSING ARC FUNCTION

#' Function for arc reversal in Dynamic Bayesian Network
#'
#' @param DBN object of class 'dbn'
#' @param from two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param to two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param cycle_OK boolean value, TRUE if cycle are admitted (default)
#'
#' @return object of class 'dbn'
#' @export
#'
#' @examples
#' DBN_example <- add.arc.dbn(DBN=DBN_example,from=c('A','t-1'),to=c('R','t'))
reverse.arc.dbn <- function(DBN, from, to, cycle_OK = TRUE) {
  if (!class(DBN) == 'dbn')
    stop("Error: DBN argument is not of class 'dbn'")
  if (!is.character(from))
    stop("Error: from is not a character")
  if (!is.character(to))
    stop("Error: to is not a character")
  from_ <- node_id(from[1], from[2])
  to_ <- node_id(to[1], to[2])
  if (!(to_ %in% DBN[['nodes']][[from[1]]][[from[2]]][['children']])) {
    stop("ERROR: An arc that does not exist could not be reversed!!!")
  }
  if (any(apply(DBN$learning$blacklist, 1, function(x)
    (x[['from']] == to_ & x[['to']] == from_)))) {
    stop("ERROR: Temporal arcs are irreversible in DBNs!!!")
  }
  DBN <-
    DynamicBayesianNetwork::delete.arc.dbn(DBN = DBN, from = from, to = to)
  DBN <-
    DynamicBayesianNetwork::add.arc.dbn(
      DBN = DBN,
      from = to,
      to = from,
      cycle_OK = TRUE
    )
  DBN
}
