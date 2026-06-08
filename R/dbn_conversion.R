# DBN TO G_TRANSITION (OF CLASS BN)


#' Transformation of a DBN in a G_transition network
#'
#' @param DBN object of class 'dbn' or "dbn.fit"
#'
#' @return object of class 'bn'
#' @export
#'
#' @examples
#' G_transition <- get.transition.net(DBN_example)
get.transition.net <- function(dbn) {
  if (class(dbn) == "dbn")
    return(from_DBN_to_G_transition(dbn))
  if (class(dbn) == "dbn.fit")
    return(from_fitted_DBN_to_fitted_G_transition(dbn))
}

#' Transformation of a DBN in a G_0 network
#'
#' @param DBN object of class 'dbn' or "dbn.fit"
#'
#' @return object of class 'bn'
#' @export
#'
#' @examples
#' G_0 <- from_DBN_to_G_0(DBN_example)
get.g0.net <- function(dbn) {
  if (class(dbn) == "dbn")
    return(from_DBN_to_G_0(dbn))
  if (class(dbn) == "dbn.fit")
    return(from_fitted_DBN_to_fitted_G_0(dbn))
}



#' Transformation of a DBN in a G_transition network
#'
#' @param DBN object of class 'dbn'
#'
#' @return object of class 'bn'
#'
#' @examples
#' G_transition <- from_DBN_to_G_transition(DBN_example)
from_DBN_to_G_transition <- function(DBN) {
  if (!class(DBN) == 'dbn')
    stop("ERROR: DBN argument is not of class 'dbn'")
  TN <- empty_bn_shell()
  for (i in 0:DBN$markov_order) {
    slice_key <- if (i == 0) 't' else paste0('t-', i)
    for (j in names(DBN$nodes)) {
      if (DBN[['nodes']][[j]][['type']] == 'Dynamic') {
        child_id <- node_id(j, slice_key)
        TN[['nodes']][[child_id]] <- DBN[['nodes']][[j]][[slice_key]]
        for (k in DBN[['nodes']][[j]][[slice_key]][['parents']]) {
          TN[['arcs']] <- rbind(TN$arcs, c(k, child_id))
        }
      }
    }
  }
  class(TN) <- "bn"
  TN
}

# DBN TO G_0 (OF CLASS BN)

#' Transformation of a DBN in a G_0 network
#'
#' @param DBN object of class 'dbn'
#'
#' @return object of class 'bn'
#'
#' @examples
#' G_0 <- from_DBN_to_G_0(DBN_example)
from_DBN_to_G_0 <- function(DBN) {
  if (!class(DBN) == 'dbn')
    stop("ERROR: DBN argument is not of class 'dbn'")
  PN <- empty_bn_shell()
  for (j in names(DBN$nodes)) {
    child_id <- node_id(j, 't_0')
    PN[['nodes']][[child_id]] <- DBN[['nodes']][[j]][['t_0']]
    for (k in DBN[['nodes']][[j]][['t_0']][['parents']]) {
      PN[['arcs']] <- rbind(PN$arcs, c(k, child_id))
    }
  }
  class(PN) <- "bn"
  PN
}



# DEFINE FUNCTION FROM dbn.fit TO G_0 bn.fit

#' Function for G_0 parameters set extraction
#'
#' @param DBN_fitted object of class 'dbn.fit'
#'
#' @return object of class 'bn.fit'
#'
#' @examples
#' fitted_0 <- from_fitted_DBN_to_fitted_G_0(fitted_DBN)
from_fitted_DBN_to_fitted_G_0 <- function(DBN_fitted) {
  if (!class(DBN_fitted) == 'dbn.fit')
    stop("ERROR: DBN_fitted argument is not of class 'dbn.fit'")
  BN_0_fitted <-
    DBN_fitted[quanteda::char_select(names(DBN_fitted), "*0", valuetype = "glob")]
  class(BN_0_fitted) <- "bn.fit"
  BN_0_fitted
}

# DEFINE FUNCTION FROM dbn.fit TO G_transition bn.fit

#' Function for G_transition parameters set extraction
#'
#' @param DBN_fitted object of class 'dbn.fit'
#'
#' @return object of class 'bn.fit'
#'
#' @examples
#' fitted_transition <- from_fitted_DBN_to_fitted_G_transition(fitted_DBN)
from_fitted_DBN_to_fitted_G_transition <- function(DBN_fitted) {
  if (!class(DBN_fitted) == 'dbn.fit')
    stop("ERROR: DBN_fitted argument is not of class 'dbn.fit'")
  BN_transition_fitted <-
    DBN_fitted[quanteda::char_select(names(DBN_fitted), "*t", valuetype = "glob")]
  class(BN_transition_fitted) <- "bn.fit"
  BN_transition_fitted
}
