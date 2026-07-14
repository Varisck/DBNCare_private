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
  if (is.dbn(dbn))
    return(from_DBN_to_G_transition(dbn))
  if (is.dbn.fit(dbn))
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
  if (is.dbn(dbn))
    return(from_DBN_to_G_0(dbn))
  if (is.dbn.fit(dbn))
    return(from_fitted_DBN_to_fitted_G_0(dbn))
}



# --------------- Internal methods --------------- 

#' This takes a fit object and renames the class names from the dbn to bn,
#' both for each node and for the network itself.
#'
#' Nodes are relabeled from 'dbn.fit.{cg,d,g}node' to 'bn.fit.{cg,d,g}node'.
#' The network class is set to c("bn.fit", "bn.fit.<net_type>") where net_type
#' is 'dnet' (all discrete nodes), 'gnet' (all Gaussian nodes) or 'cgnet'
#' (mixed / conditional Gaussian nodes).
#'
#' @param net a fitted network (list of nodes)
#'
#' @return the network with bn.fit classes
redefine_class = function(net) {
  pat <- "^dbn\\.fit\\.(cg|d|g)node$"
  node_types <- character(0)
  for (nm in names(net)) {
    cls <- class(net[[nm]])
    if (grepl(pat, cls)) {
      class(net[[nm]]) <- sub(pat, "bn.fit.\\1node", cls)
      node_types <- c(node_types, sub(pat, "\\1", cls))
    } else stop(paste("ERROR: class not recognized for node", net[[nm]]$node))
  }
  net_type <- if (all(node_types == "d")) "dnet"
              else if (all(node_types == "g")) "gnet"
              else "cgnet"
  class(net) <- c("bn.fit", paste0("bn.fit.", net_type))
  net
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
  if (!is.dbn(DBN))
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
#' @details When the DBN was built with \code{extend_g0 = TRUE}, G_0 spans
#'   several initial slices stored under keys \code{t_0}, \code{t_1}, ...,
#'   \code{t_(markov_order - 1)} (dynamic nodes only); static nodes always have
#'   only a \code{t_0} slice. Each slice \code{t_i} becomes a node \code{var_i}.
#'
#' @examples
#' G_0 <- from_DBN_to_G_0(DBN_example)
from_DBN_to_G_0 <- function(DBN) {
  if (!is.dbn(DBN))
    stop("ERROR: DBN argument is not of class 'dbn'")
  PN <- empty_bn_shell()
  for (j in names(DBN$nodes)) {
    # G_0 may span several initial slices (t_0, t_1, ...) when extend_g0 = TRUE;
    # t-i / t / type keys are excluded by the anchored pattern.
    g0_slices <- grep("^t_[0-9]+$", names(DBN$nodes[[j]]), value = TRUE)
    for (slice_key in g0_slices) {
      child_id <- paste0(j, "_", sub("^t_", "", slice_key))
      PN[['nodes']][[child_id]] <- DBN[['nodes']][[j]][[slice_key]]
      for (k in DBN[['nodes']][[j]][[slice_key]][['parents']]) {
        PN[['arcs']] <- rbind(PN$arcs, c(k, child_id))
      }
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
  if (!is.dbn.fit(DBN_fitted))
    stop("ERROR: DBN_fitted argument is not of class 'dbn.fit'")
  # G_0 nodes end in _<digits> (var_0, var_1, ...); this captures every initial
  # slice when extend_g0 = TRUE, unlike the glob "*0" which only matched var_0.
  # Transition nodes (var_t) and lagged parent refs (var_t-i) are not matched.
  g0_names <- grep("_[0-9]+$", names(DBN_fitted), value = TRUE)
  BN_0_fitted <- DBN_fitted[g0_names]
  BN_0_fitted = redefine_class(BN_0_fitted)
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
  if (!is.dbn.fit(DBN_fitted))
    stop("ERROR: DBN_fitted argument is not of class 'dbn.fit'")
  BN_transition_fitted <-
    DBN_fitted[grep("_t$", names(DBN_fitted), value = TRUE)]
  BN_transition_fitted = redefine_class(BN_transition_fitted)
  BN_transition_fitted
}
