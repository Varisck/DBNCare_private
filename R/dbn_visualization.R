# Visualisation of Dynamic Bayesian Networks.
#
# `plot()` methods for `dbn` and `dbn.fit` objects. Nodes are laid out on a grid
# of time slices: the current slice (`X_t` / `X_0`) forms the bottom row and each
# lag (`X_t-1`, `X_t-2`, ... or `X_1`, `X_2`, ...) is stacked above it, up to the
# Markov order of the network. Within a row the nodes are spread evenly across a
# fixed width, so every row spans the same horizontal extent regardless of how
# many nodes it holds.


# Nodes and arcs of the network to be plotted. Returns a list with `node_ids`
# (character vector) and `edges` (a two-column from/to character matrix).
.dbn_plot_elements <- function(x, network) {
  # Build from -> to pairs from each node's stored parents.
  parent_edges <- function(net, ids) {
    e <- do.call(rbind, lapply(ids, function(n) {
      ps <- net[[n]]$parents
      if (length(ps)) cbind(ps, n) else NULL
    }))
    if (is.null(e)) matrix(character(0), ncol = 2) else unname(e)
  }

  if (is.dbn(x)) {
    g <- if (network == "g_t") get.transition.net(x) else get.g0.net(x)
    node_ids <- names(g$nodes)
    edges <- if (nrow(g$arcs)) unname(g$arcs) else matrix(character(0), ncol = 2)
  } else if (is.dbn.fit(x)) {
    if (network == "g_t") {
      # A fitted DBN stores only the X_t nodes; the lagged parents X_t-k live as
      # parent references, so add them back as nodes.
      t_ids <- grep("_t$", names(x), value = TRUE)
      edges <- parent_edges(x, t_ids)
      node_ids <- unique(c(t_ids, edges[, 1]))
    } else {
      # Prior-network nodes are the initial slices X_0, X_1, ... (var_<digits>).
      node_ids <- grep("_[0-9]+$", names(x), value = TRUE)
      edges <- parent_edges(x, node_ids)
    }
  } else {
    stop("plot: x must be an object of class 'dbn' or 'dbn.fit'")
  }

  list(node_ids = node_ids, edges = edges)
}


# Render the time-slice grid layout with visNetwork.
.dbn_visnetwork <- function(node_ids, edges, size) {
  variables <- vapply(node_ids, get_variable_name, character(1))
  lags      <- vapply(node_ids, get_variable_time_index, numeric(1))
  max_lag   <- max(lags)

  # x: within each lag row the nodes are spread evenly across [0, size], ordered
  # by variable name so the same variable lines up in a column across rows.
  x <- numeric(length(node_ids))
  for (l in unique(lags)) {
    idx <- which(lags == l)
    idx <- idx[order(variables[idx])]
    n   <- length(idx)
    x[idx] <- if (n == 1) size / 2 else seq(0, size, length.out = n)
  }

  # y: lag 0 (X_t / X_0) forms the bottom row. visNetwork's y axis points down,
  # so the bottom row takes the largest y. Rows are spaced like the horizontal
  # node gap to keep the grid visually even.
  per_row <- max(table(lags))
  v_gap   <- if (per_row > 1) size / (per_row - 1) else size
  y <- (max_lag - lags) * v_gap

  nodes_df <- data.frame(
    id = node_ids, label = node_ids, x = x, y = y,
    physics = FALSE, stringsAsFactors = FALSE
  )
  edges_df <- data.frame(
    from = edges[, 1], to = edges[, 2],
    stringsAsFactors = FALSE
  )

  visNetwork::visNetwork(nodes_df, edges_df) %>%
    visNetwork::visNodes(shape = "ellipse", font = list(size = size / 45)) %>%
    visNetwork::visEdges(arrows = "to") %>%
    visNetwork::visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE)
}


#' Plot a Dynamic Bayesian Network
#'
#' @description
#' Draws a DBN (`dbn`) or a fitted DBN (`dbn.fit`) as a time-slice grid using
#' \pkg{visNetwork}. The current slice sits on the bottom row and each lag is
#' stacked above it, up to the Markov order of the network:
#' \itemize{
#'   \item transition network (`g_t`): row 0 holds the `X_t` nodes, row 1 the
#'     `X_t-1` nodes, ..., up to `X_t-markov_order`;
#'   \item prior network (`g_0`): row 0 holds the `X_0` nodes, row 1 the `X_1`
#'     nodes, ..., up to `X_(markov_order - 1)` when the DBN was built with
#'     `extend_g0 = TRUE`.
#' }
#' Within every row the nodes are spread evenly across the same horizontal width
#' (`size`), so a row of 5 nodes and a row of 10 nodes both span the full width;
#' shrink or grow `size` to make the graph more compact or more spread out.
#'
#' Because `plot.dbn` / `plot.dbn.fit` are S3 methods for the base
#' \code{\link[graphics]{plot}} generic, `plot(dbn)` works out of the box.
#'
#' @param x an object of class `dbn` or `dbn.fit`.
#' @param network which network to draw: `"g_t"` (transition network, the
#'   default) or `"g_0"` (prior network).
#' @param size horizontal width each row is spread across. Defaults to
#'   `num_nodes * 1000`, where `num_nodes` is the number of variables in the
#'   network.
#' @param ... further arguments (currently unused).
#'
#' @return a \pkg{visNetwork} \code{htmlwidget}, drawn when printed.
#' @name plot.dbn
#'
#' @examples
#' \dontrun{
#' dbn <- random.structure.dbn(c("A", "B", "C"), .5, .5, markov_order = 2)
#' plot(dbn)                       # transition network
#' plot(dbn, network = "g_0")      # prior network
#' plot(dbn, size = 3000)          # more compact
#' }
NULL

#' @rdname plot.dbn
#' @method plot dbn
#' @export
plot.dbn <- function(x, network = c("g_t", "g_0"), size = NULL, ...) {
  if (!is.dbn(x)) stop("plot.dbn: x must be an object of class 'dbn'")
  network <- match.arg(network)
  el <- .dbn_plot_elements(x, network)
  if (is.null(size))
    size <- length(unique(vapply(el$node_ids, get_variable_name, character(1)))) * 1000
  .dbn_visnetwork(el$node_ids, el$edges, size)
}

#' @rdname plot.dbn
#' @method plot dbn.fit
#' @export
plot.dbn.fit <- function(x, network = c("g_t", "g_0"), size = NULL, ...) {
  if (!is.dbn.fit(x)) stop("plot.dbn.fit: x must be an object of class 'dbn.fit'")
  network <- match.arg(network)
  el <- .dbn_plot_elements(x, network)
  if (is.null(size))
    size <- length(unique(vapply(el$node_ids, get_variable_name, character(1)))) * 1000
  .dbn_visnetwork(el$node_ids, el$edges, size)
}
