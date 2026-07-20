# ---- number of free parameters of a fitted node (for AIC/BIC df) ------------
dbn_node_df <- function(node) {
  if (is.dnode(node)) {
    d <- dim(node$prob)
    (d[1] - 1) * prod(d[-1])            # (levels - 1) * #parent-configurations
  } else if (is.cgnode(node)) {
    # one regression (coef rows) + one variance per discrete-parent combination
    ncol(node$coefficients) * (nrow(node$coefficients) + 1)
  } else if (is.gnode(node)) {
    length(node$regs) + 1               # intercept + slopes + residual variance
  } else {
    stop("logLik.dbn.fit: unrecognized node type (no prob/regs/coefficients)")
  }
}


# ---- per-node log-density on a reshaped frame -------------------------------
# Each returns a numeric vector with one entry per row of `df`; NA where the
# observation cannot be scored (unseen level / parent configuration, or an NA
# parameter left by replace.unidentifiable = FALSE).

# discrete: index the CPT array by (node value, parent values). dimnames(prob)
# is (node, parents...) - the same node ids used as column names in the reshaped
# frame. Observed labels are turned into integer positions (NA for a level not
# seen in training) so unscorable rows yield NA instead of an indexing error.
loglik_dnode <- function(node, df) {
  prob <- node$prob
  dn   <- dimnames(prob)
  dims <- names(dn)
  idx <- vapply(seq_along(dims),
                function(k) match(as.character(df[[dims[k]]]), dn[[k]]),
                integer(nrow(df)))
  if (!is.matrix(idx)) idx <- matrix(idx, nrow = nrow(df))
  ok <- stats::complete.cases(idx)
  out <- rep(NA_real_, nrow(df))
  if (any(ok)) out[ok] <- log(prob[idx[ok, , drop = FALSE]])
  out
}

# gaussian: mean = intercept + sum(slope * parent value); dnorm on the residual.
loglik_gnode <- function(node, df) {
  regs <- node$regs
  parents <- setdiff(names(regs), intercept_name)
  mean <- rep(unname(regs[[intercept_name]]), nrow(df))
  for (p in parents) mean <- mean + unname(regs[[p]]) * as.numeric(df[[p]])
  stats::dnorm(as.numeric(df[[node$node]]), mean = mean,
               sd = unname(node$std), log = TRUE)
}

# conditional gaussian: pick the discrete-parent combination (columns of
# $coefficients / $sd are in expand.grid(dlevels) order), then dnorm.
loglik_cgnode <- function(node, df) {
  coef <- node$coefficients
  dlevels <- node$dlevels
  dpar <- names(dlevels)
  gpar <- setdiff(rownames(coef), intercept_name)

  combos <- expand.grid(dlevels, stringsAsFactors = FALSE)  # row i <-> column i
  combo_key <- do.call(paste, c(combos[dpar], sep = "\r"))
  row_key <- do.call(paste, c(lapply(dpar, function(p) as.character(df[[p]])),
                              list(sep = "\r")))
  combo <- match(row_key, combo_key)    # NA for an unseen discrete configuration

  mean <- coef[intercept_name, combo]
  for (g in gpar) mean <- mean + coef[g, combo] * as.numeric(df[[g]])
  stats::dnorm(as.numeric(df[[node$node]]), mean = mean,
               sd = node$sd[combo], log = TRUE)
}

eval_node_loglik <- function(node, df) {
  if (is.dnode(node)) loglik_dnode(node, df)
  else if (is.cgnode(node)) loglik_cgnode(node, df)
  else if (is.gnode(node)) loglik_gnode(node, df)
  else stop(paste("logLik.dbn.fit: unrecognized node type for", node$node))
}

# per-row log-likelihood contribution of a set of nodes over one frame
eval_frame_loglik <- function(object, node_names, df, na.rm) {
  if (length(node_names) == 0 || nrow(df) == 0) return(numeric(nrow(df)))
  mat <- vapply(node_names,
                function(nm) eval_node_loglik(object[[nm]], df),
                numeric(nrow(df)))
  if (!is.matrix(mat)) mat <- matrix(mat, nrow = nrow(df))
  rowSums(mat, na.rm = na.rm)
}


#' Log-likelihood of data under a fitted DBN
#'
#' Computes the plug-in (predictive) log-likelihood of \code{data} under a
#' fitted Dynamic Bayesian Network. 
#' Parameters are taken as-is from \code{object}, passing held-out data yields an
#'  out-of-sample predictive log-likelihood.
#'
#' A fitted DBN factorizes into an initial part (the \eqn{G_0} nodes, named
#' \code{var_0, var_1, ...}) and a transition part (the \code{var_t} nodes
#' conditioned on their lagged parents):
#' \deqn{\log L = \sum_{i \in G_0} \log P(x_i^{0..mo-1} \mid pa)
#'             + \sum_{t \ge mo} \sum_{i \in trans} \log P(x_i^{t} \mid x^{t-1..t-k}).}
#' The transition sum is the one-step-ahead predictive log-likelihood.
#'
#' @param object a fitted DBN of class \code{dbn.fit}.
#' @param data a long-format \code{data.frame} with columns \code{Sample_id},
#'   \code{Time} and one column per base variable (as accepted by
#'   \code{\link{dbn.fit}}).
#' @param nodes optional character vector restricting the sum to these nodes.
#'   Entries may be node ids (\code{"A_t"}, \code{"A_0"}) or base variable names
#'   (\code{"A"}, expanded to every slice of that variable).
#' @param component one of \code{"both"} (default, initial + transition),
#'   \code{"transition"} (one-step-ahead predictive log-likelihood only) or
#'   \code{"initial"} (initial-slice log-likelihood only).
#' @param by.sample logical; if \code{TRUE} return a numeric vector of
#'   per-trajectory log-likelihoods named by \code{Sample_id} instead of the
#'   total (useful for cross-validation).
#' @param na.rm logical; if \code{TRUE} drop observations that cannot be scored
#'   (unseen level / configuration, or \code{NA} parameters); otherwise the
#'   result is \code{NA} when any such observation is present.
#' @param ... ignored, for S3 compatibility.
#'
#' @return With \code{by.sample = FALSE} a numeric of class \code{"logLik"} with
#'   \code{nobs} and \code{df} attributes, so \code{\link[stats]{AIC}} and
#'   \code{\link[stats]{BIC}} work directly. With \code{by.sample = TRUE} a named
#'   numeric vector.
#'
#' @examples
#' \dontrun{
#' fit <- dbn.fit(dbn, data = train)
#' logLik(fit, train)                              # full joint, in-sample
#' logLik(fit, test, component = "transition")     # predictive, held-out
#' AIC(logLik(fit, train)); BIC(logLik(fit, train))
#' }
#'
#' @method logLik dbn.fit
#' @export
logLik.dbn.fit <- function(object, data, nodes,
                           component = c("both", "transition", "initial"),
                           by.sample = FALSE, na.rm = FALSE, ...) {
  if (!is.dbn.fit(object))
    stop("logLik.dbn.fit: object must be of class 'dbn.fit'")
  if (!is.data.frame(data))
    stop("logLik.dbn.fit: data must be a data.frame")
  component <- match.arg(component)

  g0_names <- grep("_[0-9]+$", names(object), value = TRUE)
  t_names  <- grep("_t$",      names(object), value = TRUE)

  # optional node filter: accept node ids or base variable names
  if (!missing(nodes) && !is.null(nodes)) {
    all_names <- names(object)
    resolved <- character(0)
    for (nd in nodes) {
      if (nd %in% all_names) {
        resolved <- c(resolved, nd)
      } else {
        matched <- all_names[vapply(all_names,
                                    function(x) get_variable_name(x) == nd,
                                    logical(1))]
        if (length(matched) == 0)
          stop(paste("logLik.dbn.fit: unrecognized node or variable:", nd))
        resolved <- c(resolved, matched)
      }
    }
    g0_names <- intersect(g0_names, resolved)
    t_names  <- intersect(t_names,  resolved)
  }

  mo <- max(1L, markov_order(object))   # >= 1 so build_shifted_df never gets 0
  need_initial <- component %in% c("both", "initial")
  need_trans   <- component %in% c("both", "transition")

  ids <- character(0)
  vals <- numeric(0)
  nobs <- 0L
  df_free <- 0

  if (need_initial && length(g0_names) > 0) {
    df0 <- build_g0_df(data, markov_order = mo)
    per_row <- eval_frame_loglik(object, g0_names, df0, na.rm)
    ids  <- c(ids, as.character(df0$Sample_id))
    vals <- c(vals, per_row)
    nobs <- nobs + nrow(df0)
    df_free <- df_free + sum(vapply(g0_names,
                                    function(nm) dbn_node_df(object[[nm]]),
                                    numeric(1)))
  }

  if (need_trans && length(t_names) > 0) {
    dfT <- build_shifted_df(data, markov_order = mo, separator = "-")
    per_row <- eval_frame_loglik(object, t_names, dfT, na.rm)
    ids  <- c(ids, as.character(dfT$Sample_id))
    vals <- c(vals, per_row)
    nobs <- nobs + nrow(dfT)
    df_free <- df_free + sum(vapply(t_names,
                                    function(nm) dbn_node_df(object[[nm]]),
                                    numeric(1)))
  }

  if (by.sample) {
    if (length(vals) == 0) return(numeric(0))
    out <- tapply(vals, ids, sum, na.rm = na.rm)
    stats::setNames(as.numeric(out), names(out))
  } else {
    total <- sum(vals, na.rm = na.rm)
    attr(total, "nobs") <- as.integer(nobs)
    attr(total, "df")   <- as.numeric(df_free)
    class(total) <- "logLik"
    total
  }
}
