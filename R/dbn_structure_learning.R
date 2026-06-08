# check the data/blacklist and whitelist structure is correct
check_data_input <- function(data, nodes, blacklist, whitelist) {
  if (!all(c("Time", "Sample_id") %in% names(data)))
    stop("ERROR: Sample_id and Time are not in data columns")
  if (!is.character(data$Sample_id))
    stop("ERROR: Sample_id must be character")
  if (!is.numeric(data$Time) && !is.character(data$Time))
    stop("ERROR: Time must be numeric or character")
  if (any(is.na(data)))
    stop("ERROR: Missing data detected")
  if (!all(dplyr::count(data %>% dplyr::group_by(Sample_id, Time))$n == 1))
    stop("ERROR: One or more combinations of Sample_id and Time are repeated")
  if (any(is.na(data %>% tidyr::complete(Sample_id, Time))))
    stop("ERROR: One or more samples have incomplete temporal sequences")

  arc_pattern <- paste0("^(", paste(nodes, collapse = "|"), ")_(t|t-[0-9]+|0)$")
  to_zero_pat <- paste0("^(", paste(nodes, collapse = "|"), ")_0$")
  to_t_pat    <- paste0("^(", paste(nodes, collapse = "|"), ")_t$")
  from_t_pat  <- paste0("^(", paste(nodes, collapse = "|"), ")_(t|t-[0-9]+)$")

  validate_arc_set <- function(arcs, label) {
    if (is.null(arcs)) return(invisible())
    if (!"matrix" %in% class(arcs))
      stop(paste0("ERROR: ", label, " must be a matrix"))
    if (!all(colnames(arcs) == c("from", "to")))
      stop(paste0("ERROR: ", label, " must have columns 'from' and 'to'"))
    if (any(duplicated(as.data.frame(arcs))))
      stop(paste0("ERROR: ", label, " contains duplicates"))
    if (!all(grepl(arc_pattern, arcs)))
      stop(paste0("ERROR: ", label,
                  " entries must be of the form `var`_0, `var`_t or `var`_t-i,",
                  " where var is a column of data"))
    if (!all(grepl(to_zero_pat, arcs[, "to"]) ==
             grepl(to_zero_pat, arcs[, "from"])))
      stop(paste0("ERROR: ", label,
                  " arcs must be (var1_0, var2_0), (var1_t-i, var2_t)",
                  " or (var1_t, var2_t)"))
    if (!all(grepl(to_t_pat, arcs[, "to"]) ==
             grepl(from_t_pat, arcs[, "from"])))
      stop(paste0("ERROR: ", label,
                  " arcs into _t must originate from _t or _t-i"))
  }
  validate_arc_set(blacklist, "blacklist")
  validate_arc_set(whitelist, "whitelist")

  if (any(duplicated(as.data.frame(rbind(blacklist, whitelist)))))
    stop("ERROR: an arc cannot be both in the whitelist and in the blacklist")
}

# check the structure learning parameters 
check_sl_params <- function(method, test, score, max.sx) {
  is.naturalnumber <-
    function(x, tol = .Machine$double.eps^0.5) x > tol & abs(x - round(x)) < tol

  valid_tests  <- c(
    # discrete
    "mi", "mi-sh", "x2", "mc-mi", "smc-mi", "mi-adf", "x2-adf",
    "mc-x2", "smc-x2", "sp-mi", "sp-x2",
    # continuous (gaussian)
    "cor", "zf", "mi-g", "mi-g-sh", "mc-mi-g", "smc-mi-g",
    "mc-cor", "smc-cor", "mc-zf", "smc-zf",
    # mixed (conditional gaussian)
    "mi-cg")
  valid_scores <- c(
    # discrete
    "loglik", "aic", "bic", "ebic", "pred-loglik", "fnml", "qnml",
    "nal", "pnal", "bde", "bds", "bdj", "k2", "mbde", "bdla",
    # continuous (gaussian)
    "loglik-g", "aic-g", "bic-g", "ebic-g", "pred-loglik-g",
    "nal-g", "pnal-g", "bge",
    # mixed (conditional gaussian)
    "loglik-cg", "aic-cg", "bic-cg", "ebic-cg", "pred-loglik-cg",
    "nal-cg", "pnal-cg")

  if (!is.null(max.sx) && !is.naturalnumber(max.sx))
    stop("ERROR: max.sx must be a positive integer")

  if (!method %in% c("constraint", "score", "hybrid"))
    stop("ERROR: method must be 'score', 'constraint' or 'hybrid'")

  if (method %in% c("constraint", "hybrid") && !is.null(test) &&
      !test %in% valid_tests)
    stop(paste("ERROR: test must be one of:",
               paste(valid_tests, collapse = ", ")))

  if (method %in% c("score", "hybrid") && !is.null(score) &&
      !score %in% valid_scores)
    stop(paste("ERROR: score must be one of:",
               paste(valid_scores, collapse = ", ")))
}

# Map a detected dataset type (see dataset_type()) to bnlearn's matching
# default score / test, mirroring bnlearn's own per-type defaults:
#   discrete -> bic / mi ; gaussian -> bic-g / cor ; mixed -> bic-cg / mi-cg
dbn_default_score <- function(type) {
  switch(type, discrete = "bic", gaussian = "bic-g", mixed = "bic-cg",
         stop("ERROR: unrecognized dataset type '", type, "'"))
}
dbn_default_test <- function(type) {
  switch(type, discrete = "mi", gaussian = "cor", mixed = "mi-cg",
         stop("ERROR: unrecognized dataset type '", type, "'"))
}

# Resolve a user-supplied score / test. NULL or "auto" triggers detection of the
# data type (via dataset_type) and selection of the matching bnlearn default;
# any explicit value is passed through untouched.
resolve_score <- function(score, data) {
  if (is.null(score) || identical(score, "auto"))
    dbn_default_score(dataset_type(data)) else score
}
resolve_test <- function(test, data) {
  if (is.null(test) || identical(test, "auto"))
    dbn_default_test(dataset_type(data)) else test
}

#' Preprocessing function for Structure Learning algorithm.
#'
#' @param data object of type data.frame to be given as an input to structure learning algorithm
#' @param markov_order markov order to build the dataset
#' @param blacklist matrix with 'from' and 'to' columns of forbidden edges in the Dynamic Bayesian Network (default is NULL)
#' @param whitelist matrix with 'from' and 'to' columns of fixed edges in the Dynamic Bayesian Network (default is NULL)
#' @param allow_intraslice_edges = TRUE bool wether to allow edges of the form X_t -> Y_t
#' @param allow_t_0_edges = TRUE bool wether to allow edges of the form X_0 -> Y_0
#'
#' @return character vector made of ordered data to fit G_0, data to fit G_transition, blacklist to fit G_0, blacklist to fit G_transition, whitelist to fit G_0, whitelist to fit G_transition
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' dbn.prep.data(data)
dbn.prep.data = function(data, markov_order = 1, blacklist = NULL, whitelist = NULL, 
                          allow_intraslice_edges = TRUE,
                          allow_t_0_edges = TRUE) {

  exclude_cols = c("Sample_id", "Time")
  nodes = setdiff(names(data), exclude_cols)

  # checking input data is correct!
  check_data_input(data, nodes, blacklist, whitelist)

  # coerce character columns to factor (bnlearn requires factor for discrete
  # nodes); leave numeric columns untouched so gaussian / mixed networks work
  for (v in nodes) {
    if (is.character(data[[v]]) | is.logical(data[[v]])) data[[v]] <- as.factor(data[[v]])
  }

  # create dataframe of time 0
  df_0 = data[data$Time == 0,]
  names(df_0) = lapply(names(df_0), concat_name_post, postfix = "_0")
  # removes the cols in exclude_cols
  df_0 = select(df_0, -all_of(exclude_cols))
  
  # create dataframe with shifted cols
  df_transition = build_shifted_df(data, markov_order = markov_order, separator = "-")
  # remove Sample_id and Time form dataframe
  df_transition = select(df_transition, -all_of(exclude_cols))
  

  blacklist_new = blacklist_g0_gt(names(df_0),
                                  sapply(nodes, concat_name_post, postfix = "_t"),
                                  markov_order = markov_order,
                                  allow_intraslice_edges = allow_intraslice_edges,
                                  allow_t_0_edges = allow_t_0_edges)

  # only keep t|t-i -> t|t-i nodes and t_0 -> t_0
  blacklist_new = blacklist_new[which(arr.ind = T,
                                      (grepl("^.+(t|t-[0-9]+)$", 
                                            blacklist_new[, 'from'])
                                      & grepl("^.+(t|t-[0-9]+)$", blacklist_new[,'to'])
                                      ) |
                                        (
                                          grepl("^.+(_0)$", 
                                                blacklist_new[, 'from'])
                                          & grepl("^.+(_0)$", blacklist_new[,'to'])
                                        )
                                      ), , drop = FALSE]
    
  
  
  blacklist <- unique(rbind(blacklist, blacklist_new))
  
  # extract from blacklist both blacklist_0 and blacklist_t
  # get all pairs in blacklist where from has _0
  blacklist_0 = blacklist[which(arr.ind = T,
                                grepl("^.+(_0)$", blacklist[, "from"])), , drop = FALSE]
  # get all pairs in blacklist where from has _t(-i)
  blacklist_t = blacklist[which(arr.ind = T,
                                grepl("^.+(t|t-[0-9]+)$", blacklist[, "from"])), , drop = FALSE]

  # do the same for the whitelist
  whitelist_0 = whitelist[which(arr.ind = T,
                                grepl("^.+(_0)$", whitelist[, "from"])), , drop = FALSE]
  whitelist_t = whitelist[which(arr.ind = T,
                                grepl("^.+(t|t-[0-9]+)$", whitelist[, "from"])), , drop = FALSE]

  list(df_0, df_transition, blacklist_0, blacklist_t, whitelist_0, whitelist_t)
}

#' Generation of a Dynamic Bayesian Network given two compatible networks \eqn{G_0} and \eqn{G_{transition}}
#'
#' @param PN \eqn{G_0}, an object of class 'bn'
#' @param TN \eqn{G_{transition}}, an object of class 'bn'
#' @param markov_order markov order of the dbn to build
#'
#' @return An object of class 'dbn' with the characteristics of \eqn{G_0} and \eqn{G_{transition}}
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' c(data0, dataTR, ., blacklistTR, ., .) %<-% dbn.prep.data(data = data, markov_order = 1)
#' bn_sampled_0 <- pc.stable(data0)
#' bn_sampled_TR <- pc.stable(dataTR, blacklist = blacklistTR)
#' dbn.build.struct(PN = bn_sampled_0, TN =  bn_sampled_TR, markov_order = 1)
#' 
dbn.build.struct <- function(PN, TN, markov_order = 1) {
  if (!class(PN) == 'bn')
    stop("ERROR: PN argument is not of class 'bn'")
  if (!class(TN) == 'bn')
    stop("ERROR: TN argument is not of class 'bn'")
  
  # regenerate the full theoretical cross-slice forbidden set so that
  # DBN$learning$blacklist records every arc that is structurally forbidden
  # in this DBN class, not just what the algorithms happened to receive
  g_0_nodes = names(PN$nodes)
  g_t_nodes = names(TN$nodes)[grepl("_t$", names(TN$nodes))]
  full_blacklist = blacklist_g0_gt(g_0_nodes, g_t_nodes,
                                   markov_order = markov_order)

  DBN = list(
    learning = list(
      whitelist = rbind(PN$learning$whitelist, TN$learning$whitelist),
      blacklist = unique(rbind(full_blacklist,
                               PN$learning$blacklist,
                               TN$learning$blacklist)),
      test = list(
        G_0 = PN$learning$test, G_transition = TN$learning$test
        ),
      ntests = list(
        G_0 = PN$learning$ntests, G_transition = TN$learning$ntests
        ),
      algo = list(
        G_0 = PN$learning$algo, G_transition = TN$learning$algo
        ),
      args = list(
        G_0 = PN$learning$args, G_transition = TN$learning$args
        )
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

  variables = unlist(lapply(names(PN$nodes), get_variable_name))
  # creating DBN nodes list with the gt and g0
  for(variable in variables) {
    DBN$nodes[[variable]][['t_0']] = PN$nodes[[gsub(" ", "", 
                                                    paste(variable, "_0"))]]
    DBN$nodes[[variable]][['t']] = TN$nodes[[gsub(" ", "", 
                                                  paste(variable, "_t"))]]
    
    # copy nodes in transition net for all markov orders
    for(mo in 1:markov_order) {
      DBN$nodes[[variable]][[paste0('t-', mo)]] = TN$nodes[[gsub(" ", "",
                                                      paste(variable, "_t-", mo))]]
    }
    
    DBN$nodes[[variable]][['type']] = 'Dynamic'
  }
  DBN$arcs = rbind(PN$arcs, TN$arcs)
  class(DBN) = "dbn"
  DBN
}


#' Learn the equivalence class of DBN via PC-stable algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via PC stable algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' pc.stable.dbn(data)
#' 
pc.stable.dbn <- function(data, markov_order = 1, test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  test <- resolve_test(test, data)
  check_sl_params('constraint', test, NULL, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- pc.stable(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- pc.stable(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Grow-Shrink algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via Grow-Shrink algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' gs.dbn(data)
#' 
gs.dbn <- function(data, markov_order = 1, test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  test <- resolve_test(test, data)
  check_sl_params('constraint', test, NULL, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- gs(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- gs(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Incremental Association algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via Incremental Association algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' iamb.dbn(data)
#' 
iamb.dbn <- function(data, markov_order = 1, test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  test <- resolve_test(test, data)
  check_sl_params('constraint', test, NULL, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- iamb(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- iamb(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Fast Incremental Association algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via Fast Incremental Association algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' fast.iamb.dbn(data)
#' 
fast.iamb.dbn <- function(data, markov_order = 1, test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  test <- resolve_test(test, data)
  check_sl_params('constraint', test, NULL, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- fast.iamb(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- fast.iamb(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Interleaved Association algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via Interleaved Association algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' inter.iamb.dbn(data)
#' 
inter.iamb.dbn <- function(data, markov_order = 1, test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  test <- resolve_test(test, data)
  check_sl_params('constraint', test, NULL, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- inter.iamb(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- inter.iamb(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Incremental Association algorithm with FDR
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via Incremental Association algorithm with FDR
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' iamb.fdr.dbn(data)
#' 
iamb.fdr.dbn <- function(data, markov_order = 1, test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  test <- resolve_test(test, data)
  check_sl_params('constraint', test, NULL, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- iamb.fdr(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- iamb.fdr(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Max-Min Parents and Children algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via Max-Min Parents and Children algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' mmpc.dbn(data)
#' 
mmpc.dbn <- function(data, markov_order = 1, test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  test <- resolve_test(test, data)
  check_sl_params('constraint', test, NULL, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- mmpc(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- mmpc(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Semi-Interleaved HITON-PC algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via Semi-Interleaved HITON-PC algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' si.hiton.pc.dbn(data)
#' 
si.hiton.pc.dbn <- function(data, markov_order = 1, test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  test <- resolve_test(test, data)
  check_sl_params('constraint', test, NULL, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- si.hiton.pc(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- si.hiton.pc(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Hybrid Parents and Children algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via Hybrid Parents and Children algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' hpc.dbn(data)
#' 
hpc.dbn <- function(data, markov_order = 1, test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  test <- resolve_test(test, data)
  check_sl_params('constraint', test, NULL, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- hpc(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- hpc(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Hill Climbing algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param score a character string, the label of the scoring function to be used in the algorithm. Default `"auto"` selects the score from the data type via [dataset_type()]: `bic` (discrete), `bic-g` (gaussian) or `bic-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters of the scoring function
#'
#' @return An object of class 'dbn', which is a DBN learned via Hill Climbing algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' hc.dbn(data)
#' 
hc.dbn <- function(data, markov_order = 1, score = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  score <- resolve_score(score, data)
  check_sl_params('score', NULL, score, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- hc(data0, score = score, maxp = min(max.sx,length(names(data0))-1) , blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- hc(dataTR, score = score, maxp = min(max.sx,length(names(dataTR))-1) , blacklist = blacklistTR, whitelist = whitelistTR,...)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Tabu Search algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param score a character string, the label of the scoring function to be used in the algorithm. Default `"auto"` selects the score from the data type via [dataset_type()]: `bic` (discrete), `bic-g` (gaussian) or `bic-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters of the scoring function
#'
#' @return An object of class 'dbn', which is a DBN learned via Tabu Search algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' tabu.dbn(data)
#' 
tabu.dbn <- function(data, markov_order = 1, score = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  score <- resolve_score(score, data)
  check_sl_params('score', NULL, score, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- tabu(data0, score = score, maxp = min(max.sx,length(names(data0))-1) , blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- tabu(dataTR, score = score, maxp = min(max.sx,length(names(dataTR))-1) , blacklist = blacklistTR, whitelist = whitelistTR,...)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Hybrid HPC algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param score a character string, the label of the scoring function to be used in the algorithm. Default `"auto"` selects the score from the data type via [dataset_type()]: `bic` (discrete), `bic-g` (gaussian) or `bic-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters of the scoring function and/or chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via Hybrid HPC algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' h2pc.dbn(data)
#'
h2pc.dbn <- function(data, markov_order = 1, score = "auto", test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  score <- resolve_score(score, data)
  test  <- resolve_test(test, data)
  check_sl_params('hybrid', test, score, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- h2pc(data0, maximize.args = list(score = score, maxp = min(max.sx,length(names(data0))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(data0))-1), ...), blacklist = blacklist0, whitelist = whitelist0)
  bn_sampled_TR <- h2pc(dataTR, maximize.args = list(score = score, maxp = min(max.sx,length(names(dataTR))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(dataTR))-1), ...), blacklist = blacklistTR, whitelist = whitelistTR)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via Max-Min Hill Climbing algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param score a character string, the label of the scoring function to be used in the algorithm. Default `"auto"` selects the score from the data type via [dataset_type()]: `bic` (discrete), `bic-g` (gaussian) or `bic-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters of the scoring function and/or chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via  Max-Min Hill Climbing algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' mmhc.dbn(data)
#'
mmhc.dbn <- function(data, markov_order = 1, score = "auto", test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  score <- resolve_score(score, data)
  test  <- resolve_test(test, data)
  check_sl_params('hybrid', test, score, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- mmhc(data0, maximize.args = list(score = score, maxp = min(max.sx,length(names(data0))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(data0))-1), ...), blacklist = blacklist0, whitelist = whitelist0)
  bn_sampled_TR <- mmhc(dataTR, maximize.args = list(score = score, maxp = min(max.sx,length(names(dataTR))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(dataTR))-1), ...), blacklist = blacklistTR, whitelist = whitelistTR)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Learn the equivalence class of DBN via 2-phase Restricted Maximization algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param score a character string, the label of the scoring function to be used in the algorithm. Default `"auto"` selects the score from the data type via [dataset_type()]: `bic` (discrete), `bic-g` (gaussian) or `bic-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters of the scoring function and/or chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via 2-phase Restricted Maximization algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' rsmax2.dbn(data)
#'
rsmax2.dbn <- function(data, markov_order = 1, restrict = 'pc.stable', maximize = 'hc', score = "auto", test = "auto", max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  score <- resolve_score(score, data)
  test  <- resolve_test(test, data)
  check_sl_params('hybrid', test, score, max.sx)
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% dbn.prep.data(data = data, markov_order = markov_order, blacklist = blacklist, whitelist = whitelist)
  bn_sampled_0 <- rsmax2(data0, restrict = restrict, maximize = maximize, maximize.args = list(score = score, maxp = min(max.sx,length(names(data0))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(data0))-1), ...), blacklist = blacklist0, whitelist = whitelist0)
  bn_sampled_TR <- rsmax2(dataTR, restrict = restrict, maximize = maximize, maximize.args = list(score = score, maxp = min(max.sx,length(names(dataTR))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(dataTR))-1), ...), blacklist = blacklistTR, whitelist = whitelistTR)
  dbn.build.struct(PN = bn_sampled_0, TN = bn_sampled_TR, markov_order = markov_order)
}

#' Structure Learning function for Dynamic Bayesian Networks
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param markov_oder an integer refering to the markov order of the learning algorithm
#' @param algorithm a character string, the label of the structure learning algorithm to be used (default is `hc`)
#' @param algorithm.res a character string, the label of the structure learning algorithm to be used for the restriction phase of the hybrid algorithms (default is `pc.stable`)
#' @param algorithm.max a character string, the label of the structure learning algorithm to be used for the maximization phase of the hybrid algorithms (default is `hc`)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param test a character string, the label of the conditional independence test to be used in the algorithm. Default `"auto"` selects the test from the data type via [dataset_type()]: `mi` (discrete), `cor` (gaussian) or `mi-cg` (mixed)
#' @param score a character string, the label of the scoring function to be used in the algorithm. Default `"auto"` selects the score from the data type via [dataset_type()]: `bic` (discrete), `bic-g` (gaussian) or `bic-cg` (mixed)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param ... additional parameters of the scoring function and/or chosen test 
#'
#' @return An object of class 'dbn', which is a DBN learned via 2-phase Restricted Maximization algorithm
#' @export
#'
#' @examples
#' DBN <- random.structure.dbn(c("A","B","C"), 0.5, 0.5, markov_order = 1)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' dbn.learn.structure(data, algorithm='h2pc',score='aic') 
dbn.learn.structure <- function(data, markov_order = 1, algorithm = 'hc', algorithm.res = 'pc.stable', algorithm.max = 'hc', blacklist = NULL, whitelist = NULL, test = NULL, score = NULL, max.sx = NULL, ...){
  if (algorithm %in% c('hc','tabu')){
    get(paste(algorithm,'.dbn',sep=''))(data = data, markov_order = markov_order, score = score, blacklist = blacklist, whitelist = whitelist, max.sx = max.sx, ...)
  }
  else if (algorithm %in% c('pc.stable','gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr')){
    get(paste(algorithm,'.dbn',sep=''))(data = data, markov_order = markov_order, test = test, blacklist = blacklist, whitelist = whitelist, max.sx = max.sx, ...)
  }
  else if (algorithm %in% c('h2pc','mmhc') & algorithm.res %in% c('pc.stable','gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr') & algorithm.max %in% c('hc','tabu')){
    get(paste(algorithm,'.dbn',sep=''))(data = data, markov_order = markov_order, score = score, test = test, blacklist = blacklist, whitelist = whitelist, max.sx = max.sx, ...)
  }
  else if (algorithm == 'rsmax2' & algorithm.res %in% c('pc.stable','gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr') & algorithm.max %in% c('hc','tabu')){
    get(paste(algorithm,'.dbn',sep=''))(data = data, markov_order = markov_order, restrict = algorithm.res, maximize = algorithm.max, score = score, test = test, blacklist = blacklist, whitelist = whitelist, max.sx = max.sx, ...)
  }
  else{
    stop("ERROR: Algorithm must be one of the following: 'hc','tabu','pc.stable','gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr', 'h2pc','mmhc', 'rsmax2'")
  }
}
