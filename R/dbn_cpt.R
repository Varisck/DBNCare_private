# DEFINE CPTs FROM TERMINAL

#' Function for Conditionally Probability Table definition by terminal
#'
#' @param var_name target variable to be conditioned (optional - default: ' ')
#' @param DBN object of class 'dbn' (optional - default: an empty DBN with no nodes)
#' @param nodes_level named list with the listed levels for each node
#'
#' @return multi-dimensional numerical vector
#' @export
#'
#' @examples
#' CPT <- define_CPT() #no information about the variables
#' B_0.prob <- define_CPT(var_name='B_0', DBN=DBN_example) #definition of the CPT for target variable B_0 given the DBN
define_CPT <-
  function(var_name = ' ',
           DBN = empty.dbn(dynamic_nodes = c(), markov_order = 1),
           nodes_level = list()) {
    if (!class(DBN) == 'dbn')
      stop("ERROR: DBN argument is not of class 'dbn'")
    if (!is.character(var_name))
      stop("ERROR: var_name is not a character")
    if (var_name == ' ' & length(DBN$nodes) == 0) {
      var_name <- readline(prompt = "target variable name: ")
      n_vars_conditioning_set <-
        readline(prompt = "number of variables conditioning set: ")
      conditioning_set <- c()
      if (n_vars_conditioning_set < 0) {
        stop("ERROR: number of variables in the conditioning set cannot be negative")
      }
      else if (n_vars_conditioning_set > 0) {
        for (j in 1:max(1, n_vars_conditioning_set)) {
          conditioning <-
            readline(prompt = paste("conditional variable", j, ': '))
          if (conditioning %in% c(var_name, conditioning_set)) {
            stop("ERROR: variable name already used")
          }
          conditioning_set <- c(conditioning_set, conditioning)
        }
      }
    }
    else{
      var_star <-
        var_star <- strsplit(var_name, '_')[[1]]
      var_1 <-
        paste(var_star[1:(length(var_star) - 1)], collapse = '_')
      temp_var <-
        ifelse(var_star[length(var_star)] == '0', 't_0', var_star[length(var_star)])
      if (class(DBN) == 'dbn' &
          var_1 %in% names(DBN$nodes) &
          temp_var %in% names(DBN$nodes[[var_1]])) {
        cat('target variable name:', var_name)
        n_vars_conditioning_set <-
          length(DBN$nodes[[var_1]][[temp_var]][['parents']])
        cat('\nconditioning set:', DBN$nodes[[var_1]][[temp_var]][['parents']])
        conditioning_set <- c()
        for (v in DBN$nodes[[var_1]][[temp_var]][['parents']]) {
          conditioning_set <- c(conditioning_set, v)
        }
      }
      else{
        stop("ERROR: var_name must be a DBN node")
      }
    }
    if (var_name %in% names(nodes_level)) {
      var_n_levels <- length(nodes_level[[var_name]])
      variable_set <- list()
      variable_set[[var_name]] <- nodes_level[[var_name]]
      cat('\n', var_name, "levels:", nodes_level[[var_name]], '\n')
    }
    else{
      var_n_levels <-
        as.numeric(readline(prompt = paste(var_name, "number of levels: ")))
      if (var_n_levels <= 1) {
        stop("ERROR: target variable must have at least 2 levels")
      }
      var_levels <- c()
      for (i in 1:var_n_levels) {
        var_levels <-
          c(var_levels, readline(prompt = paste("level", i, 'of', var_name, ': ')))
      }
      variable_set <- list()
      variable_set[[var_name]] <- var_levels
    }
    conditioning_set_levels <- c(var_n_levels)
    if (n_vars_conditioning_set > 0) {
      for (j in 1:n_vars_conditioning_set) {
        conditioning <- conditioning_set[j]
        if (conditioning %in% names(nodes_level)) {
          conditioning_n_levels <- length(nodes_level[[conditioning]])
          conditioning_set_levels <-
            c(conditioning_set_levels, conditioning_n_levels)
          variable_set[[conditioning]] <-
            nodes_level[[conditioning]]
          cat(conditioning, "levels:", nodes_level[[conditioning]], '\n')
        }
        else{
          conditioning_n_levels <-
            as.numeric(readline(prompt = paste(conditioning, "number of levels: ")))
          if (conditioning_n_levels <= 1) {
            stop("ERROR: variables in the conditioning set must have at least 2 levels")
          }
          conditioning_set_levels <-
            c(conditioning_set_levels, conditioning_n_levels)
          conditioning_levels <- c()
          for (k in 1:conditioning_n_levels) {
            conditioning_levels <-
              c(conditioning_levels, readline(prompt = paste(
                "level", k, 'of', conditioning, ': '
              )))
          }
          variable_set[[conditioning]] <- conditioning_levels
        }
      }
    }
    probabilities <- c()
    for (row in 1:nrow(expand.grid(variable_set))) {
      if (numbers::rem(row, var_n_levels) == 1) {
        prob_sum = 0
      }
      for (var in 1:length(variable_set)) {
        if (var == 1) {
          string <-
            ifelse(
              n_vars_conditioning_set > 0,
              paste("P(", var_name, "=", expand.grid(variable_set)[row, 1], '|'),
              paste("P(", var_name, "=", expand.grid(variable_set)[row, 1], ') =')
            )
        }
        else if (var == length(variable_set)) {
          string <-
            paste(
              string,
              conditioning_set[var - 1],
              '=',
              expand.grid(variable_set)[row, var],
              ifelse(((prob_sum == 1) | (numbers::rem(row, var_n_levels) == 0)
              ), ')', ') =')
            )
        }
        else{
          string <-
            paste(string,
                  conditioning_set[var - 1],
                  '=',
                  expand.grid(variable_set)[row, var],
                  ',')
        }
      }
      if (prob_sum == 1) {
        cat(paste(string, 'automatically set to 0\n'))
        probabilities = c(probabilities, 0)
        next
      }
      else if (numbers::rem(row, var_n_levels) == 0) {
        cat(paste(
          string,
          'automatically set to',
          as.character(1 - prob_sum),
          '\n'
        ))
        probabilities = c(probabilities, 1 - prob_sum)
        next
      }
      prob = as.numeric(readline(prompt = string))
      if (prob < 0 | prob > 1) {
        stop("ERROR: Probabilities must vary in [0,1]")
      }
      prob_sum = prob_sum + prob
      if (prob_sum > 1) {
        stop("ERROR: The sum of probabilities given the conditioning set exceed 1")
      }
      probabilities <- c(probabilities, prob)
    }
    cat('\n')
    assign(
      gsub(' ', '', paste(var_name, '.prob')),
      array(probabilities, dim = conditioning_set_levels, dimnames = variable_set)
    )
    get(gsub(' ', '', paste(var_name, '.prob')))
  }

#' Function for Conditionally Probability Tables definition by terminal given the Dynamic Bayesian Network
#'
#' @param DBN object of class 'dbn'
#'
#' @return list of multi-dimensional vector (CPTs for each DBN node)
#' @export
#'
#' @examples
#' define_CPTs(DBN_example)
define_CPTs <-
  function(DBN = empty.dbn(dynamic_nodes = c(), markov_order = 1)) {
    static_nodes <- bnlearn::node.ordering(from_DBN_to_G_0(DBN))
    dynamic_nodes <-
      quanteda::char_select(bnlearn::node.ordering(from_DBN_to_G_transition(DBN)),
                            "*t",
                            valuetype = "glob")
    CPTs <- list()
    defined_levels <- list()
    for (i in static_nodes) {
      CPT <-
        define_CPT(var_name = i,
                   DBN = DBN,
                   nodes_level = defined_levels)
      CPTs[[i]] <- CPT
      def_lev <-
        c(dimnames(CPT),
          setNames(dimnames(CPT), array(unlist(
            sapply(names(dimnames(CPT)) , function(x) {
              gsub(" ", "", paste(substring(x, 1, (nchar(
                x
              ) - 2)), '_t')) # generalize for var_name length > 1
            })
          ))),
          setNames(dimnames(CPT), array(unlist(
            sapply(names(dimnames(CPT)) , function(x) {
              gsub(" ", "", paste(substring(x, 1, (nchar(
                x
              ) - 2)), '_t-1')) # generalize for var_name length > 1
            })
          ))))
      defined_levels <-
        c(defined_levels, def_lev[setdiff(names(def_lev), names(defined_levels))])
    }
    for (j in dynamic_nodes) {
      CPT <-
        define_CPT(var_name = j,
                   DBN = DBN,
                   nodes_level = defined_levels)
      CPTs[[j]] <- CPT
    }
    CPTs
  }
