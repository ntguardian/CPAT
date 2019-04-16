#!/usr/bin/Rscript
################################################################################
# ExpandingWindowpValComputer.R 
################################################################################
# 2019-04-10
# Curtis Miller
################################################################################
# Given a data set of serial data, computes p-values for change point statistics
# on an expanding window of data and returns a data set containing those
# p-values.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`
check_envir_has_objects <- CPAT:::check_envir_has_objects
stop_with_message <- CPAT:::stop_with_message
is.formula <- function(x) {inherits(x, "formula")}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(input, statistics, output = "out.Rda", left = 1,
                 firstright = NA, lastright = NA, cores = NA) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  suppressPackageStartupMessages(library(doParallel))
  suppressPackageStartupMessages(library(doRNG))
  suppressPackageStartupMessages(library(formula.tools))
  suppressPackageStartupMessages(library(zoo))
  suppressPackageStartupMessages(library(purrr))

  if (is.na(cores)) {
    cores <- max(1, detectCores() - 1)
  }
  registerDoParallel(cores = cores)

  input_env <- new.env()
  load(input, envir = input_env)
  temp_env <- new.env()
  load(statistics[[1]], envir = temp_env)

  # Assumption checking
  input_env_expected_objects <- c("data_set", "events", "model", "is_ts")
  temp_env_expected_objects <- c("stat_functions", "pval_functions",
                                 "plot_desc")
  check_envir_has_objects(input_env_expected_objects, envir = input_env,
                          blame_string = input)
  check_envir_has_objects(temp_env_expected_objects, envir = temp_env,
                          blame_string = statistics[[1]])

  # Load objects and check assumptions
  model <- input_env$model
  stop_with_message(is.formula(model), "Invalid model from" %s% input %s0%
                    "; must be a formula")

  data_set <- input_env$data_set
  stop_with_message((is.data.frame(data_set) |
                      (is.zoo(data_set) & length(dim(data_set)) == 2)) &
                    all(all.vars(model) %in% names(data_set)),
                    "Invalid data_set from" %s% input %s0% "; must be a" %s%
                    "data frame containing the variables listed in model")

  events <- input_env$events
  stop_with_message(is.data.frame(events) &
                    all(c("idx", "event") %in% names(events)) &
                    (is.character(events[["event"]]) | nrow(events) == 0),
                    "Invalid events from" %s% input %s0% "; must be a" %s%
                    "data frame with columns named idx and event and event" %s%
                    "is character data")

  is_ts <- input_env$is_ts
  stop_with_message(is.logical(is_ts) & length(is_ts) == 1,
                    "Invalid is_ts from" %s% input %s0% "; must be a single" %s%
                    "boolean value")

  stat_functions <- temp_env$stat_functions
  stop_with_message(is.vector(stat_functions) &
                    is.function(stat_functions[[1]]) &
                    !any(is.null(names(stat_functions))), "Invalid" %s%
                    "stat_functions from" %s% statistics[[1]] %s0% "; must" %s%
                    "be a named vector of functions")

  plot_desc <- temp_env$plot_desc
  stop_with_message(is.vector(plot_desc) & is.character(plot_desc) &
                    all(!is.null(names(plot_desc))),
                    "Invalid plot_desc from" %s% statistics[[1]] %s0% ";" %s%
                    "must be a named character vector")

  pval_functions <- temp_env$pval_functions
  stop_with_message(is.vector(pval_functions) &
                    all(sapply(pval_functions,
                               function(f) {all(c("q", "d") %in%
                                 names(formals(f)))})) &
                    all(!is.null(names(pval_functions))),
                    "Invalid pval_functions from" %s% statistics[[1]] %s0%
                    "; must be a vector of functions that all take input" %s0%
                    "'q' and 'd'")

  for (f in statistics[2:length(statistics)]) {
    if (is.na(f)) break
    load(f, envir = temp_env)
    # Import and check for errors
    check_envir_has_objects(temp_env_expected_objects, envir = temp_env,
                            blame_string = f)
    temp_stat <- temp_env$stat_functions
    stop_with_message(is.vector(temp_stat) &
                      is.function(temp_stat[[1]]) &
                      !any(is.null(names(temp_stat))), "Invalid" %s%
                      "stat_functions from" %s% f %s0% "; must be a named" %s%
                      "vector of functions")

    temp_desc <- temp_env$plot_desc
    stop_with_message(is.vector(temp_desc) & is.character(temp_desc) &
                      all(!is.null(names(temp_desc))),
                      "Invalid plot_desc from" %s% f %s0% "; must be a" %s%
                      "named character vector")

    temp_pval <- temp_env$pval_functions
    stop_with_message(is.vector(temp_pval) &
                      all(sapply(temp_pval,
                                 function(g) {all(c("q", "d") %in%
                                   names(formals(g)))})) &
                      all(!is.null(names(temp_pval))),
                      "Invalid pval_functions from" %s% f %s0% "; must be" %s%
                      "a vector of functions that all take input 'q' and" %s0%
                      "'d'")

    # Merge statistics, desc, and objects
    plot_desc <- c(plot_desc, temp_desc[which(
        !(names(temp_desc) %in% names(plot_desc)))])
    stat_functions <- c(stat_functions, temp_stat[which(
        !(names(temp_stat) %in% names(stat_functions)))])
    pval_functions <- c(pval_functions, temp_pval[which(
        !(names(temp_pval) %in% names(pval_functions)))])
  }

  # Initialization
  n <- nrow(data_set)
  if (is.na(lastright)) {
    lastright <- n
  }
  if (is.na(firstright)) {
    firstright <- min(n, left + length(get.vars(model)) + 1)
  }
  stop_with_message(firstright >= left & lastright >= left &
                    firstright <= lastright, "Must have left <= firstright" %s%
                                             "<= lastright")

  if (is_ts) {
    idx <- time(data_set)[firstright:lastright]
    starttime <- time(data_set)[left]
  } else {
    idx <- firstright:lastright
    starttime <- left
  }

  temp_reg <- CPAT:::wrapped_dynlm(formula = model,
                                   data = data_set[left:firstright, ])
  d <- ncol(model.matrix(temp_reg))

  stat_pvals <- lapply(names(stat_functions), function(s) {
    stat <- stat_functions[[s]]
    pval <- partial(pval_functions[[s]], d = d)

    foreach(t = firstright:lastright, .combine = c) %dorng% {
      stat_val <- stat(formula = model, data = data_set[left:t, ])
      pval(stat_val)
    }
  })

  names(stat_pvals) <- names(stat_functions)
  stat_pvals <- data.frame(stat_pvals)
  stat_pvals <- cbind(stat_pvals, "End" = idx)
  save(stat_pvals, plot_desc, events, file = output, ascii = TRUE)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Given a data set of serial data, computes p-values for" %s%
                  "change point statistics on an expanding window of data" %s%
                  "and returns a data set containing those p-values")
  p <- add_argument(p, "input", type = "character",
                    help = "Input file containing data set for analysis")
  p <- add_argument(p, "--output", type = "character", default = "out.Rda",
                    help = "Output file containing p-values and events")
  p <- add_argument(p, "--statistics", type = "character", 
                    nargs = Inf,
                    help = "Files containing tests to apply to data")
  p <- add_argument(p, "--left", type = "integer", default = 1,
                    help = "The left anchor of all windows of data")
  p <- add_argument(p, "--firstright", type = "integer", default = NA,
                    help = "The first right anchor of windows; by default," %s%
                           "the input to left plus the dimension of the model")
  p <- add_argument(p, "--lastright", type = "integer", default = NA,
                    help = "The last right anchor of windows; by default," %s%
                           "the last point in the data set")
  p <- add_argument(p, "--cores", type = "integer", default = NA,
                    help = "The number of cores for parallel computation;" %s%
                           "by default, all but one, or one if only one")

  cl_args <- parse_args(p)
  cl_args <- cl_args[!(names(cl_args) %in% c("opts", "help"))]

  do.call(main, cl_args[2:length(cl_args)])
}

