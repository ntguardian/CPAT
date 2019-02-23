#!/usr/bin/Rscript
################################################################################
# PowerSimRegression.R
################################################################################
# 2019-02-19
# Curtis Miller
################################################################################
# Computes simulated change point test statistics in the OLS linear regression
# context depending on how data should be simulated and what change point
# contexts are requested, saving results in a .Rda file.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# FUNCTIONS
################################################################################

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`
stop_with_message <- CPAT:::stop_with_message
check_envir_has_objects <- CPAT:::check_envir_has_objects
base_file_name <- CPAT:::base_file_name

#' Print Progress Report
#'
#' Creates a progress report string for reporting progress.
#'
#' @param n Sample size
#' @param cpt Change point string
#' @param stat Statistic string
#' @param r Regime string
#' @param n_vector Vector of all sample sizes
#' @param cpt_vector Named vector of all change point functions
#' @param stat_vector Named vector of all statistic functions
#' @param r_list Named list of all regimes
#' @return A string serving as a progress report
#' @examples
#' progress_report(50, "sqrt", "mean", "mean2", c(50, 100),
#'                 cpt_vector = c("sqrt" = sqrt),
#'                 stat_vector = c("mean" = mean, "sd" = sd),
#'                 r_list = list("mean1" = cbind(1, 1), "mean2" = cbind(1, 2)))
progress_report <- function(n, cpt, stat, r, n_vector, cpt_vector, stat_vector,
                            r_list) {
  counts <- c(length(n_vector), length(cpt_vector), length(stat_vector),
              length(r_list))
  idx <- c(which(n == n_vector), which(cpt == names(cpt_vector)),
           which(stat == names(stat_vector)), which(r == names(r_list)))
  perc <- round(100 * (idx - 1) / counts)
  names(counts) <- c("n", "cpt", "stat", "r")
  names(idx) <- names(counts)
  names(perc) <- names(counts)

  total_work_needed <- prod(counts)
  total_work_done <- prod(idx)
  total_perc <- round(100 * (total_work_done - 1) / total_work_needed)

  sprintf("Total: %%%3d | n: %%%3d | k*: %%%3d | Stat.: %%%3d | Regime: %%%3d",
          total_perc, perc["n"], perc["cpt"], perc["stat"], perc["r"])
}

################################################################################
# MAIN FUNCTION DEFINITION
################################################################################

main <- function(SIMINPUT, CONTEXTINPUT, TESTINPUT, output = NULL,
                 seed = 20190219, seedless = FALSE, replications = 5000,
                 cores = NULL, verbose = FALSE, notsafe = FALSE, help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  library(doParallel)
  library(doRNG)
  library(purrr)

  if (is.null(cores)) {
    cores = min(1, detectCores() - 1)
  }
  registerDoParallel(cores = cores)
  if (!seedless) {
    registerDoRNG(seed)
  }

  simulation_tools <- new.env()
  context_tools <- new.env()
  test_tools <- new.env()
  load(SIMINPUT, envir = simulation_tools)
  load(CONTEXTINPUT, envir = context_tools)
  load(TESTINPUT, envir = test_tools)
  simulation_tools_expected_objects <- c("df_generator", "eps_generator")
  context_tools_expected_objects <- c("n_values", "kstar_functions",
                                      "struc_models")
  test_tools_expected_objects <- c("stat_functions")
  check_envir_has_objects(simulation_tools_expected_objects,
                          envir = simulation_tools, blame_string = SIMINPUT)
  check_envir_has_objects(context_tools_expected_objects, envir = context_tools,
                          blame_string = CONTEXTINPUT)
  check_envir_has_objects(test_tools_expected_objects, envir = test_tools,
                          blame_string = TESTINPUT)

  # It should be safe now to pull what we want from the files we loaded
  # But always check assumptions
  n_values <- context_tools$n_values
  stop_with_message(is.integer(n_values) & all(n_values > 0),
                    "Invalid n_values from" %s% CONTEXTINPUT %s0% "; must" %s%
                    "be positive integers")
  kstar_functions <- context_tools$kstar_functions
  stop_with_message(is.vector(kstar_functions) &
                    is.function(kstar_functions[1] &
                      !any(is.null(names(kstar_functions)))),
                    "Invalid kstar_functions from" %s% CONTEXTINPUT %s0% ";" %s%
                    "must be named vector of functions")
  struc_models <- context_tools$struc_models
  stop_with_message(is.list(struc_models) &
                    !any(is.null(names(struc_models))) &
                    all(sapply(struc_models, is.matrix)) &
                    all(sapply(struc_models, function(m) {ncol(m) == 2})) &
                    is.numeric(as.vector(sapply(struc_models, as.vector))),
                    "Invalid struc_models from" %s% CONTEXTINPUT %s0% "; it" %s%
                    "must be a named list of numeric matrices that have two" %s%
                    "columns")

  stat_functions <- test_tools$stat_functions
  stop_with_message(is.vector(stat_functions) & is.function(stat_function[1]) &
                    !any(is.null(names(stat_functions))), "Invalid" %s%
                    "stat_functions from" %s% TESTINPUT %s0% "; must be" %s%
                    "a named vector of functions")
  
  df_generator <- simulation_tools$df_generator
  stop_with_message(is.function(df_generator) & all(c("n", "beta", "eps") %in%
                    names(formals(df_generator))), "Invalid df_generator" %s%
                    "from" %s% SIMINPUT %s0% "; must be a function with" %s%
                    "arguments n, beta, and eps")
  eps_generator <- simulation_tools$eps_generator
  stop_with_message(is.function(df_generator) & "n" %in%
                    names(formals(eps_generator)), "Invalid eps_generator" %s%
                    "from" %s% SIMINPUT %s0% "; must be a function with" %s%
                    "argument n")

  # Create list for storing output
  res_list <- sapply("n" %s0% n_values, USE.NAMES = TRUE, simplify = FALSE,
                     FUN = function(x1) {
                       sapply(names(kstar_functions), USE.NAMES = TRUE,
                         simplify = FALSE, FUN = function(x2) {
                           sapply(names(stat_functions),
                             USE.NAMES = TRUE, simplify = FALSE,
                             FUN = function(x3) {
                               sapply(names(struc_models),
                                 USE.NAMES = TRUE, simplify = FALSE,
                                 FUN = function(y) {
                                   rep(NA, times = replications)
                                 }
                               )
                             })
                         })
                     })

  if (is.null(output)) {
    output <- paste("sims", base_file_name(CONTEXTINPUT),
                    base_file_name(SIMINPUT), base_file_name(TESTINPUT),
                    sep = "_") %s0% ifelse(seedless,
                                           "", "_" %s0% seed) %s0% ".Rda"
  }

  ##############################################################################
  # MAIN LOOP
  ##############################################################################

  for (n in n_values) {
    nname <- "n" %s0% n
    for (kstar_name in names(kstar_functions)) {
      kstar <- kstar_functions[kstar_name]
      k <- floor(kstar(n))
      stop_with_message(is.numeric(k) & k >= 1 & k <= n, "The function in" %s%
                        "kstar_functions called" %s% kstar_name %s% "did" %s%
                        "not return a number between 1 and" %s% n %s% "when" %s%
                        "expected; this should happen when" %s% kstar_name %s0%
                        "(n) is called")
      for (stat_name in names(stat_functions)) {
        stat <- stat_functions[stat_name]
        stop_with_message(all(c("data", "formula") %in% names(formals(stat))),
                          "The function in stat_functions called" %s% 
                          stat_name %s% "does not have all of the" %s%
                          "expected arguments (data, formula)")
        for (r_name in names(struc_models)) {
          r <- struc_models[r_name]
          beta1 <- r[, 1]
          beta2 <- r[, 2]
          d <- nrow(r)

          simulation <- foreach(throwaway = 1:replications,
                                .combine = c) %dorng% {
            eps <- eps_generator(n = n)
            df <- rbind(df_generator(n = k, beta = beta1, eps = eps[1:k]),
                        df_generator(n = (n - k), beta = beta2,
                                     eps = eps[max((k + 1), n):n]))
            stop_with_message("y" %in% names(df), "y must be one of the" %s%
                              "columns generated by the data frame returned" %s%
                              "by df_generator(); it is the dependent variable")
            names(df[which(names(df) != "y")]) <- paste0("X", 1:d)
            f <- "y" %s% "~" %s% paste0("X", 1:d, collapse = " + ")
            f <- as.formula(f)
            
            stat(formula = f, data = df)
          }
          res_list[[nname]][[kstar_name]][[stat_name]][[r_name]] <- simulation
          if (!notsafe) {
            save(res_list, file = output)
          }

          if (verbose) {
            progress_string <- progress_report(n, kstar_name, stat_name, r_name,
                                               n_values, kstar_functions,
                                               stat_functions, struc_models)
            cat(progress_string, "\r")
          }
        }
      }
    }
  }

  cat("Output File:", output %s0% "\n")
  save(res_list, file = output)
}

################################################################################
# INTERFACE SETUP
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Computes simulated change point test statistics in" %s%
                      "OLS linear regression context.",
        option_list = list(
          make_option(c("--SIMINPUT", "-S"), type = "character",
                      help = "Name of .Rda file containing data defining" %s%
                             "how data simulations are to be performed"),
          make_option(c("--CONTEXTINPUT", "-C"), type = "character",
                      help = "Name of .Rda file defining what contexts" %s%
                             "in which to perform simulations"),
          make_option(c("--TESTINPUT", "-T"), type = "character",
                      help = "Name of .Rda file defining statistical tests" %s%
                             "that are simulated"),
          make_option(c("--output", "-o"), type = "character", default = NULL,
                      help = "Output .Rda file for saving simulated test" %s%
                             "statistics; if not specified, name" %s%
                             "automatically generated"),
          make_option(c("--seed", "-s"), type = "integer", default = 20190219,
                      help = "The seed for simulations"),
          make_option(c("--seedless", "-R"), action = "store_true",
                      default = FALSE,
                      help = "Don't set a seed (causes --seed to be ignored)"),
          make_option(c("--replications", "-N"), type = "integer",
                      default = 5000,
                      help = "Number of replications per context"),
          make_option(c("--cores", "-c"), type = "integer", default = NULL,
                      help = "Number of cores to use in simulation; set to" %s%
                             "1 to effectively disable parallelization, or" %s%
                             "leave unset to use default (all but one core)"),
          make_option(c("--verbose", "-v"), action = "store_true",
                      default = FALSE, help = "Report progress on screen"),
          make_option(c("--notsafe", "-a"), action = "store_true",
                      default = FALSE,
                      help = "Don't save frequently; only save after all" %s%
                             "simulations are done")
        )
      ))

  do.call(main, cl_args)
}

