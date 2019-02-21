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

################################################################################
# MAIN FUNCTION DEFINITION
################################################################################

main <- function(SIMINPUT, CONTEXTINPUT, output = NULL, seed = 20190219,
                 seedless = FALSE, replications = 5000, cores = NULL,
                 verbose = FALSE, help = FALSE) {
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
  load(SIMINPUT, envir = simulation_tools)
  load(CONTEXTINPUT, envir = context_tools)
  # TODO: curtis: CHECK THAT NEEDED OBJECTS ARE IN PLACE -- Tue 19 Feb 2019 05:53:30 PM MST
  simulation_tools_expected_objects <- c()
  context_tools_expected_objects <- c("n_values")
  check_envir_has_objects(simulation_tools_expected_objects,
                          envir = simulation_tools, blame_string = SIMINPUT)
  check_envir_has_objects(context_tools_expected_objects, envir = context_tools,
                          blame_string = CONTEXTINPUT)
  # It should be safe now to pull what we want from the files we loaded

  n_values <- context_tools$n_values
  stop_with_message(is.integer(n_values) & all(n_values > 0),
                    "Invalid n_values from" %s% CONTEXTINPUT %s0% "; must" %s%
                    "be positive integers")
  
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
                             "in which to perform simulations and how to" %s%
                             "compute test statistics"),
          make_option(c("--output", "-o"), type = "character", default = NULL,
                      help = "Output .Rda file for saving simulated test" %s%
                             "statistics; if not specified, name" %s%
                             "automatically generated")
          make_option(c("--seed", "-s"), type = "integer", default = 20190219,
                      help = "The seed for simulations"),
          make_option(c("--seedless", "-R"), action = "store_true",
                      default = FALSE,
                      help = "Don't set a seed (causes --seed to be ignored)")
          make_option(c("--replications", "-N"), type = "integer",
                      default = 5000,
                      help = "Number of replications per context"),
          make_option(c("--cores", "-c"), type = "integer", default = NULL,
                      help = "Number of cores to use in simulation; set to" %s%
                             "1 to effectively disable parallelization, or" %s%
                             "leave unset to use default (all but one core)"),
          make_option(c("--verbose", "-v"), action = "store_true",
                      default = FALSE, help = "Report progress on screen")
        )
      ))

  do.call(main, cl_args)
}

