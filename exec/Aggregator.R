#!/usr/bin/Rscript
################################################################################
# Aggregator.R
################################################################################
# 2019-02-23
# Curtis Miller
################################################################################
# Take simulation data and convert it to a plottable, tabular format.
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
check_envir_has_objects <- CPAT:::check_envir_has_objects
stop_with_message <- CPAT:::stop_with_message

################################################################################
# MAIN FUNCTION DEFINITION
################################################################################

main <- function(input, output = NULL, TESTINPUT, alpha = 0.05, help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  test_tools <- new.env()
  input_objects <- new.env()
  load(TESTINPUT, envir = test_tools)
  load(input, envir = input_objects)
  test_tools_expected_objects <- c("pval_functions")
  input_objects_expected_objects <- c("res_list")
  check_envir_has_objects(test_tools_expected_objects, envir = test_tools,
                          blame_string = TESTINPUT)
  check_envir_has_objects(input_objects_expected_objects, envir = input_objects,
                          blame_string = input)

  pval_functions <- test_tools$pval_functions
  # TODO: curtis: BELOW: ADD CHECK FOR NAMES OF pval_functions VECTOR -- Sat 23 Feb 2019 10:53:55 PM MST
  stop_with_message(is.vector(pval_functions) &
                    all(sapply(pval_functions,
                               function(f) {"stat" %in% names(formals(f))})),
                    "Invalid pval_functions from" %s% TESTINPUT %s0% ";" %s%
                    "must be a vector of functions that all take input 'stat'")

  # TODO: curtis: CREATE AGGREGATED DATA FRAME -- Sat 23 Feb 2019 10:52:22 PM MST
}

################################################################################
# INTERFACE SETUP
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Take simulation data and convert to tabular format.",
        option_list = list(
          make_option(c("--input", "-i"), type = "character",
                      help = "Input file with simulation data"),
          make_option(c("--output", "-o"), type = "character", default = NULL,
                      help = "Output .Rda file containing aggregated" %s%
                             "simulation data"),
          make_option(c("--alpha", "-a"), type = "double", default = 0.05,
                      help = "Statistical significance threshold"),
          make_option(c("--TESTINPUT", "-T"), type = "character",
                      help = "Name of .Rda file defining how to compute" %s%
                             "p-values for included statistical tests")
        )
      ))

  do.call(main, cl_args)
}

