#!/usr/bin/Rscript
################################################################################
# UnivariatePlotter.R
################################################################################
# 2019-02-25
# Curtis Miller
################################################################################
# Plot power simulation data when the changing variable is univariate, and save
# the plots.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

################################################################################
# FUNCTIONS
################################################################################

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(foo, bar = 0, baz = FALSE, help = FALSE, opts = NA) {
  # This function will be executed when the script is called from the command
  # line; the help and opts parameters does nothing, but are needed for
  # do.call() to work

  print("Hello!")
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

  if (sys.nframe() == 0) {
  p <- arg_parser("Plot power simulation data when the changing variable" %s%
                  "is univariate, and save the plots.")
  p <- add_argument(p, "foo", type = "character", nargs = 1,
                    help = "A positional command-line argument")
  p <- add_argument(p, "--bar", type = "integer", default = 0,
                    nargs = 1, help = "A command-line option")
  p <- add_argument(p, "--baz", flag = TRUE, help = "A command-line flag")

  cl_args <- parse_args(p)

  do.call(main, cl_args[2:length(cl_args)])
}

