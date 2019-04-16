#!/usr/bin/Rscript
################################################################################
# RealDataM1Germany.R
################################################################################
# 2019-04-15
# Curtis Miller
################################################################################
# A file containing data for data example involving the M1 supply of Germany.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(output = "RealDataM1Germany.Rda") {
  # This function will be executed when the script is called from the command
  # line

  suppressPackageStartupMessages(library(dynlm))

  data("M1Germany")

  data_set <- M1Germany
  events <- data.frame("idx" = 122, "event" = "Monetary union",
                       stringsAsFactors = FALSE)
  model <- d(logm1) ~ d(L(loggnp, 2)) + d(interest) + d(L(interest)) +
                      d(logprice) + L(logm1) + L(loggnp) + L(interest) +
                      season(logm1, ref = 4)
  is_ts <- TRUE

  save(data_set, events, model, is_ts, file = output, ascii = TRUE)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("This is a template for executable R scripts.")
  p <- add_argument(p, "output", type = "character",
                    default = "RealDataM1Germany.Rda",
                    help = "Name of output .Rda file")

  cl_args <- parse_args(p)
  cl_args <- cl_args[!(names(cl_args) %in% c("help", "opts"))]
  if (any(sapply(cl_args, is.na))) {
    # User did not specify all inputs; print help message
    print(p)
    cat("\n\nNot all needed inputs were given.\n")
    quit()
  }

  do.call(main, cl_args[2:length(cl_args)])
}

