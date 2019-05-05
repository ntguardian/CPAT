#!/usr/bin/Rscript
################################################################################
# RealDataCXWDonaldTrump.R
################################################################################
# 2019-05-04
# Curtis Miller
################################################################################
# A file containing data for data example involving the behavior of CXW around
# the election of Donald Trump.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

`%s%` <- CPAT:::`%s%`

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(output = "RealDataCXWDonaldTrump.Rda") {
  # This function will be executed when the script is called from the command
  # line

  suppressPackageStartupMessages(library(CPAT))

  data("ff")
  data("CXW")

  ff <- as.zoo(ff, order.by = as.Date(rownames(ff), format = "%Y%m%d"))

  data_set <- merge(ff, CXW)
  events <- data.frame("Time" = as.Date("2016-11-08"),
                       "Event" = "U.S. Election",
                       stringsAsFactors = FALSE)
  model <- I(CXW - RF) ~ Mkt.RF + SMB + HML + RMW + CMA
  is_ts <- TRUE

  save(data_set, events, model, is_ts, file = output, ascii = TRUE)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Generate file for data example involving the behavior" %s%
                  "of stock CXW around the election of Donald Trump")
  p <- add_argument(p, "output", type = "character",
                    default = "RealDataCXWDonaldTrump.Rda",
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

