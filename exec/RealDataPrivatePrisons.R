#!/usr/bin/Rscript
################################################################################
# RealDataPrivatePrisons.R
################################################################################
# 2019-06-12
# Curtis Miller
################################################################################
# A file containing data for data example involving the behavior of a portfolio
# of private prison stock around the time the DOJ declared discontinuation of
# private prisons in 2016.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(output = "RealDataPrivatePrisons.Rda") {
  # This function will be executed when the script is called from the command
  # line

  suppressPackageStartupMessages(library(CPAT))
  suppressPackageStartupMessages(library(xts))

  data("ff")
  data("CXW")
  data("GEO")

  GEO <- GEO[index(GEO) <= as.Date("2017-10-31"), ]

  ff <- as.zoo(ff, order.by = as.Date(rownames(ff), format = "%Y%m%d"))

  data_set <- merge(ff, CXW, GEO)
  data_set$P.RF <- 0.6 * data_set$CXW + 0.4 * data_set$GEO - data_set$RF

  events <- data.frame("Time" = as.Date("2016-08-18"),
                       "Event" = "DOJ Announcement",
                       stringsAsFactors = FALSE)
  model <- P.RF ~ Mkt.RF + SMB + HML + RMW + CMA
  is_ts <- TRUE

  save(data_set, events, model, is_ts, file = output, ascii = TRUE)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Generate file for data example involving the behavior" %s%
                  "of private prison stock around the time the DOJ declared" %s%
                  "discontinuation of private prisons in 2016.")
  p <- add_argument(p, "output", type = "character",
                    default = "RealDataPrivatePrisons.Rda",
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

