#!/usr/bin/Rscript
################################################################################
# RealDataGDPPMIPredict.R
################################################################################
# 2019-06-30
# Curtis Miller
################################################################################
# A file containing data for data example involving the prediction of YoY GDP
# using PMI.
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

main <- function(output = "RealDataGDPPMIPredict.Rda") {
  # This function will be executed when the script is called from the command
  # line

  suppressPackageStartupMessages(library(CPAT))

  data("GDPPMI")

  data_set <- GDPPMI[9:nrow(GDPPMI),]

  events <- data.frame("Time" = as.Date(character()),
                       "Event" = character(),
                       stringsAsFactors = FALSE)
  model <- GDP ~ L(PMI, 1) + L(PMI, 2) + L(PMI, 3) + L(PMI, 4)
  is_ts <- TRUE

  save(data_set, events, model, is_ts, file = output, ascii = TRUE)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Generate file for data example involving the prediction" %s%
                  "of YoY GDP using PMI")
  p <- add_argument(p, "output", type = "character",
                    default = "RealDataGDPPMIPredict.Rda",
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

