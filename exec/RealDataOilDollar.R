#!/usr/bin/Rscript
################################################################################
# RealDataOilDollar.R
################################################################################
# 2019-07-19
# Curtis Miller
################################################################################
# A file containing data for data example involving the price of oil and the
# strength of the U.S. dollar.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(output = "RealDataOilDollar.Rda") {
  # This function will be executed when the script is called from the command
  # line

  suppressPackageStartupMessages(library(CPAT))

  data("OilDollar")

  data_set <- OilDollar[49:nrow(OilDollar),]

  events <- data.frame("Time" = as.Date(character()),
                       "Event" = character(),
                       stringsAsFactors = FALSE)
  model <- Oil ~ DollarExchange
  is_ts <- TRUE

  save(data_set, events, model, is_ts, file = output, ascii = TRUE)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Creates data for example of price of oil vs. U.S. dollar")
  p <- add_argument(p, "output", type = "character",
                    default = "RealDataOilDollar.Rda",
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

