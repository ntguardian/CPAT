#!/usr/bin/Rscript
################################################################################
# MakeData.R
################################################################################
# 2018-09-22
# Curtis Miller
################################################################################
# Make the data sets used in project.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(ff_file, b_file, help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  ff <- read.csv(ff_file, row.names = 1)
  banks <- read.csv(b_file, row.names = 1)

  dput(ff)
  dput(banks)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Create data files used in project",
        option_list = list(
          make_option(c("--ff-file", "-f"), type = "character",
              help = "The location of the CSV file containing Fama-French data"),
          make_option(c("--bank-file", "-b"), type = "character",
              help = "The location of the CSV file containing stock data")        )
      ))

  names(cl_args) <- c("ff_file", "b_file", "help")
  do.call(main, cl_args)
}

