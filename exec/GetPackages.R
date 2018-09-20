#!/usr/bin/Rscript
################################################################################
# GetPackages.R
################################################################################
# 2018-09-19
# Curtis Miller
################################################################################
# This is a one-line description of the file.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(repos = "http://cran.us.r-project.org",
                 lib = NULL, destdir = NULL, help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  needed_packages <- c(
    "doParallel",
    "dplyr",
    "ggplot2",
    "gridExtra",
    "tikzDevice",
    "doParallel",
    "fGarch",
    "doRNG",
    "reshape2",
    "magrittr",
    "smoothmest",
    "foreach",
    "xts",
    "devtools",
    "Rcpp",
    "RcppArmadillo",
    "cointReg",
    "purrr",
    "Rdpack",
    "roxygen2",
    "testthat"
  )

  install.packages(needed_packages, lib = lib, destdir = destdir, repos = repos,
                   dependencies = TRUE)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Install needed packages for this project",
        option_list = list(
          make_option(c("--repos", "-r"), type = "character",
                      default = "http://cran.us.r-project.org",
                      help = "Repo from which to fetch packages"),
          make_option(c("--lib", "-l"), type = "character", default = NULL,
                      help = "Location of library to install packages"),
          make_option(c("--destdir", "-d"), type = "character", default = NULL,
                      help = "Directory where downloaded packages are stored")
        )
      ))

  do.call(main, cl_args)
}

