#!/usr/bin/Rscript
################################################################################
# RemakePackage.R
################################################################################
# 2018-09-10
# Curtis Miller
################################################################################
# Remake the CPAT package
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

`%s%` <- function(x, y) {paste(x, y)}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  pack_name <- "CPAT"
  while (grepl(pack_name, getwd()) & basename(getwd()) != "" &
         basename(getwd()) != pack_name & !grepl(":", basename(getwd()))) {
    setwd("../")
  }

  if (basename(getwd()) != pack_name) stop("I can't go back to the base" %s%
                                           "directory of the package" %s% 
                                           pack_name)
  devtools::document()
  devtools::build_vignettes()
  devtools::install()
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Remake the package with devtools"
      ))

  do.call(main, cl_args)
}

