#!/usr/bin/Rscript
################################################################################
# SimTestHS.R
################################################################################
# 2019-03-01
# Curtis Miller
################################################################################
# Define Hidalgo-Seo test statistic for simulations.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(output, linetype = "dotdash") {
  # This function will be executed when the script is called from the command
  # line

  hs_st <- function(formula, data) {
    CPAT:::stat_hs_reg(formula = formula, data = data)
  }

  ################################################################################
  # REQUIRED OBJECTS
  ################################################################################

  stat_functions <- c("HS" = hs_st)
  pval_functions <- c("HS" = function(q, d) {1 - CPAT:::phidalgo_seo(q)})
  plot_desc <- c("HS" = linetype)

  save(stat_functions, pval_functions, plot_desc, file = output, ascii = TRUE)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Define Hidalgo-Seo test statistic for simulations.")
  p <- add_argument(p, "--output", type = "character", nargs = 1,
                    default = "SimTestHS.R",
                    help = "Name of output .Rda file")
  p <- add_argument(p, "--linetype", type = "character", default = "dotdash",
                    nargs = 1,
                    help = "When plotted, type of line to plot")

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

