#!/usr/bin/Rscript
################################################################################
# SimTestCUSUM.R
################################################################################
# 2019-02-28
# Curtis Miller
################################################################################
# Define CUSUM test statistic for simulations.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# MAIN FUNCTION DEFINITION
################################################################################

main <- function(output = "SimTestCUSUM.Rda", linetype = "longdash",
                 help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  cusum_st <- function(formula, data) {
    res <- residuals(CPAT:::wrapped_dynlm(formula = formula, data = data))
    CPAT:::stat_Vn(res, use_kernel_var = TRUE, kernel = "qs", bandwidth = "and")
  }

  ##############################################################################
  # REQUIRED OBJECTS
  ##############################################################################

  stat_functions <- c("VnResid" = cusum_st)
  pval_functions <- c("VnResid" = function(q, d) {1 - CPAT:::pkolmogorov(q)})
  plot_desc <- c("VnResid" = linetype)

  save(stat_functions, pval_functions, plot_desc, file = output, ascii = TRUE)
}

################################################################################
# INTERFACE SETUP
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = paste("Define CUSUM test statistic for simulations."),
        option_list = list(
          make_option(c("--output", "-o"), type = "character",
                      default = "SimTestCUSUM.Rda",
                      help = "Name of output .Rda file"),
          make_option(c("--linetype", "-l"), type = "character",
                      default = "longdash",
                      help = "When plotted, type of line to plot")
        )
      ))

  do.call(main, cl_args)
}

