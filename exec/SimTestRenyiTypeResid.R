#!/usr/bin/Rscript
################################################################################
# SimTestRenyiTypeResid.R
################################################################################
# 2019-02-22
# Curtis Miller
################################################################################
# Define Rényi-type test statistic for simulations that uses the residuals of
# the regression model for the test.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# MAIN FUNCTION DEFINITION
################################################################################

main <- function(output = "SimTestRenyiTypeResid.Rda", linetype = "dashed",
                 help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  renyi_st <- function(formula, data) {
    res <- residuals(lm(formula = formula, data = data))
    CPAT:::stat_Zn(res, use_kernel_var = TRUE, kernel = "pa", bandwidth = "and",
                   kn = sqrt)
  }

  ##############################################################################
  # REQUIRED OBJECTS
  ##############################################################################

  stat_functions <- c("ZnResid" = renyi_st)
  pval_functions <- c("ZnResid" = function(q) {1 - CPAT:::pZn(q)})
  plot_desc <- c("ZnResid" = linetype)

  save(stat_functions, pval_functions, plot_desc, file = output, ascii = TRUE)
}

################################################################################
# INTERFACE SETUP
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = paste("Define Rényi-type test statistic for simulations",
                            "that uses the residuals of the regression model",
                            "for the test."),
        option_list = list(
          make_option(c("--output", "-o"), type = "character",
                      default = "SimTestRenyiTypeResid.Rda",
                      help = "Name of output .Rda file"),
          make_option(c("--linetype", "-l"), type = "character",
                      default = "dashed",
                      help = "When plotted, type of line to plot")
        )
      ))

  do.call(main, cl_args)
}

