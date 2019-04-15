#!/usr/bin/Rscript
################################################################################
# SimTestRenyiTypeRegNoKern.R
################################################################################
# 2019-03-07
# Curtis Miller
################################################################################
# Define the Rényi-type test statistic for simulations that uses the full
# regression model coefficients for the test and does not use kernel methods for
# long-run variance estimation.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# MAIN FUNCTION DEFINITION
################################################################################

main <- function(output = "SimTestRenyiTypeRegNoKern.Rda", linetype = "solid",
                 help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  renyi_reg_st <- function(formula, data) {
    CPAT:::stat_Zn_reg(formula = formula, data = data, kn = sqrt,
                       use_kernel_var = FALSE)
    # CPAT:::stat_Zn_reg_r(formula = formula, data = data, kn = sqrt)
  }

  ##############################################################################
  # REQUIRED OBJECTS
  ##############################################################################

  stat_functions <- c("ZnRegNoKern" = renyi_reg_st)
  pval_functions <- c("ZnRegNoKern" = function(q, d) {1 - CPAT:::pZn(q, d = d)})
  plot_desc <- c("ZnRegNoKern" = linetype)

  save(stat_functions, pval_functions, plot_desc, file = output, ascii = TRUE)
}

################################################################################
# INTERFACE SETUP
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = paste("Define Rényi-type test statistic for simulations",
                            "that uses the regression model coefficients for",
                            "the test."),
        option_list = list(
          make_option(c("--output", "-o"), type = "character",
                      default = "SimTestRenyiTypeRegNoKern.Rda",
                      help = "Name of output .Rda file"),
          make_option(c("--linetype", "-l"), type = "character",
                      default = "solid",
                      help = "When plotted, type of line to plot")
        )
      ))

  do.call(main, cl_args)
}

