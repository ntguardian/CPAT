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

main <- function(output = "SimTestRenyiTypeResid.Rda", help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  renyi_st <- function(formula, data) {
    res <- residuals(lm(formula = formula, data = data))
    CPAT:::stat_Zn(res, use_kernel_var = TRUE)
  }

  ##############################################################################
  # REQUIRED OBJECTS
  ##############################################################################

  stat_functions <- c("ZnResid" = renyi_st)

  save(stat_functions, file = output)
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
                      help = "Name of output .Rda file")
        )
      ))

  do.call(main, cl_args)
}

