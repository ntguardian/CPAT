#!/usr/bin/Rscript
################################################################################
# SimDataNormalXYHetero.R
################################################################################
# 2019-02-22
# Curtis Miller
################################################################################
# Data simulation definition file; simulate data from Normal distribution, and
# residuals from i.i.d. Normal distribution, with a change in variance a fifth
# of the way through the sample.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# MAIN FUNCTION DEFINITION
################################################################################

main <- function(output = "SimDataNormalXY.Rda", help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  ##############################################################################
  # REQUIRED OBJECTS
  ##############################################################################

  eps_generator <- function(n) {
    n1 <- ceiling(n/2)
    n2 <- n - n1
    c(rnorm(n1), rnorm(n2, sd = sqrt(10)))
  }

  df_generator <- function(n, beta, eps) {
    d <- length(beta)
    const <- rep(1, times = n)
    if (d > 1) {
      interim_mat <- matrix(rnorm(n * (d - 1), mean = 1), ncol = d - 1)
      interim_mat <- cbind(const, interim_mat)
    } else {
      interim_mat <- as.matrix(const)
    }

    y <- as.vector(interim_mat %*% as.matrix(beta)) + eps

    if (d == 1) {
      return(data.frame("y" = y))
    }
    interim_mat <- interim_mat[, 2:d, drop = FALSE]
    colnames(interim_mat) <- paste0("X", 1:(d - 1))

    as.data.frame(cbind("y" = y, interim_mat))
  }

  save(eps_generator, df_generator, file = output, ascii = TRUE)
}

################################################################################
# INTERFACE SETUP
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Defines data generation functions for simulations",
        option_list = list(
          make_option(c("--output", "-o"), type = "character",
                      default = "SimDataNormalXYHetero.Rda",
                      help = "Name of output .Rda file")
        )
      ))

  do.call(main, cl_args)
}

