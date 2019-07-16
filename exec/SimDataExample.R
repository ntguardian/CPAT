#!/usr/bin/Rscript
################################################################################
# SimDataExample.R
################################################################################
# 2019-03-01
# Curtis Miller
################################################################################
# Data simulation definition file; simulate data resembling that in example.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(output) {
  # This function will be executed when the script is called from the command
  # line

  library(zoo)
  library(dynlm)

  ##############################################################################
  # REQUIRED OBJECTS
  ##############################################################################

  eps_generator <- function(n) {
    n1 <- min(n, 62)
    n2 <- n - n1
    vec1 <- as.numeric(arima.sim(list(order = c(0, 0, 1), ma = c(0.287)),
                                 sd = sqrt(0.001), n = n1, n.start = 500))
    if (n2 == 0) {
      vec2 <- numeric()
    } else {
      vec2 <- as.numeric(rnorm(n2, sd = sqrt(0.001)))
    }
    c(vec1, vec2)
  }

  df_generator <- function(n, beta, eps) {
    d <- length(beta)
    stopifnot(d == 3)
    const <- rep(1, times = n)
    if (d > 1) {
      # Simulate RealEnergyPrice
      v1 <- as.numeric(arima.sim(n = n, model = list(
              order = c(1, 0, 12),
              ar = c(0.637),
              ma = c(0.667, 0.073, 0.434, 1.148, 0.586, 0.308, 0.265, 0.941,
                    0.797, 0.370, 0.333, 0.679)), sd = sqrt(0.0001))) + 0.650
      # Simulate IndPro
      v2 <- as.numeric(arima.sim(n = n, model = list(
              order = c(2, 0, 1),
              ar = c(0.575, 0.373),
              ma = c(0.706)
            ), sd = sqrt(0.00002)))
      interim_mat <- cbind(const, v1, v2)
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
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Defines data generation functions for simulations.")
  p <- add_argument(p, "--output", type = "character", nargs = 1,
                    default = "SimDataExample.R",
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

