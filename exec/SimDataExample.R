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
    n1 <- min(n, 156)
    n2 <- n - n1
    vec1_innov <- rnorm(n1, sd = sqrt(2.230))
    vec1 <- as.numeric(arima.sim(n = n1, model = list(
                        order = c(2, 0, 0),
                        ar = c(1.323, -0.341)
                      ), innov = vec1_innov))
    if (n2 == 0) {
      vec2 <- numeric()
    } else {
      vec2 <- as.numeric(arima.sim(n = n2, model = list(
                          order = c(2, 0, 0),
                          ar = c(1.323, -0.341)
                        ), sd = sqrt(21.941), n.start = length(vec1_innov),
                        start.innov = vec1_innov))
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
              order = c(1, 0, 2),
              ar = c(0.976),
              ma = c(0.373, 0.207)), sd = sqrt(15.609))) + 43.991
      # Simulate IndPro
      v2 <- as.numeric(arima.sim(n = n, model = list(
              order = c(1, 0, 1),
              ar = c(0.985),
              ma = c(0.367)
            ), sd = sqrt(1.323))) + 100
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

