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
    n1 <- min(n, 34)
    n2 <- n - n1
    vec1 <- as.numeric(rnorm(n1, sd = 1.856))
    if (n2 == 0) {
      vec2 <- numeric()
    } else {
      vec2 <- as.numeric(rnorm(n2, sd = 4.353))
    }
    c(vec1, vec2)
  }

  df_generator <- function(n, beta, eps) {
    d <- length(beta)
    stopifnot(d == 2)
    const <- rep(1, times = n)
    if (d > 1) {
      interim_mat <- matrix(rnorm(n * (d - 1), mean = 0.044, sd = 1.053),
                            ncol = d - 1)
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

