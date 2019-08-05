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

  ##############################################################################
  # REQUIRED OBJECTS
  ##############################################################################

  eps_generator <- function(n) {
    n1 <- min(n, 216)
    n2 <- n - n1
    x1 <- as.numeric(arima.sim(n = n1, n.start = 500, model = list(
                                 order = c(1, 0, 2), ar = c(0.941),
                                 ma = c(-1.160, 0.249)), sd = 3.438))
    if (n2 > 0) {
      x2 <- as.numeric(arima.sim(n = n2, n.start = n1, start.innov = x1,
                                 model = list(order = c(2, 0, 1),
                                              ar = c(-0.989, -0.370),
                                              ma = c(0.776)), sd = 2.581))
    } else {
      x2 <- numeric()
    }

    x <- c(x1, x2)

    return(x)
  }

  df_generator <- function(n, beta, eps) {
    d <- length(beta)
    stopifnot(d == 6)
    const <- rep(1, times = n)
    interim_mat <- cbind(const,
                         as.numeric(arima.sim(model = list(
                               ar = c(0.893, 0.059),
                               ma = c(-1.000)),
                             sd = 0.699, n = n)) + 0.0627,       # Mkt.RF
                         rnorm(n, mean =  0.0120, sd = 0.5163),  # SMB
                         as.numeric(arima.sim(model = list(
                               ar = c(0.136)
                               ), n = n)) + 0.0173,              # HML
                         rnorm(n, mean =  0.0121, sd = 0.3350),  # RMW
                         rnorm(n, mean = -0.0068, sd = 0.3594))  # CMA

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

