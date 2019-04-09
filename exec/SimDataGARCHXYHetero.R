#!/usr/bin/Rscript
################################################################################
# SimDataGARCHXYHetero.R
################################################################################
# 2019-03-01
# Curtis Miller
################################################################################
# Data simulation definition file; simulate data from Normal distribution, and
# residuals from a GARCH(1, 1) process, with a change in variance of the data
# mid-sample.
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

  spec1 <- rugarch::ugarchspec(mean.model = list(
                                 armaOrder = c(0, 0),
                                 include.mean = FALSE
                               ), fixed.pars = list(
                                 "omega" = 0.2,
                                 "alpha1" = 0.3,
                                 "beta1" = 0.5
                               ))
  spec2 <- rugarch::ugarchspec(mean.model = list(
                                 armaOrder = c(0, 0),
                                 include.mean = FALSE
                               ), fixed.pars = list(
                                 "omega" = 1,
                                 "alpha1" = 0.7,
                                 "beta1" = 0.2
                               ))

  ##############################################################################
  # REQUIRED OBJECTS
  ##############################################################################

  # GARCH(1, 1) model
  eps_generator <- function(n) {
    n1 <- ceiling(n / 2)
    n2 <- n - n1
    x_obj_1 <- rugarch::ugarchpath(spec1, n.sim = n1, n.start = 500,
                                   # Must set seed for this function; don't
                                   # trust without this
                                   rseed = sample(1:9999999, 1))
    x_obj_2 <- rugarch::ugarchpath(spec2, n.sim = n2, n.start = 500,
                                   # Must set seed for this function; don't
                                   # trust without this
                                   rseed = sample(1:9999999, 1))
    c(as.numeric(x_obj_1@path$seriesSim), as.numeric(x_obj_2@path$seriesSim))
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
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Defines data generation functions for simulations.")
  p <- add_argument(p, "--output", type = "character", nargs = 1,
                    default = "SimDataGARCHXY.R",
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

