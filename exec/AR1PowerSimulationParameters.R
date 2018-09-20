#!/usr/bin/Rscript
################################################################################
# AR1PowerSimulationParameters.R
################################################################################
# 2018-09-05
# Curtis Miller
################################################################################
# Generate a .Rda file containing objects with data for parameters to be passed
# to PowerSimulations.R; this file containes AR(1) DGP settings.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# MAIN
################################################################################

main <- function(file, help = FALSE) {
  # Main function; see definition of cl_args for parameter definitions

  # This list will contain all simulated data and will be
  # saved for later use. The structure of the list is:
  #   Level 1: Statistic being simulated
  #   Level 2: Underlying distribution of data
  #   Level 2.5 (for Z_n): The function k_n
  #   Level 3: Underlying data sample size
  #   Level 4: Location of the changepoint
  #   Level 5: Change in mean
  # The final element being stored is a vector of
  # simulated statistics.

  power_simulations <- list(
    Vn = list(
      "ar1_0.5" = list(
        n50 = list(
          c4rt = list()
        )
      )
    ),
    # Vn_weight_3rd = list(
    #   "ar1_0.5" = list(
    #     n50 = list(
    #       c4rt = list()
    #     )
    #   )
    # ),
    # Vn_trim_5perc = list(
    #   "ar1_0.5" = list(
    #     n50 = list(
    #       c4rt = list()
    #     )
    #   )
    # ),
    # Vn_trim_10perc = list(
    #   "ar1_0.5" = list(
    #     n50 = list(
    #       c4rt = list()
    #     )
    #   )
    # ),
    de = list(
      "ar1_0.5" = list(
        n50 = list(
          c4rt = list()
        )
      )
    ),
    hs = list(
      "ar1_0.5" = list(
        n50 = list(
          c4rt = list()
        )
      )
    ),
    Zn = list(
      "ar1_0.5" = list(
        log = list(
          n50 = list(
            c4rt = list()
          )
        )
      )
    )
  )

  ##############################################################################
  # PARAMETER DEFINITIONS
  ##############################################################################

  # Underlying distributions considered:
  #     norm: Standard Normal
  #     t3: t distribution, 3 degrees of freedom
  #     doubleexp: Double exponential distribution with scale parameter equal to 1
  #     exp: Exponential distribution with theta = 1, shifting by delta
  #     Various flavors of GARCH
  #     AR(1)
  #     ARMA(1,1)
  distfunc <- list(
    "ar1_0.5" = list(
        mean1 = 0,
        dist = function(n, mu) {return(as.vector(arima.sim(
                 list(order = c(1,0,0), ar = c(0.5)),
                 n.start = 200, n = n)) + mu)},
        meanparam = "mu"
    )
  )

  # Specifies whether long-run variance estimation should be used
  use_lrv <- c(
    "ar1_0.5" = TRUE
  )

  # Alternative k_n considered:
  #     logsq: (log n)^2
  #     log: log n
  #     sqrt: n^(1/2)
  knfunc <- list(
    # "logsq" = function(n) {return(ceiling(log(n)^2))},
    "log" = function(n) {return(ceiling(log(n)))}#,
    # "sqrt" = function(n) {return(ceiling(sqrt(n)))}
  )

  # Alternative n considered:
  #     n50
  #     n100
  #     n200
  #     n500
  #nval <- c(50, 100, 200, 500)
  nval <- c(50, 200, 500)

  # Alternative k* considered:
  #     c4rt: n^(1/4) (not for k_n sqrt)
  #     chalf: n / 2
  kstarfunc <- list(
    "c4rt" = function(n) {return(ceiling(n^(1/4)))}#,
    # "chalf" = function(n) {return(ceiling(n/2))}
  )

  # Alternative delta considered:
  #     d_-2, d_-1.9,...,d_1.9,d_2: every d from -2 to 2 incrementing by .1
  delta <- seq(-2, 2, by = 0.1)


  save(power_simulations, distfunc, use_lrv, knfunc, nval, kstarfunc, delta,
       file = file)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = paste("Generate .Rda file containing parameter objects",
                            "for power simulations in PopwerSimulations.R"),
        option_list = list(
          make_option(c("--file", "-f"), type = "character",
                      default = "out.Rda",
                      help = "Name of the output file")
        )
      ))

  do.call(main, cl_args)
}

