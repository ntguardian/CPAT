#!/usr/bin/Rscript
################################################################################
# LRVEstAnalysisParallel.R
################################################################################
# 2018-09-05
# Curtis Miller
################################################################################
# Perform simulations to assess the rate of convergence of the kernel-based LRV
# estimators.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

`%s%` <- CPAT:::`%s%`

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(file, replications = 10000, seed = 20180906, 
                 seedless = FALSE, help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work
  # See parameter descriptions with definition of cl_args.
  
  library(doParallel)
  library(fGarch)
  library(doRNG)

  get_lrv_vec <- CPAT:::get_lrv_vec

  registerDoParallel(cores = detectCores() - 1)
  if (!seedless) registerDoRNG(seed)

  `%s0%` <- function(x,y) {paste0(x,y)}

  N <- c(50, 100, 200, 500)
  ar_coef <- c(0, .05, .1, .2, .5, .7)

  ar_1_lrv <- function(sd, phi) {sd / (1 - phi)^2}
  theo_lrv <- sapply(ar_coef, function(phi) ar_1_lrv(1, phi))

  # AR(1) (phi > 0)
  # Bartlett kernel, Andrews bandwidth

  ar1_plotlist_pos_phi <- foreach(phi = ar_coef, .combine = rbind) %do% {
    foreach(n = N, .combine = rbind) %do% {
      print("Working for n = " %s0% n %s0% " and phi = " %s0% phi)
      lrvs <- foreach(i = 1:replications, .combine = c) %dopar% {
        get_lrv_vec(arima.sim(list(order = c(1, 0, 0), ar = phi), n = n,
                    n.start = 500))[round(n/2)]
      }
      data.frame(val = lrvs, n = rep(n, times = replications),
                 phi = rep(phi, times = replications))
    }
  }

  # AR(1) (phi < 0)
  # Bartlett kernel, Andrews bandwidth

  ar1_plotlist_neg_phi <- foreach(phi = -ar_coef, .combine = rbind) %do% {
    foreach(n = N, .combine = rbind) %do% {
      print("Working for n = " %s0% n %s0% " and phi = " %s0% phi)
      lrvs <- foreach(i = 1:replications, .combine = c) %dopar% {
        get_lrv_vec(arima.sim(list(order = c(1, 0, 0), ar = phi), n = n,
            n.start = 500))[round(n/2)]
      }
      data.frame(val = lrvs, n = rep(n, times = replications),
                 phi = rep(phi, times = replications))
    }
  }

  # GARCH(1,1)
  # Bartlett kernel, Andrews bandwidth

  garch_plotlist <- foreach(n = N, .combine = rbind) %do% {
    print("Working for n = " %s0% n)
    lrvs <- foreach(i = 1:replications, .combine = c) %dopar% {
      get_lrv_vec(garchSim(garchSpec(list(alpha = 0.1, beta = 0.1, omega = .5)),
          n.start = 500, n = n)$garch)[round(n/2)]
    }
    data.frame(val = lrvs, n = rep(n, times = replications))
  }

  # AR(1) (phi > 0)
  # Flat top kernel, Andrews bandwidth

  bartlett <- function(z) {ifelse(abs(z) < 1, 1 - abs(z), 0)}
  flattop <- function(x) {2 * bartlett(x) - bartlett(2*x)}

  ar1_ft_plotlist_pos_phi <- foreach(phi = ar_coef, .combine = rbind) %do% {
    foreach(n = N, .combine = rbind) %do% {
      print("Working for n = " %s0% n %s0% " and phi = " %s0% phi)
      lrvs <- foreach(i = 1:replications, .combine = c) %dopar% {
        get_lrv_vec(arima.sim(list(order = c(1, 0, 0), ar = phi), n = n,
            n.start = 500), kernel = flattop)[round(n/2)]
      }
      data.frame(val = lrvs, n = rep(n, times = replications),
                 phi = rep(phi, times = replications))
    }
  }

  # AR(1) (phi < 0)
  # Flat-top kernel, Andrews bandwidth

  ar1_ft_plotlist_neg_phi <- foreach(phi = -ar_coef, .combine = rbind) %do% {
    foreach(n = N, .combine = rbind) %do% {
      print("Working for n = " %s0% n %s0% " and phi = " %s0% phi)
      lrvs <- foreach(i = 1:replications, .combine = c) %dopar% {
        get_lrv_vec(arima.sim(list(order = c(1, 0, 0), ar = phi), n = n,
            n.start = 500), kernel = flattop)[round(n/2)]
      }
      data.frame(val = lrvs, n = rep(n, times = replications),
                 phi = rep(phi, times = replications))
    }
  }

  # GARCH(1,1)
  # Flat-top kernel, Andrews bandwidth

  garch_ft_plotlist <- foreach(n = N, .combine = rbind) %do% {
    print("Working for n = " %s0% n)
    lrvs <- foreach(i = 1:replications, .combine = c) %dopar% {
      get_lrv_vec(garchSim(garchSpec(list(alpha = 0.1, beta = 0.1, omega = .5)),
          n.start = 500, n = n)$garch, kernel = flattop)[round(n/2)]
    }
    data.frame(val = lrvs, n = rep(n, times = replications))
  }

  save(ar1_plotlist_pos_phi, ar1_plotlist_neg_phi, garch_plotlist, 
       ar1_ft_plotlist_pos_phi, ar1_ft_plotlist_neg_phi, garch_ft_plotlist,
       file = file)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Perform LRV estimator convergence rate simulations.",
        option_list = list(
          make_option(c("--replications", "-N"), type = "integer",
                      default = 10000,
                      help = "Number of replications to perform"),
          make_option(c("--seed", "-s"), type = "integer",
                      default = 20180906,
                      help = "RNG seed"),
          make_option(c("--seedless", "-R"), action = "store_true",
                      help = "Don't set a seed (causes --seed to be ignored)"),
          make_option(c("--file", "-f"), type = "character",
                      help = "Name of .Rda output file")
        )
      ))

  do.call(main, cl_args)
}

