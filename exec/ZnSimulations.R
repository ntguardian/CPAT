#!/usr/bin/Rscript
################################################################################
# ZnSimulations.R
################################################################################
# 2018-09-05
# Curtis Miller
################################################################################
# Generate simulations of the Rényi-type statistic to demonstrate convergence to
# the limiting distribution.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(file, replications = 100000, seed = 20180906,
                 seedless = FALSE, verbose = FALSE, help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work.
  # See definition of cl_args for parameter descriptions.
  library(smoothmest) # For rdoublex()

  dZn <- CPAT:::dZn
  sim_Zn <- CPAT:::sim_Zn

  if (!seedless) {set.seed(seed)}

  # This code is very old; I know it's bad, but it works.
  plusminusflip = function(n) {return(sample(c(-1, 1),
                                             size = n,
                                             replace = T))}

  Zn_simulations = list(norm = list(log = list(),
                                   sqrt = list(),
                                   rt4 = list(),
                                   rt16 = list()),
                       plusminusflip = list(log = list(),
                                   sqrt = list(),
                                   rt4 = list(),
                                   rt16 = list()),
                       exp = list(log = list(),
                                   sqrt = list(),
                                   rt4 = list(),
                                   rt16 = list()),
                       doubleexp = list(log = list(),
                                   sqrt = list(),
                                   rt4 = list(),
                                   rt16 = list()),
                       gamma = list(log = list(),
                                   sqrt = list(),
                                   rt4 = list(),
                                   rt16 = list()),
                       t3 = list(log = list(),
                                   sqrt = list(),
                                   rt4 = list(),
                                   rt16 = list()))

  if (verbose) {cat("Simulations for norm")}
  Zn_simulations$norm$log$n50 <- sim_Zn(replications, log, n = 50)
  Zn_simulations$norm$sqrt$n50 <- sim_Zn(replications, sqrt, n = 50)
  Zn_simulations$norm$rt4$n50 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))},
                               n = 50)
  Zn_simulations$norm$rt16$n50 <- sim_Zn(replications,
                               function(n) {return(n^(1/16))},
                               n = 50)
  Zn_simulations$norm$log$n250 <- sim_Zn(replications, log,
                                        n = 250)
  Zn_simulations$norm$sqrt$n250 <- sim_Zn(replications, sqrt,
                                         n = 250)
  Zn_simulations$norm$rt4$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 250)
  Zn_simulations$norm$rt16$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 250)
  Zn_simulations$norm$log$n500 <- sim_Zn(replications, log)
  Zn_simulations$norm$sqrt$n500 <- sim_Zn(replications, sqrt)
  Zn_simulations$norm$rt4$n500 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))})
  Zn_simulations$norm$rt16$n500 <- sim_Zn(replications,
                               function(n) {return(n^(1/16))})

  Zn_simulations$norm$log$n1000 <- sim_Zn(replications, log,
                                        n = 1000)
  Zn_simulations$norm$sqrt$n1000 <- sim_Zn(replications, sqrt,
                                         n = 1000)
  Zn_simulations$norm$rt4$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 1000)
  Zn_simulations$norm$rt16$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 1000)

  # Alternate distribution: +/- 1
  if (verbose) {cat("Simulations for plusminusflip")}
  Zn_simulations$plusminusflip$log$n50 <- sim_Zn(replications, log,
                                    gen_func = plusminusflip, n = 50)
  Zn_simulations$plusminusflip$sqrt$n50 <- sim_Zn(replications,sqrt,
                                    gen_func = plusminusflip, n = 50)
  Zn_simulations$plusminusflip$rt4$n50 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))},
                               gen_func = plusminusflip, n = 50)
  Zn_simulations$plusminusflip$rt16$n50 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              gen_func = plusminusflip, n = 50)
  Zn_simulations$plusminusflip$log$n250 <- sim_Zn(replications, log,
                                  n = 250,
                                  gen_func = plusminusflip)
  Zn_simulations$plusminusflip$sqrt$n250 <-sim_Zn(replications,sqrt,
                                  n = 250,
                                  gen_func = plusminusflip)
  Zn_simulations$plusminusflip$rt4$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 250,
                              gen_func = plusminusflip)
  Zn_simulations$plusminusflip$rt16$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 250,
                              gen_func = plusminusflip)
  Zn_simulations$plusminusflip$log$n500 <- sim_Zn(replications, log,
                                    gen_func = plusminusflip)
  Zn_simulations$plusminusflip$sqrt$n500 <- sim_Zn(replications,sqrt,
                                    gen_func = plusminusflip)
  Zn_simulations$plusminusflip$rt4$n500 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))},
                               gen_func = plusminusflip)
  Zn_simulations$plusminusflip$rt16$n500 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              gen_func = plusminusflip)
  Zn_simulations$plusminusflip$log$n1000 <- sim_Zn(replications, log,
                                  n = 1000,
                                  gen_func = plusminusflip)
  Zn_simulations$plusminusflip$sqrt$n1000 <-sim_Zn(replications,sqrt,
                                  n = 1000,
                                  gen_func = plusminusflip)
  Zn_simulations$plusminusflip$rt4$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 1000,
                              gen_func = plusminusflip)
  Zn_simulations$plusminusflip$rt16$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 1000,
                              gen_func = plusminusflip)

  # Alternate distribution: exponential
  if (verbose) {cat("Simulations for exp")}
  Zn_simulations$exp$log$n50 <- sim_Zn(replications, log,
                                    gen_func = rexp, n = 50)
  Zn_simulations$exp$sqrt$n50 <- sim_Zn(replications,sqrt,
                                    gen_func = rexp, n = 50)
  Zn_simulations$exp$rt4$n50 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))},
                               gen_func = rexp, n = 50)
  Zn_simulations$exp$rt16$n50 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              gen_func = rexp, n = 50)
  Zn_simulations$exp$log$n250 <- sim_Zn(replications, log,
                                  n = 250,
                                  gen_func = rexp)
  Zn_simulations$exp$sqrt$n250 <-sim_Zn(replications,sqrt,
                                  n = 250,
                                  gen_func = rexp)
  Zn_simulations$exp$rt4$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 250,
                              gen_func = rexp)
  Zn_simulations$exp$rt16$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 250,
                              gen_func = rexp)
  Zn_simulations$exp$log$n500 <- sim_Zn(replications, log,
                                    gen_func = rexp)
  Zn_simulations$exp$sqrt$n500 <- sim_Zn(replications,sqrt,
                                    gen_func = rexp)
  Zn_simulations$exp$rt4$n500 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))},
                               gen_func = rexp)
  Zn_simulations$exp$rt16$n500 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              gen_func = rexp)

  Zn_simulations$exp$log$n1000 <- sim_Zn(replications, log,
                                  n = 1000,
                                  gen_func = rexp)
  Zn_simulations$exp$sqrt$n1000 <-sim_Zn(replications,sqrt,
                                  n = 1000,
                                  gen_func = rexp)
  Zn_simulations$exp$rt4$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 1000,
                              gen_func = rexp)
  Zn_simulations$exp$rt16$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 1000,
                              gen_func = rexp)

  # Alternate distribution: gamma
  if (verbose) {cat("Simulations for gamma")}
  Zn_simulations$gamma$log$n50 <- sim_Zn(replications, log,
                                    gen_func = rgamma,
                                    args = list(shape = 1), n = 50)
  Zn_simulations$gamma$sqrt$n50 <- sim_Zn(replications,sqrt,
                                    gen_func = rgamma,
                                    args = list(shape = 1), n = 50)
  Zn_simulations$gamma$rt4$n50 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))},
                               gen_func = rgamma,
                               args = list(shape = 1), n = 50)
  Zn_simulations$gamma$rt16$n50 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              gen_func = rgamma,
                              args = list(shape = 1), n = 50)
  Zn_simulations$gamma$log$n250 <- sim_Zn(replications, log,
                                  n = 250,
                                  gen_func = rgamma,
                                  args = list(shape = 1))
  Zn_simulations$gamma$sqrt$n250 <-sim_Zn(replications,sqrt,
                                  n = 250,
                                  gen_func = rgamma,
                                  args = list(shape = 1))
  Zn_simulations$gamma$rt4$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 250,
                              gen_func = rgamma,
                              args = list(shape = 1))
  Zn_simulations$gamma$rt16$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 250,
                              gen_func = rgamma,
                              args = list(shape = 1))
  Zn_simulations$gamma$log$n500 <- sim_Zn(replications, log,
                                    gen_func = rgamma,
                                    args = list(shape = 1))
  Zn_simulations$gamma$sqrt$n500 <- sim_Zn(replications,sqrt,
                                    gen_func = rgamma,
                                    args = list(shape = 1))
  Zn_simulations$gamma$rt4$n500 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))},
                               gen_func = rgamma,
                               args = list(shape = 1))
  Zn_simulations$gamma$rt16$n500 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              gen_func = rgamma,
                              args = list(shape = 1))

  Zn_simulations$gamma$log$n1000 <- sim_Zn(replications, log,
                                  n = 1000,
                                  gen_func = rgamma,
                                  args = list(shape = 1))
  Zn_simulations$gamma$sqrt$n1000 <-sim_Zn(replications,sqrt,
                                  n = 1000,
                                  gen_func = rgamma,
                                  args = list(shape = 1))
  Zn_simulations$gamma$rt4$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 1000,
                              gen_func = rgamma,
                              args = list(shape = 1))
  Zn_simulations$gamma$rt16$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 1000,
                              gen_func = rgamma,
                              args = list(shape = 1))


  # Alternate distribution: double exponential
  if (verbose) {cat("Simulations for doubleexp")}
  Zn_simulations$doubleexp$log$n50 <- sim_Zn(replications, log,
                                    gen_func = rdoublex,
                                    sd = 1/sqrt(2), n = 50)
  Zn_simulations$doubleexp$sqrt$n50 <- sim_Zn(replications,sqrt,
                                    gen_func = rdoublex,
                                    sd = 1/sqrt(2), n = 50)
  Zn_simulations$doubleexp$rt4$n50 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))},
                               gen_func = rdoublex,
                                    sd = 1/sqrt(2), n = 50)
  Zn_simulations$doubleexp$rt16$n50 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              gen_func = rdoublex,
                                    sd = 1/sqrt(2), n = 50)
  Zn_simulations$doubleexp$log$n250 <- sim_Zn(replications, log,
                                  n = 250,
                                  gen_func = rdoublex,
                                    sd = 1/sqrt(2))
  Zn_simulations$doubleexp$sqrt$n250 <-sim_Zn(replications,sqrt,
                                  n = 250,
                                  gen_func = rdoublex,
                                    sd = 1/sqrt(2))
  Zn_simulations$doubleexp$rt4$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 250,
                              gen_func = rdoublex,
                                    sd = 1/sqrt(2))
  Zn_simulations$doubleexp$rt16$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 250,
                              gen_func = rdoublex,
                                    sd = 1/sqrt(2))
  Zn_simulations$doubleexp$log$n500 <- sim_Zn(replications, log,
                                    gen_func = rdoublex,
                                    sd = 1/sqrt(2))
  Zn_simulations$doubleexp$sqrt$n500 <- sim_Zn(replications,sqrt,
                                    gen_func = rdoublex,
                                    sd = 1/sqrt(2))
  Zn_simulations$doubleexp$rt4$n500 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))},
                               gen_func = rdoublex,
                                    sd = 1/sqrt(2))
  Zn_simulations$doubleexp$rt16$n500 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              gen_func = rdoublex,
                                    sd = 1/sqrt(2))
  Zn_simulations$doubleexp$log$n1000 <- sim_Zn(replications, log,
                                  n = 1000,
                                  gen_func = rdoublex,
                                    sd = 1/sqrt(2))
  Zn_simulations$doubleexp$sqrt$n1000 <-sim_Zn(replications,sqrt,
                                  n = 1000,
                                  gen_func = rdoublex,
                                    sd = 1/sqrt(2))
  Zn_simulations$doubleexp$rt4$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 1000,
                              gen_func = rdoublex,
                                    sd = 1/sqrt(2))
  Zn_simulations$doubleexp$rt16$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 1000,
                              gen_func = rdoublex,
                                    sd = 1/sqrt(2))


  # Alternate distribution: t(3)
  if (verbose) {cat("Simulations for t3")}
  Zn_simulations$t3$log$n50 <- sim_Zn(replications, log,
                                    gen_func = rt,
                                    args = list(df = 3),
                                    sd = sqrt(3), n = 50)
  Zn_simulations$t3$sqrt$n50 <- sim_Zn(replications,sqrt,
                                    gen_func = rt,
                                    args = list(df = 3),
                                    sd = sqrt(3), n = 50)
  Zn_simulations$t3$rt4$n50 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))},
                               gen_func = rt,
                               args = list(df = 3),
                               sd = sqrt(3), n = 50)
  Zn_simulations$t3$rt16$n50 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              gen_func = rt,
                              args = list(df = 3),
                              sd = sqrt(3), n = 50)
  Zn_simulations$t3$log$n250 <- sim_Zn(replications, log,
                                  n = 250,
                                  gen_func = rt,
                                  args = list(df = 3),
                                  sd = sqrt(3))
  Zn_simulations$t3$sqrt$n250 <-sim_Zn(replications, sqrt,
                                  n = 250,
                                  gen_func = rt,
                                  args = list(df = 3),
                                  sd = sqrt(3))
  Zn_simulations$t3$rt4$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 250,
                              gen_func = rt,
                              args = list(df = 3),
                              sd = sqrt(3))
  Zn_simulations$t3$rt16$n250 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 250,
                              gen_func = rt,
                              args = list(df = 3),
                              sd = sqrt(3))
  Zn_simulations$t3$log$n500 <- sim_Zn(replications, log,
                                    gen_func = rt,
                                    args = list(df = 3),
                                    sd = sqrt(3))
  Zn_simulations$t3$sqrt$n500 <- sim_Zn(replications,sqrt,
                                    gen_func = rt,
                                    args = list(df = 3),
                                    sd = sqrt(3))
  Zn_simulations$t3$rt4$n500 <- sim_Zn(replications,
                               function(n) {return(n^(1/4))},
                               gen_func = rt,
                               args = list(df = 3),
                               sd = sqrt(3))
  Zn_simulations$t3$rt16$n500 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              gen_func = rt,
                              args = list(df = 3),
                              sd = sqrt(3))
  Zn_simulations$t3$log$n1000 <- sim_Zn(replications, log,
                                  n = 1000,
                                  gen_func = rt,
                                  args = list(df = 3),
                                  sd = sqrt(3))
  Zn_simulations$t3$sqrt$n1000 <-sim_Zn(replications, sqrt,
                                  n = 1000,
                                  gen_func = rt,
                                  args = list(df = 3),
                                  sd = sqrt(3))
  Zn_simulations$t3$rt4$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/4))},
                              n = 1000,
                              gen_func = rt,
                              args = list(df = 3),
                              sd = sqrt(3))
  Zn_simulations$t3$rt16$n1000 <- sim_Zn(replications,
                              function(n) {return(n^(1/16))},
                              n = 1000,
                              gen_func = rt,
                              args = list(df = 3),
                              sd = sqrt(3))

  save(Zn_simulations, file = file)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = paste("Generates simulated Rényi-type statistics under",
                            "different data generating processes and saves to",
                            "a file"),
        option_list = list(
          make_option(c("--file", "-f"), type = "character",
                      default = "ZnSimulations.Rda",
                      help = "Name of .Rda file to save Zn_simulations object"),
          make_option(c("--seed", "-s"), type = "integer", default = 20180906,
                      help = "Seed for RNG"),
          make_option(c("--seedless", "-R"), action = "store_true",
                      help = "Don't set a seed for random number generation"),
          make_option(c("--verbose", "-v"), action = "store_true",
                      help = "Messages during simulations"),
          make_option(c("--replications", "-r"), type = "integer",
                      default = 100000,
                      help = "Number of simulation replications")

        )
      ))

  do.call(main, cl_args)
}

