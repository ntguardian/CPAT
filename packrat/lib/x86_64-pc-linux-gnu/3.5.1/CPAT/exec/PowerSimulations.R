#!/usr/bin/Rscript
################################################################################
# PowerSimulations.R
################################################################################
# 2018-09-04
# Curtis Miller
################################################################################
# Create simulated statistic values using simulated data
################################################################################

# This is code for which I don't want to worry about complete compliance with
# the style guide or necessarily good R practice, simply because this code is
# old and it works. I could try to functionalize and vectorize and create a
# main() function and separate parameters from execution and so on, but this
# works and it's not difficult to maintain. If I ever need to make a big change
# (such as add another test statistic), I will do all that work, but not today.
# So be careful when using source() on this file.

# TODO: curtis: ALLOW FOR COMMAND LINE STATISTIC SELECTION -- Wed 05 Sep 2018

# Command line functionality
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# COMMAND LINE INTERFACE
################################################################################

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`

cl_args <- parse_args(OptionParser(
    description = "Test statistic simulation script",
    option_list = list(
      make_option(c("--replications", "-N"), type = "integer", default = 5000,
                  help = "Number of replications per point"),
      make_option(c("--seed", "-s"), type = "integer", default = 20180904,
                  help = "The seed of the simulations"),
      make_option(c("--seedless", "-R"), action = "store_true",
                  help = "Don't set a seed (causes --seed to be ignored)"),
      make_option(c("--prefix", "-p"), type = "character", default = "",
                  help = "Prefix to attach to save names (useful for saving" %s%
                    "in directories)"),
      make_option(c("--infile", "-i"), type = "character",
                  help = "Input .Rda file with parameter definitions"),
      make_option(c("--outfile", "-o"), type = "character", default = "out.csv",
                  help = "Location to save CSV table with file names and" %s%
                    "statistic labels")
    )
))

################################################################################
# SETUP
################################################################################

# library(smoothmest)  # For rdoublex()
library(fGarch)
library(doParallel)
library(doRNG)

sim_Vn_stat <- CPAT:::sim_Vn_stat
sim_Zn_stat <- CPAT:::sim_Zn_stat
sim_de_stat <- CPAT:::sim_de_stat
sim_hs_stat <- CPAT:::sim_hs_stat

registerDoParallel(cores = detectCores())

if (!cl_args$seedless) {
  registerDoRNG(cl_args$seed)
}

sims <- cl_args$replications

out_df <- data.frame(file = character(0), statistic = character(0),
                     stringsAsFactors = FALSE)

prefix <- cl_args$prefix
infile <- cl_args$infile

load(infile)

for (f in c("power_simulations", "distfunc", "use_lrv", "nval", "delta",
            "kstarfunc")) {
  if (!exists(f)) stop("Object" %s% f %s% "was not defined in" %s% infile)
}

################################################################################
# LOOP
################################################################################

for (distname in names(distfunc)) {
  distargs <- distfunc[[distname]]
  # Prefixes will be applied to these later
  Vn_filename <- paste(distname, "Vn.Rda", sep = "_")
  # Vn_w3_filename <- paste(distname, "Vn_weight_3rd.Rda", sep = "_")
  # Vn_t5_filename <- paste(distname, "Vn_trim_5perc.Rda", sep = "_")
  # Vn_t10_filename <- paste(distname, "Vn_trim_10perc.Rda", sep = "_")
  Zn_filename <- paste(distname, "Zn.Rda", sep = "_")
  de_filename <- paste(distname, "de.Rda", sep = "_")
  hs_filename <- paste(distname, "hs.Rda", sep = "_")
  
  for (n in nval) {
    nname <- paste("n", n, sep = "")
        
    for (kstarname in names(kstarfunc)) {
      kstar <- kstarfunc[[kstarname]]
      distargs$changepoint <- kstar(n)
          
      for (d in delta) {
        savename <- paste("d", d, collapse = "", sep = "_")
        if (length(power_simulations$Vn[[distname]][[nname]][[
            kstarname]][[savename]]) < sims) {
          distargs$mean2 <- distargs$mean1 + d
          print(paste("Simulating: Vn", distname, nname, kstarname,
                savename, sep = "$"))
          ptm <- proc.time()["elapsed"]
          power_simulations$Vn[[distname]][[nname]][[kstarname]][[
                savename]] <- sim_Vn_stat(
                size = sims, n = n, gen_func = rchangepoint,
                args = distargs, use_kernel_var = use_lrv[distname],
                parallel = TRUE)
          print(paste("Completed in", proc.time()["elapsed"] - ptm,
                "seconds"))
          saveobj <- power_simulations$Vn
          save(saveobj, file = prefix %s0% Vn_filename)
        }
        # if (length(power_simulations$Vn_weight_3rd[[distname]][[
        #     nname]][[kstarname]][[savename]]) < sims) {
        #   distargs$mean2 <- distargs$mean1 + d
        #   print(paste("Simulating: Vn_weight_3rd", distname, nname,
        #         kstarname, savename, sep = "$"))
        #   ptm <- proc.time()["elapsed"]
        #   power_simulations$Vn_weight_3rd[[distname]][[nname]][[
        #         kstarname]][[savename]] <- sim_Vn_stat(
        #         size = sims, n = n, gen_func = rchangepoint,
        #         args = distargs, tau = 1/3, 
        #         use_kernel_var = use_lrv[distname],
        #         parallel = TRUE)
        #   print(paste("Completed in", proc.time()["elapsed"] - ptm,
        #         "seconds"))
        #   saveobj <- power_simulations$Vn_weight_3rd
        #   save(saveobj, file = prefix %s0% Vn_w3_filename)
        # }
        # if (length(power_simulations$Vn_trim_5perc[[distname]][[
        #     nname]][[kstarname]][[savename]]) < sims) {
        #   distargs$mean2 <- distargs$mean1 + d
        #   print(paste("Simulating: Vn_trim_5perc", distname, nname,
        #     kstarname, savename, sep = "$"))
        #   ptm <- proc.time()["elapsed"]
        #   power_simulations$Vn_trim_5perc[[distname]][[nname]][[
        #     kstarname]][[savename]] <- sim_Vn_stat(
        #     size = sims, n = n, gen_func = rchangepoint,
        #     args = distargs, tau = 1/2, kn = function(n) {
        #     ceiling(0.025*n)}, use_kernel_var = use_lrv[distname],
        #     parallel = TRUE)
        #   print(paste("Completed in", proc.time()["elapsed"] - ptm,
        #     "seconds"))
        #   saveobj <- power_simulations$Vn_trim_5perc
        #   save(saveobj, file = prefix %s0% Vn_t5_filename)
        # }
        # if (length(power_simulations$Vn_trim_10perc[[distname]][[
        #     nname]][[kstarname]][[savename]]) < sims) {
        #   distargs$mean2 <- distargs$mean1 + d
        #   print(paste("Simulating: Vn_trim_10perc", distname, nname,
        #     kstarname, savename, sep = "$"))
        #   ptm <- proc.time()["elapsed"]
        #   power_simulations$Vn_trim_10perc[[distname]][[nname]][[
        #     kstarname]][[savename]] <- sim_Vn_stat(
        #     size = sims, n = n, gen_func = rchangepoint,
        #     args = distargs, tau = 1/2,
        #     kn = function(n) {ceiling(0.05*n)},
        #     use_kernel_var = use_lrv[distname],
        #     parallel = TRUE)
        #   print(paste("Completed in", proc.time()["elapsed"] - ptm,
        #     "seconds"))
        #   saveobj <- power_simulations$Vn_trim_10perc
        #   save(saveobj, file = prefix %s0% Vn_t10_filename)
        # }
        if (length(power_simulations$de[[distname]][[nname]][[
            kstarname]][[savename]]) < sims) {
          distargs$mean2 <- distargs$mean1 + d
          print(paste("Simulating: de", distname, nname, kstarname,
              savename, sep = "$"))
          ptm <- proc.time()["elapsed"]
          power_simulations$de[[distname]][[nname]][[kstarname]][[
              savename]] <- sim_de_stat(
              size = sims, a = function(n) {log(n/(log(n)^(3/2)))},
              b = function(n) {log(n/(log(n)^(3/2)))},
              n = n, gen_func = rchangepoint, args = distargs,
              use_kernel_var = use_lrv[distname],
              parallel = TRUE)
          print(paste("Completed in", proc.time()["elapsed"] - ptm,
              "seconds"))
          saveobj <- power_simulations$de
          save(saveobj, file = prefix %s0% de_filename)
        }
        if (length(power_simulations$hs[[distname]][[nname]][[
            kstarname]][[savename]]) < sims) {
          distargs$mean2 <- distargs$mean1 + d
          print(paste("Simulating: hs", distname, nname, kstarname,
              savename, sep = "$"))
          ptm <- proc.time()["elapsed"]
          power_simulations$hs[[distname]][[nname]][[kstarname]][[
              savename]] <- sim_hs_stat(
            size = sims, n = n, gen_func = rchangepoint,
            args = distargs, corr = use_lrv[distname],
            parallel = TRUE)
          print(paste("Completed in", proc.time()["elapsed"] - ptm,
              "seconds"))
          saveobj <- power_simulations$hs
          save(saveobj,file = prefix %s0% hs_filename)
        }
        
        for (knname in names(knfunc)) {
          kn <- knfunc[[knname]]
          if (length(power_simulations$Zn[[distname]][[knname]][[
              nname]][[kstarname]][[savename]]) < sims &
            !(knname == "sqrt" & kstarname == "c4rt")) {
            distargs$mean2 <- distargs$mean1 + d
            print(paste("Simulating: Zn", distname, knname, nname,
                  kstarname, savename, sep = "$"))
            ptm <- proc.time()["elapsed"]
            power_simulations$Zn[[distname]][[knname]][[nname]][[
                kstarname]][[savename]] <-
                sim_Zn_stat(size = sims, kn = kn,
                      n = n, gen_func = rchangepoint,
                      args = distargs,
                      use_kernel_var = use_lrv[distname],
                      parallel = TRUE)
            print(paste("Completed in",
                proc.time()["elapsed"] - ptm, "seconds"))
            saveobj <- power_simulations$Zn
            save(saveobj, file = prefix %s0% Zn_filename)
          }
        }
      }
    }
  }

  power_simulations$Vn[distname] <- NULL
  #power_simulations$Vn_weight_3rd[distname] <- NULL
  #power_simulations$Vn_trim_5perc[distname] <- NULL
  #power_simulations$Vn_trim_10perc[distname] <- NULL
  power_simulations$Zn[distname] <- NULL
  power_simulations$de[distname] <- NULL
  power_simulations$hs[distname] <- NULL

  out_df <- rbind(out_df, data.frame(
    "file" = c(Vn_filename, Zn_filename, de_filename, hs_filename),
    "statistic" = c("Vn", "Zn", "de", "hs"),
    stringsAsFactors = FALSE))
  write.csv(out_df, file = cl_args$outfile)
}
