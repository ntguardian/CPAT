#!/usr/bin/Rscript
################################################################################
# PowerPlot.R
################################################################################
# 2018-09-07
# Curtis Miller
################################################################################
# Construct plots for power simulations.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(file, outfile = "power_sim", by_stat = FALSE, verbose = FALSE, 
                 statlines = "c('Zn' = 'solid', 'Vn' = 'longdash'," %s%
                             "'de' = 'dotdash', 'hs' = 'F24242')",
                 help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work
  # The help parameter is not used here, but is needed; see cl_args definition
  # for more information about parameters

  dat <- read.csv(file)
  if (!identical(sort(names(dat)), sort(c("power", "stat","dist", "kn", "n",
          "cpt", "delta")))) {
    stop("Invalid column names for" %s% file %s0% "; they should be power," %s%
         "stat, dist, kn, n, cpt, and delta")
  }
  if (!is.numeric(dat$power)) {
    stop("The column power in" %s% file %s% "could not be interpreted as" %s%
         "values")
  }
  if (!is.integer(dat$n)) stop("The n column in" %s% file %s% "could not be" %s%
                               "understood as integers")

  power_plot_tikz_by_n <- CPAT:::power_plot_tikz_by_n
  power_plot_tikz <- CPAT:::power_plot_tikz

  statlines <- eval(parse(text = statlines))

  power_dist <- unique(dat$dist)
  power_cpt <- unique(dat$cpt)
  power_stat <- unique(dat$stat)
  power_kn <- unique(dat$kn)
  power_n <- unique(dat$n)

  for (c in power_cpt) {
    for (d in power_dist) {
      for (k in power_kn) {
        if (by_stat) {
          for (s in power_stat) {
            power_plot_tikz(dat, d, k, c, s,
                            filename = paste(outfile, d, s, k, c, sep = "_"),
                            verbose = verbose)
          }
        } else {
          for (n in power_n) {
            power_plot_tikz_by_n(dat, d, k, c, n, statlines = statlines,
                                 filename = paste(outfile, d, "n" %s0% n, k, c,
                                                  sep = "_"),
                                 verbose = verbose)
          }
        }
      }
    }
  }
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Construct plots for power simulations.",
        option_list = list(
          make_option(c("--file", "-f"), type = "character",
                      help = "The CSV file from which to load data"),
          make_option(c("--outfile", "-o"), type = "character",
                      default = "power_plot",
                      help = "Prefix for output file names (could include" %s%
                             "path information)"),
          make_option(c("--by-stat", "-s"), action = "store_true",
                      help = "Plots should be split by statistic, not" %s%
                             "sample size"),
          make_option(c("--verbose", "-v"), action = "store_true",
                      help = "Print progress reports"),
          make_option(c("--statlines", "-t"), type = "character",
                      default = "c('Zn' = 'solid', 'Vn' = 'longdash'," %s%
                                "'de' = 'dotdash', 'hs' = 'F24242')",
                      help = "An R-parseable string creating a named" %s%
                             "character vector where names are labels for" %s%
                             "statistics and values are identifiers for" %s%
                             "line types ggplot2 can understand")
        )
      ))

  names(cl_args)[which(names(cl_args) == c("by-stat"))] <- "by_stat"
  do.call(main, cl_args)
}

