#!/usr/bin/Rscript
################################################################################
# DistConvPlot.R
################################################################################
# 2018-09-07
# Curtis Miller
################################################################################
# Create Tikz/PDF plots demonstrating convergence in distribution of the
# Rényi-type statistic
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}
if (!requireNamespace("CPAT", quietly = TRUE)) stop("Install CPAT")

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(file, outfile = "dist_conv", verbose = FALSE, width = 4,
                 height = 3, help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work
  # help parameter does nothing

  load(file)
  if (!exists(Zn_simulations)) stop("The object Zn_simulations must be in" %s%
                                    "loaded file")

  dZn <- CPAT:::dZn
  dist_conv_plot_tikz <- CPAT:::dist_conv_plot_tikz

  # filename <- paste0("dist_conv_", dist, "_n", size, "_", trim)
  filename <- gsub(".", "", filename, fixed = TRUE)
  filename.tex <- paste0(filename, ".tex")
  filename.pdf <- paste0(filename, ".pdf")

  conv_T <- as.integer(substr(names(Zn_simulations[[1]][[1]]), 2, 100))
  conv_dist <- names(Zn_simulations)
  conv_trim <- names(Zn_simulations[[1]])

  for (T in conv_T) {
    for (d in conv_dist) {
      for (t in conv_trim) {
        # dist_conv_plot_tikz(d, t, T, title = paste(paste0("$T = ", T, "$"),
        # conv_dstring[d], conv_tstring[t], sep = ", "))
        dist_conv_plot_tikz(d, t, T, title = paste0("$T = ", T, "$"),
                            verbose = verbose,
                            filename = paste0(outfile, "_", d, "_n", T, "_", t),
                            makePDF = TRUE, width = width, height = height)
      }
    }
  }
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Create PDF plots for demonstrating the Rényi-type" %s%
                      "statistic's convergence in distribution under the" %s%
                      "null hypothesis.",
        option_list = list(
          make_option(c("--file", "-f"), type = "character",
                      help = ".Rda file to load data from"),
          make_option(c("--outfile", "-o"), type = "character",
                      help = "Prefix string of output file" %s%
                             "(not including extension; will be of the form" %s%
                             "outfile_dist_n##_trim.pdf; can include file" %s%
                             "path information)"),
          make_option(c("--verbose", "-v"), action = "store_true",
                      help = "Print to output progress reports"),
          make_option(c("--width", "-w"), type = "double", default = 4,
                      help = "Width of the resulting plot (inches)"),
          make_option(c("--height", "-H"), type = "double", default = 3,
                      help = "Height of the resulting plot")
        )
      ))

  do.call(main, cl_args)
}

