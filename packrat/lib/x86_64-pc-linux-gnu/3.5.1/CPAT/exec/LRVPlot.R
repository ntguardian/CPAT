#!/usr/bin/Rscript
################################################################################
# LRVPlot.R
################################################################################
# 2018-09-07
# Curtis Miller
################################################################################
# Create plots for simulations of kernel-based long-run variance estimates.
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

main <- function(file, outfile = "lrv_est_plot", verbose = FALSE, width = 4.5,
                 height = 3.5, help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work
  # help parameter is ignored; see cl_args for parameter description

  library(dplyr)

  lrv_plot_tikz <- CPAT:::lrv_plot_tikz

  load(file)

  if (!(exists("ar1_plotlist_pos_phi") &
        exists("ar1_plotlist_neg_phi") &
        exists("ar1_ft_plotlist_pos_phi") &
        exists("ar1_ft_plotlist_neg_phi") &
        exists("garch_plotlist") &
        exists("garch_ft_plotlist"))) {
    stop("Not all objects expected in loaded .Rda file were found! Must" %s%
         "have ar1_plotlist_pos_phi, ar1_plotlist_neg_phi, garch_plotlist," %s%
         "ar1_ft_plotlist_pos_phi, ar1_ft_plotlist_neg_phi, garch_ft_plotlist")
  }
  
  lrv_n <- unique(ar1_plotlist_pos_phi$n)
  lrv_phi <- unique(ar1_plotlist_pos_phi$phi)
  lrv_phi_neg <- unique(ar1_plotlist_neg_phi$phi)

  for (n in lrv_n) {
    for (p in lrv_phi) {
      lrv_plot_tikz(ar1_plotlist_pos_phi, n, p, ker_name = "bartlett",
                    true_lrv = get_theo_lrv(1, p), xrange = c(0, 20),
                    filename = paste(outfile, "bartlett",
                      "ar" %s0% "_" %s0% p, n, sep = "_"))
      lrv_plot_tikz(ar1_ft_plotlist_pos_phi, n, p, ker_name = "flattop",
                    true_lrv = get_theo_lrv(1, p), xrange = c(0, 20),
                    filename = paste(outfile, "flattop",
                      "ar" %s0% "_" %s0% p, n, sep = "_"))
    }
  }

  for (n in lrv_n) {
    for (p in lrv_phi_neg) {
      lrv_plot_tikz(ar1_plotlist_neg_phi, n, p, ker_name = "bartlett",
                    true_lrv = get_theo_lrv(1, p), xrange = c(0, 5),
                    filename = paste(outfile, "bartlett",
                      "ar" %s0% "_" %s0% p, n, sep = "_"))
      lrv_plot_tikz(ar1_ft_plotlist_neg_phi, n, p, ker_name = "flattop",
                    true_lrv = get_theo_lrv(1, p), xrange = c(0, 5),
                    filename = paste(outfile, "flattop",
                      "ar" %s0% "_" %s0% p, n, sep = "_"))
    }
  }

  for (n in lrv_n) {
    lrv_plot_tikz(garch_plotlist, n, ker_name = "bartlett",
                  true_lrv = garch_1_1_lrv(.1, .1, .5), xrange = c(0, 2),
                    filename = paste(outfile, "bartlett", "garch", n, sep = "_"))
    lrv_plot_tikz(garch_ft_plotlist, n, ker_name = "flattop",
                  true_lrv = garch_1_1_lrv(.1, .1, .5), xrange = c(0, 2),
                    filename = paste(outfile, "flattop", "garch", n, sep = "_"))
  }
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Create plots for simulations of kernel-based" %s% 
                      "long-run variance estimates.",
        option_list = list(
          make_option(c("--file", "-f"), type = "character",
                      help = "The name of the .Rda file containing LRV" %s%
                             "simuation data objects"),
          make_option(c("--outfile", "-o"), type = "character",
                      help = "Stem of output file names; can include path" %s%
                             "information"),
          make_option(c("--verbose", "-v"), action = "store_true",
                      help = "Print updates about progress"),
          make_option(c("--width", "-w"), type = "double", default = 4.5,
                      help = "Width of the resulting plot"),
          make_option(c("--height", "-H"), type = "double", default = 3.5,
                      help = "Height of the resulting plot")
        )
      ))

  do.call(main, cl_args)
}

