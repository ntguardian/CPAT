#!/usr/bin/Rscript
################################################################################
# LRVEstAnalysisPlot.R
################################################################################
# 2018-09-05
# Curtis Miller
################################################################################
# Plot the results of LRV estimation simulations
################################################################################

# XXX: curtis: DEPRECATED!!! -- Fri 07 Sep 2018 01:47:11 PM MDT

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(file, prefix = "lrv_est_", pdf = FALSE, width = 5,
                 height = 7, help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work
  # See cl_args for argument descriptions

  library(ggplot2)
  library(gridExtra)
  library(tikzDevice)

  `%s0%` <- CPAT:::`%s0%`

  load(file)

  for (f in c("ar1_plotlist_pos_phi", "ar1_plotlist_neg_phi", "garch_plotlist",
              "ar1_ft_plotlist_pos_phi", "ar1_ft_plotlist_neg_phi",
              "garch_ft_plotlist")) {
    if (!exists(f)) stop("Object" %s0% f %s0% "was not defined in" %s0% file)
  }

  ar_1_lrv <- function(sd, phi) {sd / (1 - phi)^2}
  garch_1_1_lrv <- function(alpha, beta, omega) {omega / (1 - alpha - beta)}

  plot_objs <- list(
    "ar1_plot_pos_phi" = ggplot(ar1_plotlist_pos_phi, aes(x = val)) +
      geom_density(fill = "red", alpha = .5) +
      facet_grid(phi ~ n, scales = "free_y") +
      geom_vline(aes(xintercept = ar_1_lrv(1, phi))) +
      xlab("Estimated Long Run Variance (theoretical as line)") +
      ylab("Density") +
      ggtitle("Estimated LRV for AR(1) processes
              depending on sample size (columns) and AR coefficient (rows)
              using Bartlett kernel and Andrews method for selecting bandwidth") +
      theme_bw(),
    "ar1_plot_neg_phi" = ggplot(ar1_plotlist_neg_phi, aes(x = val)) +
      geom_density(fill = "red", alpha = .5) +
      facet_grid(phi ~ n, scales = "free_y") +
      geom_vline(aes(xintercept = ar_1_lrv(1, phi))) +
      xlab("Estimated Long Run Variance (theoretical as line)") +
      ylab("Density") +
      xlim(0, 2) +
      #coord_cartesian(xlim = c(0, 2)) +
      ggtitle("Estimated LRV for AR(1) processes
              depending on sample size (columns) and AR coefficient (rows)
              using Bartlett kernel and Andrews method for selecting bandwidth") +
      theme_bw(),
    "garch_plot" = ggplot(garch_plotlist, aes(x = val)) +
      geom_density(fill = "red", alpha = .5) +
      facet_grid(. ~ n) +
      coord_cartesian(xlim = c(0, 2)) +
      geom_vline(aes(xintercept = garch_1_1_lrv(.1, .1, .5))) +
      xlab("Estimated Long Run Variance (theoretical as line)") +
      ylab("Density") +
      ggtitle("Estimated LRV for GARCH(1) processes
              depending on sample size (columns)
              using Bartlett kernel and Andrews method for selecting bandwidth") +
      theme_bw(),
    "ar1_ft_plot_pos_phi" = ggplot(ar1_ft_plotlist_pos_phi, aes(x = val)) +
      geom_density(fill = "red", alpha = .5) +
      facet_grid(phi ~ n, scales = "free_y") +
      geom_vline(aes(xintercept = ar_1_lrv(1, phi))) +
      xlab("Estimated Long Run Variance (theoretical as line)") +
      ylab("Density") +
      ggtitle("Estimated LRV for AR(1) processes
              depending on sample size (columns) and AR coefficient (rows)
              using flat-top kernel and Andrews method for selecting bandwidth") +
      theme_bw(),
    "ar1_ft_plot_neg_phi" = ggplot(fplotlist2, aes(x = val)) +
      geom_density(fill = "red", alpha = .5) +
      facet_grid(phi ~ n, scales = "free_y") +
      geom_vline(aes(xintercept = ar_1_lrv(1, phi))) +
      xlab("Estimated Long Run Variance (theoretical as line)") +
      ylab("Density") +
      xlim(0, 2) +
      #coord_cartesian(xlim = c(0, 2)) +
      ggtitle("Estimated LRV for AR(1) processes
              depending on sample size (columns) and AR coefficient (rows)
              using flat-top kernel and Andrews method for selecting bandwidth") +
      theme_bw(),
    "garch_ft_plot" = ggplot(fplotlist_garch, aes(x = val)) +
      geom_density(fill = "red", alpha = .5) +
      facet_grid(. ~ n) +
      coord_cartesian(xlim = c(0, 2)) +
      geom_vline(aes(xintercept = garch_1_1_lrv(.1, .1, .5))) +
      xlab("Estimated Long Run Variance (theoretical as line)") +
      ylab("Density") +
      ggtitle("Estimated LRV for GARCH(1) processes
              depending on sample size (columns)
              using flat-top kernel and Andrews method for selecting bandwidth") +
      theme_bw()
  )

  plot_suffix <- names(plot_objs) %s0% ifelse(pdf, ".tex", ".png")
  filenames <- prefix %s0% plot_suffix

  for (i in 1:length(plot_objs)) {
    if (pdf) {
      tikz(filenames[i], width = width, height = height, standAlone = TRUE)
    } else {
      png(filenames[i], width = width, height = height, units = "in")
    }
    print(plot_objs[[i]])
    dev.off()
    if (pdf) {texi2pdf(filenames[i])}
  }
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Plot the results of LRV estimation simulations.",
        option_list = list(
          make_option(c("--file", "-f"), type = "character",
                      help = "Input .Rda file"),
          make_option(c("--prefix", "-e"), type = "character",
                      default = "lrv_est_",
                      help = "Prefix for output files"),
          make_option(c("--pdf", "-p"), action = "store_true",
                      help = "Save output as pdf files"),
          make_option(c("--width", "-w"), type = "double", default = 5,
                      help = "Width of the resulting image"),
          make_option(c("--height", "-H"), type = "double", default = 7,
                      help = "Height of the resulting image")
        )
      ))

  do.call(main, cl_args)
}
