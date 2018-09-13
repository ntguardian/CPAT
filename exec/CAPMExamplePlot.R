#!/usr/bin/Rscript
################################################################################
# CAPMExamplePlot.R
################################################################################
# 2018-09-07
# Curtis Miller
################################################################################
# Create plots of the CAPM example comparing test statistics.
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

main <- function(file, outfile = "BankCAPMChange", verbose = FALSE,
                 kernelless = FALSE, width = 6, height = 4, help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work
  # See cl_args definition for parameter details

  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tikzDevice)
  library(reshape2)
  library(magrittr)

  load(file)

  if (!exists("test_p_vals")) stop("file did not contain test_p_vals")

  test_p_vals <- test_p_vals[148:211,]
  dat_df <- cbind(data.frame("i" = 1:length(test_p_vals[,1]),
                             "Date" = as.vector(rownames(test_p_vals))),
                  test_p_vals) %>%
      melt(id.vars = c("i", "Date"))

  dat_df$Date <- as.Date(dat_df$Date, "%Y-%m-%d")

  sig1_bk <- ggplot(dat_df %>%
                        filter(variable %in% c("cusum", "de", "hr", "andrew")),
                    aes(x = Date, y = ifelse(round(value, 15)==0,
                                             value_for_log0,
                                             -log10(value)), group = variable,
                        linetype = variable)) +
      geom_line(size = 0.85) +
      scale_linetype_manual(name = "Statistic",
                            labels = c("CUSUM", 'Darling-Erd\\"{o}s',
                                       "R\\'{e}nyi-type"),
                            values = c("longdash", "dotdash", "solid",
                                       "dotted")) +
      geom_hline(yintercept = -log10(0.05), color = "grey25",
                 linetype = "dashed") +
      geom_text(label = "Lehman Bros. Bankruptcy", x = as.Date("2008-09-15"),
                y = 3, angle = 90, vjust = -0.5, size = rel(3)) +
      geom_vline(xintercept = as.Date("2008-09-15"), linetype = "dotted") +
      scale_y_continuous(breaks = c(-log10(0.05), 5, 10, 15, 17),
                         labels = c("$-\\log_{10}(\\alpha)$", 5, 10, 15, 17)) +
      # xlab("End Year") + ylab("$-\\log_{10} p$") +
      xlab("") + ylab("") +
      theme_bw() +
      theme(plot.title = element_text(size = rel(2)),
            legend.position = "none")
            # legend.position = c(.15,.85))

  if (!kernelless) {
    sig2_bk <- sig1_bk %+% (dat_df %>%
                            filter(variable %in%
                                     c("cusum_ker", "de_ker", "hr_ker",
                                       "andrew")))
  }

  if (verbose) {cat("Creating" %s% outfile %s0% "\n")}
  tikz(outfile %s0% ".tex", width = 6, height = 4, standAlone = TRUE)
  if (kernelless) {
    print(sig1_bk)
  } else {
    print(sig2_bk)
  }
  dev.off()
  texi2pdf(outfile %s0% ".tex")
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Create plots of the CAPM example comparing test" %s%
                      "statistics.",
        option_list = list(
          make_option(c("--file", "-f"), type = "character",
                      help = ".Rda file containing test_p_vals object"),
          make_option(c("--outfile", "-o"), type = "character",
                      default = "BankCAPMChange",
                      help = "Stem of output plot file (can contain" %s%
                             "path information)"),
          make_option(c("--verbose", "-v"), action = "store_true",
                      help = "Show extra information"),
          make_option(c("--kernelless", "-k"), action = "store_true",
                      help = "Do not use kernel-based long-run variance" %s%
                             "estimation techniques"),
          make_option(c("--width", "-w"), type = "double", default = 6,
                      help = "Width of the plot"),
          make_option(c("--height", "-H"), type = "double", default = 4,
                      help = "Height of the plot")
        )
      ))

  do.call(main, cl_args)
}

