#!/usr/bin/Rscript
################################################################################
# PowerSimStatDataGenerator.R
################################################################################
# 2018-09-04
# Curtis Miller
################################################################################
# Generates a data file containing results of power simulations
################################################################################


if (suppressPackageStartupMessages(!require("optparse"))) {
  install.packages("optparse")
  suppressPackageStartupMessages(optparse)
}

# Main function; executes the program
#
# See cl_args definition for description of arguments
main <- function(filemeta, statmeta, prefix = "", outfile = "out.csv",
  alpha = 0.05, help = FALSE) {
  power_sim_stat_df_creator <- CPAT:::power_sim_stat_df_creator
  perc_Zn_theo <- CPAT:::qZn(1 - alpha)
  perc_Vn_theo <- CPAT:::qkolmogorov(1 - alpha)
  perc_de_theo <- CPAT:::qdarling_erdos(1 - alpha)
  perc_hs_theo <- CPAT:::qhidalgo_seo(1 - alpha)

  power_sim_stat_data <- power_sim_stat_df_creator(filemeta, statmeta, prefix)
  write.csv(power_sim_stat_data, file = outfile, row.names = FALSE)
}

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
      description = paste("Generates a data file containing the results of",
                          "power simulations, stored in .Rda."),
      option_list = list(
          make_option(c("--filemeta", "-f"), type = "character",
                      help = paste("CSV file containing the names of files and",
                                   "which statistics they represent")),
          make_option(c("--statmeta", "-s"), type = "character",
                      help = "CSV file containing statistic metadata"),
          make_option(c("--prefix", "-p"), type = "character", default = "",
                      help = "Prefix to prepend to file names (say, location)"),
          make_option(c("--outfile", "-o"), type = "character",
                      default = "out.csv",
                      help = "Name of the destination CSV file"),
          make_option(c("--alpha", "-a"), type = "double", default = 0.05,
                      help = "Significance level of tests")
          )
      ))

  do.call(main, cl_args)
}
