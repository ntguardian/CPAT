#!/usr/bin/Rscript
################################################################################
# SizeTableGenerator.R
################################################################################
# 2019-06-09
# Curtis Miller
################################################################################
# Generate a table for the "size" (Type I error rate) of simulated statistics
# when residuals are iid Normal in the presence of homo- or heteroskedasticity
# for differing sample sizes and for differing dimension.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

`%s%` <- CPAT:::`%s%`

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(basedir, output = "out.tex") {
  # This function will be executed when the script is called from the command
  # line

  library(dplyr)
  library(reshape2)
  library(xtable)

  setwd(basedir)

  files <- c("tb_norm_2d" = "Normal/SimsNormal.Rda",
             "tb_norm_hetero_2d" = "Normal/SimsNormalHetero.Rda",
             "tb_norm_4d" = "4DNormal/Sims4DNormal.Rda",
             "tb_norm_hetero_4d" = "4DNormal/Sims4DNormalHetero.Rda",
             "tb_norm_6d" = "6DNormal/Sims6DNormal.Rda",
             "tb_norm_hetero_6d" = "6DNormal/Sims6DNormalHetero.Rda")

  size_tables <- lapply(files, function(f) {
                                 load(f)
                
                                 tb <- power_sim_stat_data %>%
                                         filter(d == 0) %>%
                                         select(power, stat, n) %>%
                                         acast(stat ~ n,
                                               value.var = "power")
                                 tb <- tb[c('ZnRegNoKern', 'ZnResid',
                                            'VnResid', 'HS'), ]
                                 rownames(tb) <- c(
                                   "R\\'enyi-Type for Regression",
                                   "Univariate R\\'enyi-Type",
                                   "Univariate",
                                   "Hidalgo-Seo")

                                 tb
                               })

  size_tables_nonhetero <- size_tables[c("tb_norm_2d", "tb_norm_4d",
                                         "tb_norm_6d")]
  size_tables_hetero <- size_tables[c("tb_norm_hetero_2d", "tb_norm_hetero_4d",
                                      "tb_norm_hetero_6d")]

  arr <- array(c(unlist(size_tables_nonhetero), unlist(size_tables_hetero)),
               dim = c(nrow(size_tables_nonhetero[[1]]),
                       ncol(size_tables_nonhetero[[1]]),
                       length(size_tables_nonhetero), 2))
  dimnames(arr) <- list(rownames(size_tables_nonhetero[[1]]),
                        colnames(size_tables_nonhetero[[1]]),
                        c("2D", "4D", "6D"),
                        c("Homoskedastic", "Heteroskedastic"))

  flatarr <- ftable(arr, row.vars = c(4, 1), col.vars = c(3, 2))

  print(xtableFtable(flatarr, label = "tbl:statsize",
                     caption = paste("Simulated Type I error rates of test",
                                     "statistics depending on sample size,",
                                     "dimension, and whether the data is",
                                     "homoskedastic/heteroskedastic.")),
                                 file = output)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Generate a table for the \"size\" (Type I error rate) of" %s%
                  "simulated statistics when residuals are iid Normal in" %s%
                  "the presence of homo- or heteroskedasticity for" %s%
                  "differing sample sizes and for differing dimension.")
  p <- add_argument(p, "basedir", type = "character", nargs = 1,
                    help = "Base directory from which to search for" %s%
                           "properly formatted directories Normal/," %s%
                           "4DNormal/, and 6DNormal/")
  p <- add_argument(p, "--output", type = "character", default = "out.tex",
                    nargs = 1, help = "Output tex file containing table")

  cl_args <- parse_args(p)
  cl_args <- cl_args[!(names(cl_args) %in% c("help", "opts"))]
  if (any(sapply(cl_args, is.na))) {
    # User did not specify all inputs; print help message
    print(p)
    cat("\n\nNot all needed inputs were given.\n")
    quit()
  }

  do.call(main, cl_args[2:length(cl_args)])
}

