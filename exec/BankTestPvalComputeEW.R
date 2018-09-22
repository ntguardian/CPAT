#! /usr/bin/Rscript
################################################################################
# BankTestPvalComputeEW.R
################################################################################
# 2018-09-07
# Curtis Miller
################################################################################
# Generate p-values for various statistical tests applied to bank portfolio
# data, when the CAPM model's parameters are computed
################################################################################

if (suppressPackageStartupMessages(!require("optparse"))) {
  install.packages("optparse")
  suppressPackageStartupMessages(optparse)
}

# Main function; executes the program
#
# See cl_args definition for description
main <- function(ff_file, b_file, out = "BankCAPMPValues.Rda", help = FALSE) {
  library(doParallel)
  library(foreach)
  library(xts)
  library(CPAT)

  data(ff)
  data(banks)

  ff <- xts::xts(ff, order.by = as.Date(rownames(ff), format = "%Y%m%d"))
  banks <- xts::xts(banks, order.by = as.Date(rownames(banks),
                                              format = "%Y-%m-%d"))

  model_banks_df <- merge(banks["2005/2008"], ff, join = "inner")[-1,]
  names(model_banks_df)[1] <- "Return"
  tail(model_banks_df)

  capm_model <- lm(I(Return - RF) ~ Mkt.RF + HML + RMW + SMB + CMA,
                   data = model_banks_df)

  test_p_vals <- get_expanding_window_pvals_reg(I(Return - RF) ~ Mkt.RF + HML + 
                                                  RMW + SMB + CMA, 
                                                data = as.data.frame(
                                                  model_banks_df),
                                                min_n = 754, m = 931)

  rownames(test_p_vals) <- rownames(as.data.frame(
    model_banks_df[test_p_vals[, "n"],]
  ))

  save(test_p_vals, file = out)
}

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
      description = paste("Applied the CUSUM, Darling-Erdos, Hidalgo-Seo,",
                          "Andrews, and Renyi-type tests to data for the",
                          "Fama-French model applied to a portfolio of returns",
                          "of banking sector stocks, and creates a file",
                          "containing the results of the tests"),
      option_list = list(
        make_option(c("--ff-file", "-f"), type = "character",
            help = "The location of the CSV file containing Fama-French data"),
        make_option(c("--bank-file", "-b"), type = "character",
            help = "The location of the CSV file containing stock data"),
        make_option(c("--out", "-o"), type = "character",
            default = "BankCAPMPValues.Rda",
            help = "The name of the output file")
      )
  ))

  names(cl_args) <- c("ff_file", "b_file", "out", "help")
  do.call(main, cl_args)
}
