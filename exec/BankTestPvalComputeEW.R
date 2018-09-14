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

  ff <- read.csv(ff_file, skip = 3, row.names = 1)
  portf <- read.csv(b_file, skip = 9,
                    stringsAsFactors = FALSE)
  portf <- portf[-nrow(portf),]
  portf$X <- as.character(portf$X)

  get_expanding_window_pvals_reg <- CPAT:::get_expanding_window_pvals_reg

  cl <- makeCluster(max(detectCores() - 1, min(2, detectCores())), outfile = "",
                    type = "FORK")
  registerDoParallel(cl)

  portf <- xts::xts(as.data.frame(lapply(portf[1:24099, -1], as.numeric)),
               order.by = as.Date(portf$X[1:24099], format = "%Y%m%d"))
  portf <- xts::xts(as.data.frame(lapply(portf, function(l) {
                 l[l == -99.99 | l == -999.99] <- NA
                 return(l)
               })), order.by = as.Date(rownames(as.data.frame(portf)),
                                       format = "%Y-%m-%d"))
  banks <- portf[,"Banks"]
  ff <- xts::xts(ff, order.by = as.Date(rownames(ff), format = "%Y%m%d"))

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
