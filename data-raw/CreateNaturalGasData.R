#!/usr/bin/Rscript
################################################################################
# CreateNaturalGasData.R
################################################################################
# 2019-04-10
# Curtis Miller
################################################################################
# Generate a data set containing natural gas price, demand, and GDP.
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

main <- function(gas = "EIANaturalGasConsumptionPrice.csv",
                 gdp = "UnadjustedGDP.csv", output = "natural_gas_demand.rda") {
  # This function will be executed when the script is called from the command
  # line

  library(dplyr)

  gas_raw <- read.csv(gas)
  gdp_raw <- read.csv(gdp)
  
  ptemp1 <- gas_raw %>% filter(!is.na(Price)) %>% .$Price %>% ts(start = 1976,
                                                                 end = 2012 +
                                                                       11/12,
                                                                 frequency = 12)
  ctemp1 <- gas_raw$Consumption %>% ts(start = 1973, end = 2018 + 11/12,
                                       frequency = 12)
  gtemp1 <- gdp_raw$NA000334Q %>% ts(start = 1947, end = 2018 + 3/4,
                                     frequency = 4)

  ptemp2 <- aggregate(ptemp1, nfrequency = 4, mean)
  ctemp2 <- aggregate(ctemp1, nfrequency = 4, sum)
  jtemp <- ts.intersect(gtemp1, ptemp2, ctemp2)
  colnames(jtemp) <- c("gdp", "price", "consumption")
  natural_gas_demand <- jtemp
  save(natural_gas_demand, file = output)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Generate a data set containing natural gas price," %s%
                  "demand, GDP")
  p <- add_argument(p, "--gas", type = "character",
                    default = "EIANaturalGasConsumptionPrice.csv",
                    help = "Data set containing energy prices and demand")
  p <- add_argument(p, "--gdp", type = "character",
                    default = "UnadjustedGDP.csv",
                    help = "Data set containing U.S. GDP")
  p <- add_argument(p, "--output", type = "character",
                    default = "natural_gas_demand.rda",
                    help = "Name of output .rda file")

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

