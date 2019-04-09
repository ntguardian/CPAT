#!/usr/bin/Rscript
################################################################################
# CreateEnergyData.R
################################################################################
# 2019-04-08
# curtis
################################################################################
# Generate a data set containing energy data, energy_demand.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(price = "USEnergyPrices.csv",
                 consumption = "USEnergyConsumption.csv",
                 gdp = "GDP.csv", output = "energy_demand.rda") {
  # This function will be executed when the script is called from the command
  # line

  library(dplyr)
  library(reshape2)

  prices_raw <- read.csv(price)
  consumption_raw <- read.csv(consumption)
  gdp_raw <- read.csv(gdp)
  ptemp1 <- prices_raw$CPIENGSL %>% ts(start = 1957, end = 2019 + 1/12,
                                       frequency = 12)
  ctemp1 <- consumption_raw %>% select(Month, Primary.Energy.Consumption.Total)
  colnames(ctemp1) <- c("Month", "Consumption")
  ctemp1$Month <- as.character(ctemp1$Month)
  ctemp2 <- ctemp1$Consumption %>% ts(start = 1973, end = 2018 + 11/12,
                                      frequency = 12)
  gtemp1 <- gdp_raw$GDP %>% ts(start = 1947, end = 2018.75, frequency = 4)
  ptemp2 <- aggregate(ptemp1, nfrequency = 4, mean)
  ctemp3 <- aggregate(ctemp2, nfrequency = 4, sum)
  jtemp <- ts.intersect(gtemp1, ptemp2, ctemp3)
  colnames(jtemp) <- c("gdp", "price", "consumption")
  energy_demand <- jtemp
  save(energy_demand, file = output)
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Generate a data set containing energy data, energy_demand")
  p <- add_argument(p, "--price", type = "character",
                    default = "USEnergyPrices.csv",
                    help = "Data set containing energy prices")
  p <- add_argument(p, "--consumption", type = "character",
                    default = "USEnergyConsumption.csv",
                    help = "Data set containing energy consumption")
  p <- add_argument(p, "--gdp", type = "character",
                    default = "GDP.csv", help = "Data set containing U.S. GDP")
  p <- add_argument(p, "--output", type = "character",
                    default = "energy_demand.rda",
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

