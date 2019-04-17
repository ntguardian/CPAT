#!/usr/bin/Rscript
################################################################################
# TimeSeriespValPlotter.R
################################################################################
# 2019-04-17
# Curtis
################################################################################
# Plots p-values that are in sequential order, along with markers for
# significant events. Files are expected to be of the format produced by
# ExpandingWindowpValComputer.R
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`
check_envir_has_objects <- CPAT:::check_envir_has_objects
stop_with_message <- CPAT:::stop_with_message

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(input, output = "out", variable = "", level = 0.05,
                 levellinetype = "dotted", width = 6, height = 4) {
  # This function will be executed when the script is called from the command
  # line

  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(tikzDevice))
  suppressPackageStartupMessages(library(tools))

  input_env <- new.env()
  load(input, envir = input_env)

  # Assumption checking
  input_env_expected_objects <- c("events", "plot_desc", "stat_pvals")
  ckeck_envir_has_objects(input_env_expected_objects, envir = input_env,
                          blame_string = input)
  
  events <- input_env$events
  stop_with_message(is.data.frame(events) &
                    all(c("idx", "event") %in% names(events)) & 
                    class(events$idx) == numeric &
                    class(events$event) == "character",
                    "Invalid events from" %s% input %s0% "; must be data" %s%
                    "frame with character column 'event' and numeric column" %s%
                    "idx")

  plot_desc <- input_env$plot_desc
  stop_with_message(is.character(plot_desc) & all(!is.null(names(plot_desc))),
                    "Invalid plot_desc from" %s% input %s% "; must be" %s%
                    "named character vector")

  stat_pvals <- input_env$stat_pvals
  stop_with_message(is.data.frame(stat_pvals) &
                    length(stat_pvals) == (length(plot_desc) + 1) &
                    all(names(plot_desc) %in% names(stat_pvals)) &
                    all(sapply(stat_pvals[names(plot_desc)], is.numeric)),
                    "Invalid stat_pvals from" %s% input)

  if (variable == "") {
    variable = setdiff(names(stat_pvals), names(plot_desc))[[1]]
  }
  names(stat_pvals)[which(names(stat_pvals) == variable)] <- "Time"
  stop_with_message(is.numeric(stat_pvals$Time), "The variable" %s% variable %s%
                    "is interpreted as \"time\" and needs to be numeric")

  # Create data frame for plot
  # TODO: curtis: FINISH ME -- Wed 17 Apr 2019 05:02:06 PM MDT
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Plots p-values that are in sequential order, along with" %s%
                  "markers for significant events. Files are expected to" %s%
                  "be of the format produced by ExpandingWindowpValComputer.R")
  p <- add_argument(p, "input", type = "character", nargs = 1,
                    help = ".Rda file containing data to plot")
  p <- add_argument(p, "--output", type = "character", default = "out",
                    help = "Name of output file, without extension")
  p <- add_argument(p, "--variable", type = "character", default = "",
                    help = "The \"time\" variable (if not set," %s%
                           "automatically determined)")
  p <- add_argument(p, "--level", type = "double", default = 0.05,
                    help = "Where to plot a horizontal line representing" %s%
                           "the significance level of the test")
  p <- add_argument(p, "--levellinetype", type = "character",
                    default = "dotted",
                    help = "The type of line to plot for the level (set to" %s%
                           "\"blank\" for no line)")
  p <- add_argument(p, "--width", type = "double", default = 6,
                    help = "Width of each plot")
  p <- add_argument(p, "--height", type = "double", default = 4,
                    help = "Height of each plot")

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

