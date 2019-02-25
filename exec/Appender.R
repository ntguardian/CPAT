#!/usr/bin/Rscript
################################################################################
# Appender.R
################################################################################
# 2019-02-25
# Curtis Miller
################################################################################
# Given files containing data frames created by the aggregator, combine all
# objects (data frames, vectors, etc.) into a unified file.
################################################################################

# argparser: A package for handling command line arguments
if (!suppessPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

################################################################################
# FUNCTIONS
################################################################################

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`
check_envir_has_objects <- CPAT:::check_envir_has_objects
stop_with_message <- CPAT:::stop_with_message

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(inputs, output = "Appended.Rda", help = FALSE, opts = NA) {
  # This function will be executed when the script is called from the command
  # line; the help and opts parameters does nothing, but are needed for
  # do.call() to work

  temp_env <- new.env()
  load(inputs[[1]], envir = temp_env)

  # Assumption checking
  temp_env_expected_objects <- c("power_sim_stat_data", "plot_desc")
  check_envir_has_objects(temp_env_expected_objects, envir = temp_env,
                          blame_string = inputs[[1]])
  power_sim_stat_data <- temp_env$power_sim_stat_data
  plot_desc <- temp_env$plot_desc
  stop_with_message(is.character(plot_desc) & all(!is.null(names(plot_desc))),
                    "Invalid plot_desc from" %s% inputs[[1]] %s0% "; must" %s%
                    "be named character vector")
  power_df_names <- names(power_sim_stat_data)
  stop_with_message(is.data.frame(power_sim_stat_data) &
                    all(!is.null(power_df_names)) &
                    all(c("power", "stat", "cpt", "n") %in% power_df_names) &
                    length(power_df_names) > 4 &  # 4 comes from length of above
                    all(sapply(
                        power_sim_stat_data[c("power", "stat", "cpt", "n")],
                        class) == c("numeric", "character", "character",
                                    "integer")) &
                    all(sapply(subset(
                        power_sim_stat_data, select = -c(power, stat, cpt, n)),
                        class) == "numeric") &
                    all(unique(power_sim_stat_data$stat) %in% names(plot_desc)),
                    "Invalid power_sim_stat_data from" %s% inputs[[1]])

  # Now for actual loop
  if (length(inputs) == 1) {
    save(power_sim_stat_data, plot_desc, file = output)
    return()
  }
  
  for (f in inputs[2:length(inputs)]) {
    load(f, envir = temp_env)
    check_envir_has_objects(temp_env_expected_objects, envir = temp_env,
                            blame_string = f)
    temp_df <- temp_env$power_sim_stat_data
    temp_desc <- temp_env$plot_desc
    stop_with_message(is.character(temp_desc) & all(!is.null(names(plot_desc))),
                      "Invalid plot_desc from" %s% f %s0% "; must be named" %s%
                      "character vector")
    stop_with_message(is.data.frame(temp_df) & 
                      length(temp_df) == length(power_sim_stat_data) &
                      all(!is.null(names(temp_df))) &
                      all(names(temp_df) %in% power_df_names) &
                      all(sapply(temp_df, class) ==
                        sapply(power_sim_stat_df, class)) &
                      all(unique(temp_df$stat) %in% names(temp_desc)),
                      "Invalid power_sim_stat_data from" %s% f)
  }
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Given files containing data frames and vectors created" %s%
                  "by the aggregator, combine all objects into a unified file.")
  p <- add_argument(p, "--inputs", type = "character", 
                    nargs = Inf, help = "Input files to combine")
  p <- add_argument(p, "--output", type = "character", default = "Appended.Rda",
                    help = "Output file for data")

  cl_args <- parse_args(p)

  do.call(main, cl_args[2:length(cl_args)])
}

