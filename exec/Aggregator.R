#!/usr/bin/Rscript
################################################################################
# Aggregator.R
################################################################################
# 2019-02-23
# Curtis Miller
################################################################################
# Take simulation data and convert it to a plottable, tabular format.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# FUNCTIONS
################################################################################

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`
check_envir_has_objects <- CPAT:::check_envir_has_objects
stop_with_message <- CPAT:::stop_with_message

################################################################################
# MAIN FUNCTION DEFINITION
################################################################################

main <- function(input, output = NULL, TESTINPUT, CONTEXTINPUT, alpha = 0.05,
                 help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  suppressPackageStartupMessages(library(dplyr))

  test_tools <- new.env()
  input_objects <- new.env()
  context_tools <- new.env()
  load(TESTINPUT, envir = test_tools)
  load(input, envir = input_objects)
  load(CONTEXTINPUT, envir = context_tools)
  test_tools_expected_objects <- c("pval_functions", "plot_desc")
  input_objects_expected_objects <- c("res_list")
  context_tools_expected_objects <- c("struc_name_conversion")
  check_envir_has_objects(test_tools_expected_objects, envir = test_tools,
                          blame_string = TESTINPUT)
  check_envir_has_objects(input_objects_expected_objects, envir = input_objects,
                          blame_string = input)
  check_envir_has_objects(context_tools_expected_objects, envir = context_tools,
                          blame_string = CONTEXTINPUT)

  # Condition checking and object "importation"
  struc_name_conversion <- context_tools$struc_name_conversion
  stop_with_message(is.data.frame(struc_name_conversion) &
                    !any(is.null(rownames(struc_name_conversion))) &
                    all(sapply(struc_name_conversion, is.numeric)),
                    "Invalid struc_name_conversion from" %s% CONTEXTINPUT %s0%
                    "; must be a data.frame with numeric columns and" %s%
                    "all rows named, with names corresponding to names" %s0%
                    "struc_models")
  res_list <- input_objects$res_list
  stop_with_message(all(sapply(res_list, is.list)) &
                    all(sapply(res_list,
                               function(l) {all(sapply(l, is.list))})) &
                    all(sapply(res_list,
                               function(l) {all(sapply(l,
                                   function(m) {
                                     all(sapply(m, is.list))
                               }))})) &
                    all(sapply(res_list,
                               function(l) {all(sapply(l,
                                   function(m) {
                                     all(sapply(m, function(q) {
                                         all(sapply(q, is.numeric)) &
                                           !any(is.null(names(q))) &
                                           all(names(q) %in%
                                             rownames(struc_name_conversion))
                               }))}))})),
                    "Invalid res_list from" %s% input %s0% "; structure" %s%
                    "wrong (or perhaps it's struc_name_conversion from" %s%
                    CONTEXTINPUT %s0% ")")
  plot_desc <- test_tools$plot_desc
  stop_with_message(is.vector(plot_desc) & is.character(plot_desc) &
                    all(!is.null(names(plot_desc))),
                    "Invalid plot_desc from" %s% TESTINPUT %s0% "; must" %s%
                    "be a named character vector")
  pval_functions <- test_tools$pval_functions
  stop_with_message(is.vector(pval_functions) &
                    all(sapply(pval_functions,
                               function(f) {"q" %in% names(formals(f))})) &
                    all(!is.null(names(pval_functions))) &
                    all(sapply(res_list,
                               function(l) {
                                 all(sapply(l, function(m) {
                                       all(names(m) %in%
                                         names(pval_functions)) &
                                       all(names(m) %in% names(plot_desc))
                                     }))})),
                    "Invalid pval_functions from" %s% TESTINPUT %s0% ";" %s%
                    "must be a vector of functions that all take input 'q'" %s0%
                    ", and all test statistics in res_list must be present" %s%
                    "in the names of pval_functions")

  # Create long-form data frame
  res_list_collapsed <- lapply(res_list,
                               function(l) {
                                 # Now at n level
                                 lapply(l,
                                        function(m) {
                                          # Now at change point function level
                                          t <- lapply(names(m), function(q) {
                                            vals <- m[[q]]  # Change level
                                            pval_f <- Vectorize(
                                              pval_functions[[q]])
                                            pvals <- lapply(vals, pval_f)
                                            sig_prop <- sapply(pvals,
                                                               function(v) {
                                                                 mean(v <= alpha)
                                                               })
                                            names(sig_prop) <- names(vals)
                                            sig_prop
                                 })
                                 names(t) <- names(m)
                                 t
                               })})
  res_list_df <- lapply(res_list_collapsed,
                        function(l) {
                          # Now at n level
                          t <- lapply(names(l), function(m) {
                                        # Now at change point level
                                        cpt_list <- l[[m]]
                                        df <- stack(cpt_list)
                                        df$power <- df$values
                                        df$stat <- df$ind
                                        df$cpt <- m
                                        df$change <- rownames(df)
                                        rownames(df) <- NULL
                                        df[c("power", "stat", "cpt", "change")]
                                      })
                          names(t) <- names(l)
                          t
                        })
  res_list_df <- lapply(res_list_df,
                        function(l) {
                          names(l) <- NULL
                          l
                        })
  res_list_df <- lapply(res_list_df, function(l) {do.call(rbind, l)})
  res_list_df <- lapply(names(res_list_df),
                        function(l) {
                          df <- res_list_df[[l]]
                          df$n <- l
                          df
                        })
  merged_data_frame <- do.call(rbind, res_list_df)
  
  # Alter data in merged_data_frame
  merged_data_frame$n <- as.integer(substr(merged_data_frame$n, 2, 1000000000))
  merged_data_frame <- left_join(merged_data_frame,
                                 cbind("change" =
                                         rownames(struc_name_conversion),
                                       struc_name_conversion))
  power_sim_stat_data <- select(merged_data_frame, select = -change)
  power_sim_stat_data$stat <- as.character(power_sim_stat_data$stat)

  if (is.null(output)) {
    base <- CPAT:::base_file_name(input)
    output <- base %s0% "DataFrame.Rda"
  }

  save(power_sim_stat_data, plot_desc, file = output)
}

################################################################################
# INTERFACE SETUP
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Take simulation data and convert to tabular format.",
        option_list = list(
          make_option(c("--input", "-i"), type = "character",
                      help = "Input file with simulation data"),
          make_option(c("--output", "-o"), type = "character", default = NULL,
                      help = "Output .Rda file containing aggregated" %s%
                             "simulation data"),
          make_option(c("--alpha", "-a"), type = "double", default = 0.05,
                      help = "Statistical significance threshold"),
          make_option(c("--TESTINPUT", "-T"), type = "character",
                      help = "Name of .Rda file defining how to compute" %s%
                             "p-values for included statistical tests"),
          make_option(c("--CONTEXTINPUT", "-C"), type = "character",
                      help = "Name of .Rda file defining how to translate" %s%
                             "regime difference labels to numbers")
        )
      ))

  do.call(main, cl_args)
}

