#!/usr/bin/Rscript
################################################################################
# UnivariatePlotter.R
################################################################################
# 2019-02-25
# Curtis Miller
################################################################################
# Plot power simulation data when the changing variable is univariate, and save
# the plots.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
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

main <- function(input, prefix = "", variable = "", level = 0.05,
                 levellinetype = "dotted", width = 3, height = 2,
                 verbose = FALSE, help = FALSE, opts = NA) {
  # This function will be executed when the script is called from the command
  # line; the help and opts parameters does nothing, but are needed for
  # do.call() to work

  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(tikzDevice))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tools))

  input_env <- new.env()
  load(input, envir = input_env)

  # Assumption checking
  input_env_expected_objects <- c("power_sim_stat_data", "plot_desc")
  check_envir_has_objects(input_env_expected_objects, envir = input_env,
                          blame_string = input)
  power_sim_stat_data <- input_env$power_sim_stat_data
  plot_desc <- input_env$plot_desc
  stop_with_message(is.character(plot_desc) & all(!is.null(names(plot_desc))),
                    "Invalid plot_desc from" %s% input %s% "; must be" %s%
                    "named character vector")
  power_df_names <- names(power_sim_stat_data)
  if (variable == "") {
    variable <- power_df_names[which(!(power_df_names %in%
                               c("power", "stat", "cpt", "n")))[1]]
  }
  stop_with_message(is.data.frame(power_sim_stat_data) &
                    all(!is.null(power_df_names)) &
                    all(c("power", "stat", "cpt", "n") %in% power_df_names) &
                    length(power_df_names) > 4 &
                    variable %in% power_df_names &
                    all(sapply(
                        power_sim_stat_data[c("power", "stat", "cpt", "n",
                                              variable)],
                        class) == c("numeric", "character", "character",
                                    "integer", "numeric")) &
                    all(unique(power_sim_stat_data$stat) %in% names(plot_desc)),
                    "Invalid power_sim_stat_data from" %s% input)

  # Create data frame for plot
  plot_df <- power_sim_stat_data[c("power", "stat", "cpt", "n", variable)]
  names(plot_df)[5] <- "var"

  # Create basis plot
  n_start <- unique(plot_df$n)[1]  # This is just initialization
  p_orig <- ggplot(data = plot_df %>% filter(n == n_start),
                     aes(x = var, y = power, group = stat, linetype = stat)) +
                   geom_line() +
                   scale_linetype_manual(values = plot_desc) +
                   ylim(0, 1) +
                   theme_bw() +
                   geom_hline(yintercept = level, linetype = levellinetype) +
                   theme(legend.position = "none",
                         plot.title = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank())

  # Actual loop
  n_levels <- unique(plot_df$n)
  filenames_base <- prefix %s0% "n" %s0% n_levels
  names(filenames_base) <- n_levels
  filenames_tex <- filenames_base %s0% ".tex"
  filenames_pdf <- filenames_base %s0% ".pdf"
  names(filenames_tex) <- names(filenames_base)
  names(filenames_pdf) <- names(filenames_base)
  cwd <- getwd()  # Get the current working directory
  for (nc in n_levels) {
    nch <- as.character(nc)
    df <- plot_df %>% filter(n == nc)
    if (verbose) {
      cat("Creating plot for n =", nc, "\nCorresponding files:",
          filenames_tex[nch], filenames_pdf[nch], "\n")
    }
    p <- p_orig %+% df
    tikz(filenames_tex[nch], width = width, height = height, standAlone = TRUE,
         sanitize = TRUE)
    print(p)
    dev.off()
    setwd(dirname(filenames_tex[[nch]]))
    texi2pdf(basename(filenames_tex[nch]), clean = TRUE, quiet = verbose)
    setwd(cwd)
  }

  if (verbose) {
    cat("Loop complete; files created:\n")
    print(data.frame("n" = n_levels, "TeX File" = filenames_tex,
                     "PDF File" = filenames_pdf))
    cat("\nStatistic Lines Key:\n")
    print(plot_desc)
    cat("\n\nLevel Line:", levellinetype, "\n\n")
  }
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Plot power simulation data when the changing variable" %s%
                  "is univariate, and save the plots.")
  p <- add_argument(p, "input", type = "character", nargs = 1,
                    help = ".Rda file containing data to plot")
  p <- add_argument(p, "--prefix", type = "character", default = "",
                    help = "Prefix for output file names (including directory)")
  p <- add_argument(p, "--variable", type = "character", default = "",
                    help = "The x-axis variable (if not set, automatically" %s%
                           "determined)")
  p <- add_argument(p, "--level", type = "double", default = 0.05,
                    help = "Where to plot a horizontal line representing" %s%
                           "the significance level of the test")
  p <- add_argument(p, "--levellinetype", type = "character",
                    default = "dotted",
                    help = "The type of line to plot for the level (set to" %s%
                           "\"blank\" for no line)")
  p <- add_argument(p, "--width", type = "double", default = 3,
                    help = "Width of each plot")
  p <- add_argument(p, "--height", type = "double", default = 2,
                    help = "Height of each plot")
  p <- add_argument(p, "--verbose", flag = TRUE, help = "Extra messages")

  cl_args <- parse_args(p)

  do.call(main, cl_args[2:length(cl_args)])
}

