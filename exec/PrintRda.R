#!/usr/bin/Rscript
################################################################################
# PrintRda.R
################################################################################
# 2019-03-04
# Curtis Miller
################################################################################
# Print the contents of a .Rda file to the screen in a readable format.
################################################################################

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

################################################################################
# EXECUTABLE SCRIPT MAIN FUNCTIONALITY
################################################################################

main <- function(file) {
  # This function will be executed when the script is called from the command
  # line

  file_env <- new.env()
  load(file, envir = file_env)
  for (obj in ls(file_env)) {
    cat(paste0("\n\n", paste0(rep("-", nchar(obj)), collapse = ""), "\n", obj, "\n",
        paste0(rep("-", nchar(obj)), collapse = ""), "\n\n"))
    str(file_env[[obj]])
  }
}

################################################################################
# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION
################################################################################

if (sys.nframe() == 0) {
  p <- arg_parser("Print the contents of a .Rda file in a readable format.")
  p <- add_argument(p, "file", type = "character", nargs = 1,
                    help = "The file to load")

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

