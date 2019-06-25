#!/usr/bin/Rscript
################################################################################
# ContextExample.R
################################################################################
# 2019-06-08
# Curtis Miller
################################################################################
# Context of simulations for structural change in regression models resembling
# context of data example.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# FUNCTION DEFINITIONS
################################################################################

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`
base_file_name <- CPAT:::base_file_name

#' Generate Matrix Containing Regime Change Information
#'
#' This function generates a \eqn{d \times 2} matrix where the first column is
#' the entry to \code{base} and the second column is the first plus a vector
#' that is \code{delta} times a length-one vector with all entries identical.
#'
#'
#' @param base The "base" model, pre-regime change
#' @param delta The amount by which to change the "base" model after the change
#' @return A \code{\link[base]{matrix}} with two columns and the same number of
#'         rows as the length of \code{base}
#' @examples
#' gen_regime_mat(1:2, 0.1)
gen_regime_mat <- function(base, delta = 0) {
  stopifnot(length(base) == 2)
  d <- 2
  delta_vec <- (c(-0.671, -0.722) - 
                c( 0    , -0.098)) * delta
  cbind(base, base + delta_vec)
}

################################################################################
# MAIN FUNCTION DEFINITION
################################################################################

main <- function(output = "ContextExample.Rda", help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  library(purrr)

  grm_md <- partial(gen_regime_mat, base = c(0, -0.098))
  delta <- ((0):20)/10

  ##############################################################################
  # REQUIRED OBJECTS
  ##############################################################################
  
  # n_values <- as.integer(c(216 + 9, 216 + 18, 216 + 36))
  n_values <- as.integer(c(250))
  kstar_functions <- c("ElectionDay" = function(n) {min(n, 216)})
  struc_models <- lapply(delta, grm_md)
  names(struc_models) <- "d" %s0% delta
  struc_name_conversion <- data.frame("d" = delta)
  rownames(struc_name_conversion) <- names(struc_models)

  save(n_values, kstar_functions, struc_models, struc_name_conversion,
       file = output, ascii = TRUE)
}

################################################################################
# INTERFACE SETUP
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = paste("Context of simulations for structural",
                            "change in regression models resembling example"),
        option_list = list(
          make_option(c("--output", "-o"), type = "character",
                      default = "ContextExample.Rda", help = "Output .Rda file")
        )
      ))

  do.call(main, cl_args)
}

