#!/usr/bin/Rscript
################################################################################
# ContextSecondary.R
################################################################################
# 2019-09-15
# Curtis Miller
################################################################################
# Secondary context of simulations for structural change in regression models.
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
  d <- length(base)
  delta_vec <- delta * (10)^(0:(d - 1))

  cbind(base, base + delta_vec)
}

################################################################################
# MAIN FUNCTION DEFINITION
################################################################################

main <- function(output = "ContextSecondary.Rda", help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  library(purrr)

  grm_12 <- partial(gen_regime_mat, base = 1:2)
  delta <- ((-20):20)/10

  ##############################################################################
  # REQUIRED OBJECTS
  ##############################################################################
  
  n_values <- as.integer(c(50, 250, 500, 750, 1000))
  kstar_functions <- c("c35rd" = function(n) {n^{3/5}})
  struc_models <- lapply(delta, grm_12)
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
        description = paste("Primary context of simulations for structural",
                            "change in regression models"),
        option_list = list(
          make_option(c("--output", "-o"), type = "character",
                      default = "ContextSecondary.Rda", help = "Output .Rda file")
        )
      ))

  do.call(main, cl_args)
}

