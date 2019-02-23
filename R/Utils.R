################################################################################
# Utils.R
################################################################################
# 2018-08-27
# Curtis Miller
################################################################################
# Functions and other extras used throughout the package
################################################################################

################################################################################
# ROXYGEN2 TAGS
################################################################################

#' @useDynLib CPAT
#' @importFrom Rcpp sourceCpp
NULL

################################################################################
# OPERATORS
################################################################################

#' Concatenate (With Space)
#'
#' Concatenate and form strings (with space separation)
#'
#' @param x One object
#' @param y Another object
#' @return A string combining \code{x} and \code{y} with a space separating them
#' @examples
#' `%s%` <- CPAT:::`%s%`
#' "Hello" %s% "world"
`%s%` <- function(x, y) {paste(x, y)}

#' Concatenate (Without Space)
#'
#' Concatenate and form strings (no space separation)
#'
#' @inheritParams %s%
#' @return A string combining \code{x} and \code{y}
#' @examples
#' `%s0%` <- CPAT:::`%s0%`
#' "Hello" %s0% "world"
`%s0%` <- function(x, y) {paste0(x, y)}

################################################################################
# CLASS CHECKING
################################################################################

#' Check For Formulas
#'
#' Checks if an object is a formula.
#'
#' @param x Object to check
#' @return \code{TRUE} if \code{x} is a \code{\link[stats]{formula}},
#'         \code{FALSE} otherwise
#' @examples
#' CPAT:::is.formula(y ~ x)
#' CPAT:::is.formula(2)
is.formula <- function(x) {
  is(x, "formula")
}

################################################################################
# ERROR CHECKING
################################################################################

#' Check An Environment for Objects
#'
#' Check that an environment has expected objects, and stop if it does not.
#'
#' @param objects A character vector listing what objects to expect in
#'                \code{envir}
#' @param envir The environment to check for objects
#' @param blame_string A string that gives more detailed output in error message
#'                     if not all files are found; default is the environment
#'                     passed to envir
#' @examples
#' x <- 1
#' CPAT:::check_envir_has_objects(c("x"))
check_envir_has_objects <- function(objects, envir = globalenv(),
                                    blame_string = NULL) {
  if (is.null(blame_string)) {
    blame_string <- deparse(substitute(envir))
  }
  if (!all(objects %in% ls(envir))) {
    stop(blame_string %s% "does not have all expected objects; must have" %s%
         paste(objects, collapse = ", "))
  }
}

#' Check For Condition and Stop With Message
#'
#' Check if \code{bool} is \code{TRUE}; if not, stop and report \code{message}
#'
#' @param bool Condition to check; if \code{FALSE}, \code{\link[base]{stop}} is
#'             called
#' @param message Message to report if \code{\link[base]{stop}} is called
#' @examples
#' x <- 1
#' CPAT:::stop_with_message(x == 1, message = "x is not 1")
stop_with_message <- function(bool, message = NULL) {
  if (is.null(message)) {
    message <- "all(" %s0% deparse(substitute(bool)) %s0% ") is FALSE"
  }

  if (!all(bool)) stop(message)
}

################################################################################
# SYSTEM UTILITIES
################################################################################

#' Extract Base File Name
#'
#' Extract the base name of the file without path or extension.
#'
#' @param x String from which to extract base name
#' @return A string containing the base file name without extension
#' @examples
#' CPAT:::base_file_name("~/Documents/test.txt")
base_file_name <- function(x) {
  tools::file_path_sans_ext(basename(x))
}

