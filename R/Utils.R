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
