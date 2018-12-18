################################################################################
# MathFunctions.R
################################################################################
# 2018-12-17
# Curtis Miller
################################################################################
# Mathematical functions (that are not probability functions)
################################################################################

################################################################################
# MATHEMATICS
################################################################################

#' Compute Zeros of the Bessel Function of the First Kind
#'
#' Returns the zeros of the Bessel function of the first kind, \eqn{J_\nu}.
#'
#' This function is an interface to the function \code{besselJ_zeros_cpp}, a
#' function written in C++ and serves effectively as an interface to a Boost C++
#' function \code{cyl_bessel_j_zero}. Thus this function does nothing other than
#' make the Boost function available to R.
#'
#' See the references of \code{\link[base]{besselJ}} for more about bessel
#' functions.
#'
#' @param b The (one-based) index of the last zero to return
#' @param a The (one-based) index of the first zero to return (so \code{a = 1}
#'          represents the first positive zero)
#' @param nu The order of the Bessel function
#' @return A vector containing the zeros of the Bessel function
#' @examples
#' CPAT:::besselJ_zeros(4)
#' CPAT:::besselJ_zeros(a = 3, b = 10, nu = 3.5)
besselJ_zeros <- function(b, a = 1, nu = 1) {
  a <- floor(a)
  b <- floor(b)
  if (a < 1) {stop("Must have a >= 1")}
  if (b < a) {stop("Must have a <= b")}

  besselJ_zeros_cpp(nu, a, b)
}

################################################################################
# MATH UTILITIES
################################################################################

#' Find Number of Summands Needed for Numerical Accuracy of \code{pBst}
#'
#' Find the number of summands needed to achieve numerical accuracy of the sum
#' involved in \code{\link{pBst}}.
#'
#' The number of summands needed is found by using the root finder
#' \code{\link[stats]{uniroot}} to determine when one of the summands involed in
#' the infinite sum that partly describes the CDF computed by
#' \code{\link{pBst}} is equal to numerical accuracy.
#'
#' @param q Quantile input to CDF
#' @param b Point in space Bessel process hits
#' @param nu The parameter \eqn{\nu > -1} of the Bessel process
#' @param error The desired numerical error of the sum
#' @param lower Lower end point to search; will be the first entry of the
#'              \code{interval} vector passed to \code{\link[stats]{uniroot}}
#' @param upper Upper end point to search; will be the second entry of the
#'              \code{interval} vector passed to \code{\link[stats]{uniroot}}
#' @return Integer for number of summands
#' @examples
#' CPAT:::pBst_summand_solver(1, 1)
pBst_summand_solver <- function(q, b, nu = -1/2,
                                error = .Machine$double.eps,
                                lower = 1, upper = 1000) {
  if (upper <= lower) {stop("Must have lower < upper")}
  if ((q <= 0) || (b <= 0)) {
    return(1)
  }

  f <- function(k) {
    jnk <- besselJ_zeros(k, k, nu = nu)
    jnk^(nu - 1)/(besselJ(jnk, nu = nu + 1)) * exp(-jnk^2/(2 * b^2) * q) - error
  }

  ceiling(uniroot(f, c(lower, upper))$root)
}
