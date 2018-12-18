################################################################################
# MathFunctions.R
################################################################################
# 2018-12-17
# Curtis Miller
################################################################################
# Mathematical functions (that are not probability functions)
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
