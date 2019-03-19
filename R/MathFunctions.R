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

#' Sequence \eqn{b_n} of the Darling-Erdös Law
#'
#' Computes \eqn{b_n(m) = (2 \log \log (n) + (m \log \log \log n) / 2 -
#' \log(\Gamma(m/2)))^2/(2 \log \log n)}
#'
#' @param n The parameter \eqn{n}
#' @param m The parameter \eqn{m}
#' @return The number \eqn{b_n(m)}
#' @examples
#' CPAT:::b_n(5, 2)
b_n <- function(n, m) {
  x <- (2 * log(log(n)) + (m * log(log(log(n))))/2 - log(gamma(m/2)))^2 /
    (2 * log(log(n)))
  x <- ifelse(!is.finite(x), 0, x)
  x
}

#' Sequence \eqn{a_n} of the Darling-Erdös Law
#'
#' Computes \eqn{a_n(m) = \sqrt{b_n(m)/(2 \log \log n)}}, with \eqn{b_n(m)} as
#' described by \code{\link{b_n}}.
#'
#' @param n The parameter \eqn{n}
#' @param m The parameter \eqn{m}
#' @return The number \eqn{a_n(m)}
#' @examples
#' CPAT:::a_n(5, 2)
a_n <- function(n, m) {
  x <- sqrt(b_n(n, m)/(2 * log(log(n))))
  x <- ifelse(!is.finite(x), 0, x)
  x
}

#' Bartlett Kernel
#'
#' Pure R implementation of the Bartlett kernel.
#'
#' \deqn{k(x) = \begin{cases} 1 - |x| & |x| \leq 1 \\ 0 & \text{o.w.}
#' \end{cases}}
#'
#' @param x Input to \eqn{k(x)}
#' @return Value of kernel
#' @examples
#' CPAT:::ba_kernel_r(0.5)
ba_kernel_r <- function(x) {
  x <- abs(x)
  ifelse(x <= 1, 1 - x, 0)
}

#' Truncated Kernel
#'
#' Pure R implementation of the truncated kernel.
#'
#' \deqn{k(x) = \begin{cases} 1 & |x| \leq 1 \\ 0 & \text{o.w.}
#' \end{cases}}
#'
#' @param x Input to \eqn{k(x)}
#' @return Value of kernel
#' @examples
#' CPAT:::tr_kernel_r(0.5)
tr_kernel_r <- function(x) {
  x <- abs(x)
  ifelse(x <= 1, 1, 0)
}

#' Parzen Kernel
#'
#' Pure R implementation of the Parzen kernel.
#'
#' \deqn{k(x) = \begin{cases} 1 - 6x^2 + 6|x|^3 & |x| \leq 1/2 \\ 
#' 2(1 - |x|)^3 & 1/2 < |x| \leq 1 \\ 0 & \text{o.w.}
#' \end{cases}}
#'
#' @param x Input to \eqn{k(x)}
#' @return Value of the kernel
#' @examples
#' CPAT:::pa_kernel_r(0)
#' CPAT:::pa_kernel_r(0.4)
#' CPAT:::pa_kernel_r(0.6)
#' CPAT:::pa_kernel_r(1.4)
pa_kernel_r <- function(x) {
  x <- abs(x)
  x <- ifelse(x <= 1/2, -x, x)
  x <- ifelse(x > 1, Inf, x)
  x <- ifelse(x <= 0, 1 - 6 * x^2 - 6 * x^3, -x)
  x <- ifelse(x <= 0, 2 * (1 + x)^3, x)
  x <- ifelse(is.infinite(x), 0, x)
  x
}

#' Tukey-Hanning Kernel
#'
#' Pure R implementation of the Tukey-Hanning kernel.
#'
#' \deqn{k(x) = \begin{cases} (1 + \cos(\pi x))/2 & |x| \leq 1 \\ 0 &
#' \text{o.w.} \end{cases}}
#'
#' @param x Input to \eqn{k(x)}
#' @return Value of the kernel
#' @examples
#' CPAT:::th_kernel_r(0.5)
th_kernel_r <- function(x) {
  x <- abs(x)
  ifelse(x <= 1, (1 + cos(pi * x))/2, 0)
}

#' Quadratic Spectral Kernel
#'
#' Pure R implementation of the quadratic spectral kernel.
#'
#' \deqn{k(x) = \frac{25}{12\pi^2 x^2} \left(\frac{\sin(6\pi x/5)}{6\pi x/5} -
#' \cos(6\pi x / 5)\right)}
#'
#' @param x Input to \eqn{k(x)}
#' @return Value of the kernel
#' @examples
#' CPAT:::qs_kernel_r(0.4)
qs_kernel_r <- function(x) {
  spxd5 <- 6 * pi * x / 5
  x <- ifelse(x == -Inf, Inf, x)
  x <- ifelse(x == 0, -Inf, x)
  x <- ifelse(is.finite(x), 3 * (sin(spxd5)/spxd5 - cos(spxd5))/(spxd5^2), x)
  x <- ifelse(x == Inf, 0, x)
  x <- ifelse(x == -Inf, 1, x)
  x
}

################################################################################
# MATH UTILITIES
################################################################################

#' Find Number of Summands Needed for Numerical Accuracy of \code{pBst}
#'
#' Find the number of summands needed to achieve numerical accuracy of the sum
#' involved in \code{\link{pBst}}.
#'
#' The number of summands needed is determined by using a loop that runs over
#' the summands until it encounters a summand that is not greater than the
#' specified level of numerical accuracy. The index of that last summand is then
#' returned.
#'
#' @param q Quantile input to CDF
#' @param b Point in space Bessel process hits
#' @param nu The parameter \eqn{\nu > -1} of the Bessel process
#' @param error The desired numerical error of the sum
#' @return Integer for number of summands
#' @export
#' @examples
#' pBst_summand_solver(1, 1)
pBst_summand_solver <- function(q, b, nu = -1/2,
                                error = .Machine$double.eps) {
  if ((q <= 0) || (b <= 0)) {
    return(1)
  }

  f <- function(k) {
    jnk <- besselJ_zeros(k, k, nu = nu)
    jnk^(nu - 1)/(besselJ(jnk, nu = nu + 1)) * exp(-jnk^2/(2 * b^2) * q) - error
  }

  k <- 1
  while (abs(f(k)) > error) {
    k <- k + 1
  }

  k
}

#' Find Number of Summands Needed for Numerical Accuracy of \code{dBst}
#'
#' Find the number of summands needed to achieve numerical accuracy of the sum
#' involved in \code{\link{dBst}}.
#'
#' The number of summands needed is determined by using a loop that runs over
#' the summands until it encounters a summand that is not greater than the
#' specified level of numerical accuracy. The index of that last summand is then
#' returned.
#'
#' @param x Quantile input to PDF
#' @param b Point in space Bessel process hits
#' @param nu The parameter \eqn{\nu > -1} of the Bessel process
#' @param error The desired numerical error of the sum
#' @return Integer for number of summands
#' @export
#' @examples
#' dBst_summand_solver(1, 1)
dBst_summand_solver <- function(x, b, nu = -1/2,
                                error = .Machine$double.eps) {
  if ((x <= 0) || (b <= 0)) {
    return(1)
  }

  f <- function(k) {
    jnk <- besselJ_zeros(k, k, nu = nu)
    jnk^(nu + 1)/(besselJ(jnk, nu = nu + 1)) * exp(-jnk^2/(2 * b^2) * x) - error
  }

  k <- 1
  while (abs(f(k)) > error) {
    k <- k + 1
  }

  k
}
