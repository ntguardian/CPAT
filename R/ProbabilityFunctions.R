################################################################################
# ProbabilityFunctions.R
################################################################################
# 2018-08-10
# Curtis Miller
################################################################################
# Functions for probability distributions and simulating random variables
################################################################################

################################################################################
# DENSITY FUNCTIONS
################################################################################

#' Density Function of the First Hitting Time of a Bessel Process
#'
#' Density function of the distribution of the first time a Bessel process with
#' parameter \eqn{\nu > 1} hits \eqn{b > 0}.
#'
#' Let \eqn{\tau_b^{(\nu)}} be the first time a Bessel process with parameter
#' \eqn{\nu} hits \eqn{b > 0}. Let \eqn{J_{\nu}(x)} be the Bessel function
#' (of the first kind) with order \eqn{\nu}, and let \eqn{j_{\nu, k}} be the
#' \eqn{k}th zero of \eqn{J_{\nu}(x)}. Let \eqn{\Gamma(x)} be the gamma
#' function. Then the density function of \eqn{\tau_b^{(\nu)}} is
#'
#' \deqn{\frac{1}{2^\nu b^2 \Gamma(\nu + 1)} \sum_{k = 1}^{\infty} \frac{j_{\nu,
#' k}^{\nu + 1}}{J_{\nu + 1}(j_{\nu, k})} e^{-\frac{j_{\nu, k}^2}{2b^2}t}}
#'
#' This was found by differentiating the CDF computed by \code{\link{pBst}}.
#'
#' @param x Points at which to evaluate the density function
#' @param summands Number of summands to use in summation; default is to pick
#'                 the number of summands with
#'                 \code{\link{dBst_summand_solver}} (it could be slow, so for
#'                 performance it may be best to pick a fixed number)
#' @inheritParams pBst
#' @return The value of the density function at \code{x}
#' @examples
#' CPAT:::dBst(0.1, 1)
dBst <- function(x, b, nu = -1/2, summands = NULL) {
  if (is.null(summands)) {
    summands <- tryCatch({
      dBst_summand_solver(x = x, b = b, nu = nu)
    }, error = function(e) {
      warning("Solver to determine number of summands failed; defaulting to" %s%
              "500")
      500
    })
  }
  if (nu <= -1) {stop("Function valid only for nu > -1")}
  # Trivial cases
  if (b < 0) {
    return(0)
  } else if ((b == 0) & (x >= 0)) {
    return(0)
  }
  if (x < 0) {
    return(0)
  }

  jnk <- besselJ_zeros(summands, nu = nu)
  terms <- jnk^(nu + 1)/(besselJ(jnk, nu = nu + 1)) * exp(-jnk^2/(2 * b^2) * x)
  K <- (2^(nu) * gamma(nu + 1) * b^2)^(-1)

  return(K * sum(terms))
}
dBst <- Vectorize(dBst, "x")

#' Rényi-Type Statistic Limiting Distribution Density Function
#'
#' Function for computing the value of the density function of the limiting
#' distribution of the Rényi-type statistic.
#'
#' The density function was found by differentiating the CDF, as described by
#' \code{\link{pZn}}.
#'
#' @param x Point at which to evaluate the density function (note that this
#'          parameter is not vectorized)
#' @param d Dimension parameter
#' @param summands Number of summands to use in summation (the default should be
#'                 machine accurate)
#' @return Value of the density function at \eqn{x}
#' @examples
#' CPAT:::dZn(1)
dZn <- function(x, d = 1, summands = NULL) {
  nu <- d/2 - 1
  if (is.null(summands)) {
    summands <- tryCatch({
      dBst_summand_solver(x = 1, b = x, nu = nu)
    }, error = function(e) {
      warning("Solver to determine number of summands failed; defaulting to" %s%
              "500")
      500
    })
  }

  if (d < 1) {stop("Invalid input to d")}
  if (x <= 0) {
    return(0)
  }

  jnk <- besselJ_zeros(summands, nu = nu)
  terms <- jnk^(nu + 1)/(besselJ(jnk, nu = nu + 1)) * exp(-jnk^2/(2 * x^2))

  2 * (1 - pBst(q = 1, b = x, nu = nu, summands = summands)) * (2^(nu - 1) * 
    gamma(nu + 1) * x^3)^(-1) * sum(terms)
}
dZn <- Vectorize(dZn, "x")
  
################################################################################
# CDF FUNCTIONS
################################################################################

#' CDF of First Hitting Time of Bessel Process
#'
#' CDF of the distribution of the first time a Bessel process with parameter
#' \eqn{\nu > -1} hits \eqn{b > 0}.
#'
#' Let \eqn{\tau_b^{(\nu)}} be the first time a Bessel process with parameter
#' \eqn{\nu} hits \eqn{b > 0}. Let \eqn{J_{\nu}(x)} be the Bessel function
#' (of the first kind) with order \eqn{\nu}, and let \eqn{j_{\nu, k}} be the
#' \eqn{k}th zero of \eqn{J_{\nu}(x)}. Let \eqn{\Gamma(x)} be the gamma
#' function. Then the CDF of \eqn{\tau_b^{(\nu)}} is
#'
#' \deqn{1 - \frac{1}{2^{\nu - 1} \Gamma(\nu + 1)} \sum_{k = 1}^{\infty}
#' \frac{j_{\nu, k}^{\nu - 1}}{J_{\nu + 1}(j_{\nu, k})} e^{-\frac{j_{\nu,
#' k}^2}{2b^2} t}}
#'
#' (This was obtained in \insertCite{kent80}{CPAT}, but the formula above was
#' given in \insertCite{hamanamatsumoto13}{CPAT}.)
#'
#' @param q Quantile input to CDF
#' @param b Point in space Bessel process hits
#' @param nu The parameter \eqn{\nu > -1} of the Bessel process
#' @param summands Number of summands to use in summation; default is to pick
#'                 the number of summands with
#'                 \code{\link{pBst_summand_solver}} (it could be slow, so for
#'                 performance it may be best to pick a fixed number)
#' @return If \eqn{T} is the random variable as described, \eqn{P(T \leq q)}
#' @references
#'   \insertAllCited{}
#' @examples
#' CPAT:::pBst(1, 1)
pBst <- function(q, b, nu = -1/2, summands = NULL) {
  if (is.null(summands)) {
    summands <- tryCatch({
      pBst_summand_solver(q = q, b = b, nu = nu)
    }, error = function(e) {
      warning("Solver to determine number of summands failed; defaulting to" %s%
              "500")
      500
    })
  }
  if (nu <= -1) {stop("Function valid only for nu > -1")}
  # Trivial cases
  if (b < 0) {
    return(0)
  } else if ((b == 0) & (q >= 0)) {
    return(1)
  }
  if (q < 0) {
    return(0)
  }

  jnk <- besselJ_zeros(summands, nu = nu)
  terms <- jnk^(nu - 1)/(besselJ(jnk, nu = nu + 1)) * exp(-jnk^2/(2 * b^2) * q)
  K <- (2^(nu - 1) * gamma(nu + 1))^(-1)

  return(1 - K * sum(terms))
}
pBst <- Vectorize(pBst, "q")

#' Kolmogorov CDF
#'
#' CDF of the Kolmogorov distribution.
#'
#' @param q Quantile input to CDF
#' @param summands Number of summands for infinite sum (the default should have
#'        machine accuracy)
#' @return If \eqn{Z} is the random variable following the Kolmogorov
#'         distribution, the quantity \eqn{P(Z \leq q)}
#' @examples
#' CPAT:::pkolmogorov(0.1)
pkolmogorov <- function(q, summands = ceiling(q * sqrt(72) + 3/2)) {
  # Formerly called pKolmogorov
  sqrt(2 * pi) * sapply(q, function(x) { if (x > 0) {
    sum(exp(-(2 * (1:summands) - 1)^2 * pi^2/(8 * x^2)))/x
  } else {
    0
  }})
}
pkolmogorov <- Vectorize(pkolmogorov, "q")

#' Darling-Erdös Statistic CDF
#'
#' CDF for the limiting distribution of the Darling-Erdös statistic.
#'
#' @param q Quantile input to CDF
#' @return If \eqn{Z} is the random variable with this distribution, the
#'         quantity \eqn{P(Z \leq q)}
#' @examples
#' CPAT:::pdarling_erdos(0.1)
pdarling_erdos <- function(q) {
  # Formerly called pDarlingErdos
  exp(-2 * exp(-q))
}
pdarling_erdos <- Vectorize(pdarling_erdos, "q")

#' Rènyi-Type Statistic CDF
#'
#' CDF for the limiting distribution of the Rènyi-type statistic.
#'
#' If \eqn{G_{\nu, b}(x)} is the CDF of the first time a Bessel process with
#' parameter \eqn{\nu} hits \eqn{b > 0} (as described by \code{\link{pBst}})
#' then the CDF of the Rényi-type statistic when the null hypothesis is true is
#' \eqn{F(x) = (1 - G_{d/2 - 1, x}(1))^2}, where \eqn{d} is the dimensionality
#' parameter of the statistic. (This comes from combining the limiting
#' distribution of the statistic described in
#' \insertCite{horvathricemiller19}{CPAT} with the expression for the CDF of the
#' hitting time of the Bessel process described in
#' \insertCite{hamanamatsumoto13}{CPAT}.)
#'
#' @param q Quantile input to CDF
#' @param d Dimension parameter
#' @param summands Number of summands for infinite sum; if \code{NULL},
#'                 automatically determined using
#'                 \code{\link{pBst_summand_solver}} (which isn't necessarily
#'                 fast, so consider picking a fixed number if speed is
#'                 important)
#' @return If \eqn{Z} is the random variable following the limiting
#'         distribution, the quantity \eqn{P(Z \leq q)}
#' @references
#'   \insertAllCited{}
#' @examples
#' CPAT:::pZn(0.1)
pZn <- function(q, d = 1, summands = NULL) {
  if (d < 1) {stop("Invalid input to d")}
  if (q <= 0) {return(0)}
  (1 - pBst(q = 1, b = q, nu = d/2 - 1, summands = summands))^2
}
pZn <- Vectorize(pZn, "q")

#' Hidalgo-Seo Statistic CDF
#'
#' CDF of the limiting distribution of the Hidalgo-Seo statistic
#'
#' @param q Quantile input to CDF
#' @return If \eqn{Z} is the random variable following the limiting
#'         distribution, the quantity \eqn{P(Z \leq q)}
#' @examples
#' CPAT:::phidalgo_seo(0.1)
phidalgo_seo <- function(q) {
  # Formerly called pHidalgoSeo
  pdarling_erdos(q/2)
}

################################################################################
# QUANTILE FUNCTIONS
################################################################################

#' Bessel Process First Hitting Time Quantile Function
#'
#' Quantile function of the distribution of the first time a Bessel process with
#' parameter \eqn{\nu > -1} hits \eqn{b > 0}.
#'
#' This function uses \code{\link[stats]{uniroot}} for finding this quantity,
#' and many of the the accepted parameters are arguments for that function; see
#' its documentation for more details.
#'
#' @param p The probability associated with the desired quantile
#' @param interval,tol,... Arguments to be passed to
#'        \code{\link[stats]{uniroot}}
#' @inheritParams pBst
#' @return The quantile associated with \code{p}
#' @examples
#' CPAT:::qBst(0.5, b = 1)
qBst <- function(p, b, nu = -1/2, summands = NULL, interval = c(0, 100),
                 tol = .Machine$double.eps, ...) {
  if (p == 1) return(Inf)
  if (p == 0) return(0)
  if (p < 0 | p > 1) return(NaN)
  objective <- function(q) {pBst(q, b = b, nu = nu, summands = summands) - p}
  # Set up arguments for uniroot()
  args <- list(...); args$tol <- tol; args$interval <- interval
  args$f <- objective
  res <- do.call(uniroot, args)

  res$root
}
qBst <- Vectorize(qBst, "p")

#' Darling-Erdös Statistic Limiting Distribution Quantile Function
#'
#' Quantile function for the limiting distribution of the Darling-Erdös
#' statistic.
#'
#' @param p The probability associated with the desired quantile
#' @return The quantile associated with \code{p}
#' @examples
#' CPAT:::qdarling_erdos(0.5)
qdarling_erdos <- function(p) {
  # Formerly called qDarlingErdos
  -log(log(1/sqrt(p)))
}
qdarling_erdos <- Vectorize(qdarling_erdos, "p")

#' Hidalgo-Seo Statistic Limiting Distribution Quantile Function
#'
#' Quantile function for the limiting distribution of the Hidalgo-Seo statistic
#'
#' @param p The probability associated with the desired quantile
#' @return A The quantile associated with \code{p}
#' @examples
#' CPAT:::qhidalgo_seo(0.5)
qhidalgo_seo <- function(p) {
  # Formerly called qHidalgoSeo
  2 * qdarling_erdos(p)
}

#' Rènyi-Type Statistic Quantile Function
#'
#' Quantile function for the limiting distribution of the Rènyi-type statistic.
#'
#' This function uses \code{\link[stats]{uniroot}} for finding this quantity,
#' and many of the the accepted parameters are arguments for that function; see
#' its documentation for more details.
#'
#' @param p Value of the CDF at the quantile
#' @param d Dimension parameter
#' @param summands Number of summands for infinite sum
#' @param interval,tol,... Arguments to be passed to
#'        \code{\link[stats]{uniroot}}
#' @return The quantile associated with \code{p}
#' @examples
#' CPAT:::qZn(0.5)
qZn <- function(p, d = 1, summands = 500, interval = c(0, 100),
                tol = .Machine$double.eps, ...) {
  if (p == 1) return(Inf)
  if (p == 0) return(0)
  if (p < 0 | p > 1) return(NaN)
  objective <- function(q) {pZn(q, d = d, summands = summands) - p}
  # Set up arguments for uniroot()
  args <- list(...); args$tol <- tol; args$interval <- interval
  args$f <- objective
  res <- do.call(uniroot, args)

  res$root
}
qZn <- Vectorize(qZn, "p")

#' Kolmogorov Distribution Quantile Function
#'
#' Quantile function for the Kolmogorov distribution.
#'
#' This function uses \code{\link[stats]{uniroot}} for finding this quantity,
#' and many of the the accepted parameters are arguments for that function; see
#' its documentation for more details.
#'
#' @param p Value of the CDF at the quantile
#' @param summands Number of summands for infinite sum
#' @param interval,tol,... Arguments to be passed to 
#'        \code{\link[stats]{uniroot}}
#' @return The quantile associated with \code{p}
#' @examples
#' CPAT:::qkolmogorov(0.5)
qkolmogorov <- function(p, summands = 500, interval = c(0, 100),
                tol = .Machine$double.eps, ...) {
  if (p == 1) return(Inf)
  if (p == 0) return(0)
  if (p < 0 | p > 1) return(NaN)
  objective <- function(q) {pkolmogorov(q, summands = summands) - p}
  # Set up arguments for uniroot()
  args <- list(...); args$tol <- tol; args$interval <- interval
  args$f <- objective
  res <- do.call(uniroot, args)

  res$root
}
qkolmogorov <- Vectorize(qkolmogorov, "p")

################################################################################
# SIMULATION FUNCTIONS
################################################################################

#' Rènyi-Type Statistic Simulation (Assuming Variance)
#'
#' Simulates multiple realizations of the Rènyi-type statistic when the long-run
#' variance of the data is known.
#'
#' @param size Number of realizations to simulate
#' @param kn A function returning a positive integer that is used in the
#'           definition of the Rènyi-type statistic effectively setting the
#'           bounds over which the maximum is taken
#' @param n The sample size for each realization
#' @param gen_func The function generating the random sample from which the
#'                 statistic is computed
#' @param args A list of arguments to be passed to \code{gen_func}
#' @param sd The square root of the second moment of the data
#' @return A vector of simulated realizations of the Rènyi-type statistic
#' @examples
#' CPAT:::sim_Zn(100, kn = function(n) {floor(log(n))})
#' CPAT:::sim_Zn(100, kn = function(n) {floor(log(n))},
#'               gen_func = CPAT:::rchangepoint, args = list(changepoint = 250,
#'                                                           mean2 = 1))
sim_Zn <- function(size, kn, n = 500, gen_func = rnorm, args = NULL, sd = 1) {
  # Formerly called simZn
  Zn_realization <- function() {
    # Generate data set
    if (!is.list(args)) {
      dataset <- do.call(gen_func, list(n = n))
    } else {
      dataset <- do.call(gen_func, c(list(n = n), args))
    }

    max(sapply(floor(kn(n)):(n - floor(kn(n))), function(k)
          abs(1/k*sum(dataset[1:k]) -
            1/(n-k)*sum(dataset[(k+1):n]))
        ))
  }

  sqrt(kn(n))/sd * sapply(1:size, function(throwaway) Zn_realization())
}

#' CUSUM Statistic Simulation (Assuming Variance)
#'
#' Simulates multiple realizations of the CUSUM statistic when the long-run
#' variance of the data is known.
#'
#' @param size Number of realizations to simulate
#' @param n The sample size for each realization
#' @param gen_func The function generating the random sample from which the
#'                 statistic is computed
#' @param sd The square root of the second moment of the data
#' @param args A list of arguments to be passed to \code{gen_func}
#' @return A vector of simulated realizations of the CUSUM statistic
#' @examples
#' CPAT:::sim_Vn(100)
#' CPAT:::sim_Vn(100, gen_func = CPAT:::rchangepoint,
#'               args = list(changepoint = 250, mean2 = 1))
sim_Vn <- function(size, n = 500, gen_func = rnorm, sd = 1, args = NULL) {
  # Formerly called simVn
  Vn_realization <- function() {
    # Generate data set
    if (!is.list(args)) {
      dataset <- do.call(gen_func, list(n = n))
    } else {
      dataset <- do.call(gen_func, c(list(n = n),
                                     args))
    }

    dataset_mean <- mean(dataset)
    max(sapply(1:n, function(k) abs(sum(dataset[1:k]) - k * dataset_mean)))
  }

  (1/(sd*sqrt(n))) * sapply(1:size, function(throwaway) Vn_realization())
}

#' CUSUM Statistic Simulation
#'
#' Simulates multiple realizations of the CUSUM statistic.
#'
#' This differs from \code{sim_Vn()} in that the long-run variance is estimated
#' with this function, while \code{sim_Vn()} assumes the long-run variance is
#' known. Estimation can be done in a variety of ways. If \code{use_kernel_var}
#' is set to \code{TRUE}, long-run variance estimation using kernel-based
#' techniques will be employed; otherwise, a technique resembling standard
#' variance estimation will be employed. Any technique employed, though, will
#' account for the potential break points, as described in
#' \insertCite{horvathricemiller19;textual}{CPAT}. See the documentation for
#' \code{\link{stat_Vn}} for more details.
#' 
#' The parameters \code{kernel} and \code{bandwidth} control parameters for
#' long-run variance estimation using kernel methods. These parameters will be
#' passed directly to \code{\link{stat_Vn}}.
#'
#' Versions of the CUSUM statistic, such as the weighted or trimmed statistics,
#' can be simulated with the function by passing values to \code{kn} and
#' \code{tau}; again, see the documentation for \code{\link{stat_Vn}}.
#'
#' @param size Number of realizations to simulate
#' @param kn A function returning a positive integer that is used in the
#'           definition of the trimmed CUSUSM statistic effectively setting the
#'           bounds over which the maximum is taken
#' @param tau The weighting parameter for the weighted CUSUM statistic (defaults
#'            to zero for no weighting)
#' @param use_kernel_var Set to \code{TRUE} to use kernel-based long-run
#'                       variance estimation (\code{FALSE} means this is not
#'                       employed)
#' @param kernel If character, the identifier of the kernel function as used in
#'               the \pkg{cointReg} (see documentation for
#'               \code{cointReg::getLongRunVar}); if function, the kernel
#'               function to be used for long-run variance estimation (default
#'               is the Bartlett kernel in \pkg{cointReg}); this parameter
#'               has no effect if \code{use_kernel_var} is \code{FALSE}
#' @param bandwidth If character, the identifier of how to compute the bandwidth
#'                  as defined in the \pkg{cointReg} package (see
#'                  documentation for \code{cointReg::getLongRunVar}); if
#'                  function, a function to use for computing the bandwidth; if
#'                  numeric, the bandwidth to use (the default behavior is to
#'                  use the method described in \insertCite{andrews91b}{CPAT},
#'                  as used in \pkg{cointReg}); this parameter has no effect if
#'                  \code{use_kernel_var} is \code{FALSE}
#' @param n The sample size for each realization
#' @param gen_func The function generating the random sample from which the
#'                 statistic is computed
#' @param args A list of arguments to be passed to \code{gen_func}
#' @param parallel Whether to use the \pkg{foreach} and \pkg{doParallel}
#'                 packages to parallelize simulation (which needs to be
#'                 initialized in the global namespace before use)
#' @return A vector of simulated realizations of the CUSUM statistic
#' @references
#'  \insertAllCited{}
#' @examples
#' CPAT:::sim_Vn_stat(100)
#' CPAT:::sim_Vn_stat(100, kn = function(n) {floor(0.1 * n)}, tau = 1/3,
#'                    use_kernel_var = TRUE, gen_func = CPAT:::rchangepoint,
#'                    args = list(changepoint = 250, mean2 = 1))
sim_Vn_stat <- function(size, kn = function(n) {1}, tau = 0,
                        use_kernel_var = FALSE, kernel = "ba",
                        bandwidth = "and", n = 500, gen_func = rnorm,
                        args = NULL, parallel = FALSE) {
  # Formerly called simVnStat
  Vn_realization <- function() {
    # Generate data set
    if (!is.list(args)) {
      dataset <- do.call(gen_func, list(n = n))
    } else {
      dataset <- do.call(gen_func, c(list(n = n),
                                     args))
    }

    stat_Vn(dataset, kn = kn, tau = tau, use_kernel_var = use_kernel_var,
      kernel = kernel, bandwidth = bandwidth)
    # The following should be equivalent (yet slow) R code
    #n = length(dataset)
    #dataset_mean = mean(dataset)
    #return(max(sapply(
    #  1:n, function(k)
    #    abs(sum(dataset[1:k]) -
    #          k*dataset_mean)/sqrt((sum((dataset[1:k] -
    #           mean(dataset[1:k]))^2)+sum((dataset[(k+1):n] -
    #           mean(dataset[(k+1):n]))^2))/n)))/sqrt(n))
  }

  has_parallel <- requireNamespace("foreach", quietly = TRUE) &&
                  requireNamespace("doParallel", quietly = TRUE)
  if (parallel & !has_parallel) {
    warning("Either foreach or doParallel is not available; defaulting" %s%
            "to non-parallel implementation")
    parallel <- FALSE
  } else {
    times <- foreach::times
    `%dopar%` <- foreach::`%dopar%`
  }
  if (parallel) {
    foreach::foreach(i = 1:size, .combine = 'c') %dopar% Vn_realization()
  } else {
    sapply(1:size, function(throwaway) Vn_realization())
  }
}

#' Rènyi-Type Statistic Simulation
#'
#' Simulates multiple realizations of the Rènyi-type statistic.
#'
#' This differs from \code{sim_Zn()} in that the long-run variance is estimated
#' with this function, while \code{sim_Zn()} assumes the long-run variance is
#' known. Estimation can be done in a variety of ways. If \code{use_kernel_var}
#' is set to \code{TRUE}, long-run variance estimation using kernel-based
#' techniques will be employed; otherwise, a technique resembling standard
#' variance estimation will be employed. Any technique employed, though, will
#' account for the potential break points, as described in
#' \insertCite{horvathricemiller19;textual}{CPAT}. See the documentation for
#' \code{\link{stat_Zn}} for more details.
#' 
#' The parameters \code{kernel} and \code{bandwidth} control parameters for
#' long-run variance estimation using kernel methods. These parameters will be
#' passed directly to \code{\link{stat_Zn}}.
#'
#' @param size Number of realizations to simulate
#' @param kn A function returning a positive integer that is used in the
#'           definition of the Rènyi-type statistic effectively setting the
#'           bounds over which the maximum is taken
#' @param use_kernel_var Set to \code{TRUE} to use kernel-based long-run
#'                       variance estimation (\code{FALSE} means this is not
#'                       employed)
#' @param kernel If character, the identifier of the kernel function as used in
#'               the \pkg{cointReg} (see documentation for
#'               \code{cointReg::getLongRunVar}); if function, the kernel
#'               function to be used for long-run variance estimation (default
#'               is the Bartlett kernel in \pkg{cointReg}); this parameter
#'               has no effect if \code{use_kernel_var} is \code{FALSE}
#' @param bandwidth If character, the identifier of how to compute the bandwidth
#'                  as defined in the \pkg{cointReg} package (see
#'                  documentation for \code{cointReg::getLongRunVar}); if
#'                  function, a function to use for computing the bandwidth; if
#'                  numeric, the bandwidth to use (the default behavior is to
#'                  use the \insertCite{andrews91b;textual}{CPAT} method, as
#'                  used in \pkg{cointReg}); this parameter has no effect if
#'                  \code{use_kernel_var} is \code{FALSE}
#' @param n The sample size for each realization
#' @param gen_func The function generating the random sample from which the
#'                 statistic is computed
#' @param args A list of arguments to be passed to \code{gen_func}
#' @param parallel Whether to use the \pkg{foreach} and \pkg{doParallel}
#'                 packages to parallelize simulation (which needs to be
#'                 initialized in the global namespace before use)
#' @return A vector of simulated realizations of the Rènyi-type statistic
#' @references
#'  \insertAllCited{}
#' @examples
#' CPAT:::sim_Zn_stat(100)
#' CPAT:::sim_Zn_stat(100, kn = function(n) {floor(log(n))},
#'             use_kernel_var = TRUE, gen_func = CPAT:::rchangepoint,
#'             args = list(changepoint = 250, mean2 = 1))
sim_Zn_stat <- function(size, kn = function(n) {floor(sqrt(n))},
                        use_kernel_var = FALSE, kernel = "ba",
                        bandwidth = "and", n = 500, gen_func = rnorm,
                        args = NULL, parallel = FALSE) {
  # Formerly called simZnStat
  Zn_realization <- function() {
    # Generate data set
    if (!is.list(args)) {
      dataset <- do.call(gen_func, list(n = n))
    } else {
      dataset <- do.call(gen_func, c(list(n = n),
                                     args))
    }
    stat_Zn(dataset, kn, use_kernel_var = use_kernel_var,
      kernel = kernel, bandwidth = bandwidth)
  }

  has_parallel <- requireNamespace("foreach", quietly = TRUE) &&
                  requireNamespace("doParallel", quietly = TRUE)
  if (parallel & !has_parallel) {
    warning("Either foreach or doParallel is not available; defaulting" %s%
            "to non-parallel implementation")
    parallel <- FALSE
  } else {
    times <- foreach::times
    `%dopar%` <- foreach::`%dopar%`
  }
  if (parallel) {
    foreach::foreach(i = 1:size, .combine = 'c') %dopar% Zn_realization()
  } else {
    sapply(1:size, function(throwaway) Zn_realization())
  }
}

#' Darling-Erdös Statistic Simulation
#'
#' Simulates multiple realizations of the Darling-Erdös statistic.
#'
#' If \code{use_kernel_var} is set to \code{TRUE}, long-run variance estimation
#' using kernel-based techniques will be employed; otherwise, a technique
#' resembling standard variance estimation will be employed. Any technique
#' employed, though, will account for the potential break points, as described
#' in \insertCite{horvathricemiller19;textual}{CPAT}. See the documentation for
#' \code{\link{stat_de}} for more details.
#' 
#' The parameters \code{kernel} and \code{bandwidth} control parameters for
#' long-run variance estimation using kernel methods. These parameters will be
#' passed directly to \code{\link{stat_de}}.
#'
#' @param size Number of realizations to simulate
#' @param a The function that will be composed wit
#'          \eqn{l(x) = (2 \log(x))^{1/2}}
#' @param b The function that will be composed with
#'          \eqn{u(x) = 2 \log(x) + \frac{1}{2} \log(\log(x)) -
#'          \frac{1}{2}\log(pi)}
#' @param use_kernel_var Set to \code{TRUE} to use kernel-based long-run
#'                       variance estimation (\code{FALSE} means this is not
#'                       employed)
#' @param kernel If character, the identifier of the kernel function as used in
#'               the \pkg{cointReg} (see documentation for
#'               \code{cointReg::getLongRunVar}); if function, the kernel
#'               function to be used for long-run variance estimation (default
#'               is the Bartlett kernel in \pkg{cointReg}); this parameter
#'               has no effect if \code{use_kernel_var} is \code{FALSE}
#' @param bandwidth If character, the identifier of how to compute the bandwidth
#'                  as defined in the \pkg{cointReg} package (see
#'                  documentation for \code{cointReg::getLongRunVar}); if
#'                  function, a function to use for computing the bandwidth; if
#'                  numeric, the bandwidth to use (the default behavior is to
#'                  use the \insertCite{andrews91b;textual}{CPAT} method, as
#'                  used in \pkg{cointReg}); this parameter has no effect if
#'                  \code{use_kernel_var} is \code{FALSE}
#' @param n The sample size for each realization
#' @param gen_func The function generating the random sample from which the
#'                 statistic is computed
#' @param args A list of arguments to be passed to \code{gen_func}
#' @param parallel Whether to use the \pkg{foreach} and \pkg{doParallel}
#'                 packages to parallelize simulation (which needs to be
#'                 initialized in the global namespace before use)
#' @return A vector of simulated realizations of the Darling-Erdös statistic
#' @references
#'  \insertAllCited{}
#' @examples
#' CPAT:::sim_de_stat(100)
#' CPAT:::sim_de_stat(100, use_kernel_var = TRUE,
#'                    gen_func = CPAT:::rchangepoint,
#'                    args = list(changepoint = 250, mean2 = 1))
sim_de_stat <- function(size, a = log, b = log, use_kernel_var = FALSE,
                        kernel = "ba", bandwidth = "and", n = 500,
                        gen_func = rnorm, args = NULL, parallel = FALSE) {
  # Formerly called simDEStat
  de_realization <- function() {
    # Generate data set
    if (!is.list(args)) {
      dataset <- do.call(gen_func, list(n = n))
    } else {
      dataset <- do.call(gen_func, c(list(n = n), args))
    }

    stat_de(dataset, a = a, b = b, use_kernel_var = use_kernel_var,
      kernel = kernel, bandwidth = bandwidth)
  }

  has_parallel <- requireNamespace("foreach", quietly = TRUE) &&
                  requireNamespace("doParallel", quietly = TRUE)
  if (parallel & !has_parallel) {
    warning("Either foreach or doParallel is not available; defaulting" %s%
            "to non-parallel implementation")
    parallel <- FALSE
  } else {
    times <- foreach::times
    `%dopar%` <- foreach::`%dopar%`
  }
  if (parallel) {
    foreach::foreach(i = 1:size, .combine = 'c') %dopar% de_realization()
  } else {
    sapply(1:size, function(throwaway) de_realization())
  }
}

#' Hidalgo-Seo Statistic Simulation
#'
#' Simulates multiple realizations of the Hidalgo-Seo statistic.
#'
#' If \code{corr} is \code{TRUE}, then the residuals of the data-generating
#' process are assumed to be correlated and the test accounts for this in
#' long-run variance estimation; see the documentation for \code{\link{stat_hs}}
#' for more details. Otherwise, the sample variance is the estimate for the
#' long-run variance, as described in \insertCite{hidalgoseo13;textual}{CPAT}.
#'
#' @param size Number of realizations to simulate
#' @param corr Whether long-run variance should be computed under the assumption
#'             of correlated residuals
#' @param use_kernel_var Set to \code{TRUE} to use kernel-based long-run
#'                       variance estimation (\code{FALSE} means this is not
#'                       employed); \emph{TODO: NOT CURRENTLY IMPLEMENTED}
#' @param kernel If character, the identifier of the kernel function as used in
#'               the \pkg{cointReg} (see documentation for
#'               \code{cointReg::getLongRunVar}); if function, the kernel
#'               function to be used for long-run variance estimation (default
#'               is the Bartlett kernel in \pkg{cointReg}); this parameter
#'               has no effect if \code{use_kernel_var} is \code{FALSE};
#'               \emph{TODO: NOT CURRENTLY IMPLEMENTED}
#' @param bandwidth If character, the identifier of how to compute the bandwidth
#'                  as defined in the \pkg{cointReg} package (see
#'                  documentation for \code{cointReg::getLongRunVar}); if
#'                  function, a function to use for computing the bandwidth; if
#'                  numeric, the bandwidth to use (the default behavior is to
#'                  use the \insertCite{andrews91b;textual}{CPAT} method, as
#'                  used in \pkg{cointReg}); this parameter has no effect if
#'                  \code{use_kernel_var} is \code{FALSE}; \emph{TODO: NOT
#'                  CURRENTLY IMPLEMENTED}
#' @param n The sample size for each realization
#' @param gen_func The function generating the random sample from which the
#'                 statistic is computed
#' @param args A list of arguments to be passed to \code{gen_func}
#' @param parallel Whether to use the \pkg{foreach} and \pkg{doParallel}
#'                 packages to parallelize simulation (which needs to be
#'                 initialized in the global namespace before use)
#' @return A vector of simulated realizations of the Hidalgo-Seo statistic
#' @references
#'  \insertAllCited{}
#' @examples
#' CPAT:::sim_hs_stat(100)
#' CPAT:::sim_hs_stat(100, gen_func = CPAT:::rchangepoint, 
#'                    args = list(changepoint = 250, mean2 = 1))
sim_hs_stat <- function(size, corr = TRUE, gen_func = rnorm, args = NULL,
                        n = 500, parallel = FALSE, use_kernel_var = FALSE,
                        kernel = "ba", bandwidth = "and") {
  # Formerly known as simHSStat
  hs_realization <- function() {
    # Generate data set
    if (!is.list(args)) {
      dataset <- do.call(gen_func, list(n = n))
    } else {
      dataset <- do.call(gen_func, c(list(n = n),
                                     args))
    }

    stat_hs(dataset, corr = corr)
  }

  has_parallel <- requireNamespace("foreach", quietly = TRUE) &&
                  requireNamespace("doParallel", quietly = TRUE)
  if (parallel & !has_parallel) {
    warning("Either foreach or doParallel is not available; defaulting" %s%
            "to non-parallel implementation")
    parallel <- FALSE
  } else {
    times <- foreach::times
    `%dopar%` <- foreach::`%dopar%`
  }
  if (parallel) {
    foreach::foreach(i = 1:size, .combine = 'c') %dopar% hs_realization()
  } else {
    sapply(1:size, function(throwaway) hs_realization())
  }
}

#' Simulate Univariate Data With a Single Change Point
#'
#' This function simulates univariate data with a structural change.
#'
#' This function generates artificial change point data, where up to the
#' specified change point the data has one mean, and after the point it has a
#' different mean. By default, the function simulates standard Normal data with
#' no change. If \code{changepoint} is \code{NULL}, then by default the change
#' point will be at about the middle of the data.
#'
#' @param n An integer for the data set's sample size
#' @param changepoint An integer for where the change point occurs
#' @param mean1 The mean prior to the change point
#' @param mean2 The mean after the change point
#' @param dist The function with which random data will be generated
#' @param meanparam A string for the parameter in \code{dist} representing the
#'                  mean
#' @param ... Other arguments to be passed to dist
#' @return A vector of the simulated data
#' @examples
#' CPAT:::rchangepoint(500)
#' CPAT:::rchangepoint(500, changepoint = 10, mean2 = 2, sd = 2)
#' CPAT:::rchangepoint(500, changepoint = 250, dist = rexp, meanparam = "rate",
#'                     mean1 = 1, mean2 = 2)
rchangepoint <- function(n, changepoint = NULL, mean1 = 0, mean2 = 0,
                         dist = rnorm, meanparam = "mean", ...) {
  if (is.null(changepoint)) {
    changepoint <- ceiling(n/2)
  } else if (!is.integer(changepoint)) {
    changepoint <- round(changepoint)
  }

  # Arguments to be passed to dist prior to changepoint
  distargs_pre <- list(...)
  distargs_pre$n <- changepoint
  distargs_pre[[meanparam]] <- mean1

  # Arguments to be passed to dist after changepoint
  distargs_post <- list(...)
  distargs_post$n <- n - changepoint
  distargs_post[[meanparam]] <- mean2

  # Get and return data
  datavec1 <- do.call(dist, distargs_pre)
  datavec2 <- do.call(dist, distargs_post)

  c(datavec1, datavec2)
}

