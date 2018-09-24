################################################################################
# ChangePointTests.R
################################################################################
# 2018-08-10
# Curtis Miller
################################################################################
# Definition of all change point statistical tests.
################################################################################

################################################################################
# GENERAL STATISTICAL FUNCTIONS
################################################################################

#' Variance Estimation Consistent Under Change
#'
#' Estimate the variance (using the sum of squared errors) with an estimator
#' that is consistent when the mean changes at a known point.
#'
#' This is the estimator
#'
#' \deqn{\hat{\sigma}^2_{T,t} = T^{-1}\left(\sum_{s = 1}^t \left(X_s -
#'       \bar{X}_t\right)^2 + \sum_{s = t + 1}^{T}\left(X_s - \tilde{X}_{T - t}
#'       \right)^2\right)}
#'
#' where \eqn{\bar{X}_t = t^{-1}\sum_{s = 1}^t X_s} and \eqn{\tilde{X}_{T - t} =
#' (T - t)^{-1} \sum_{s = t + 1}^{T} X_s}. In this implementation, \eqn{T} is
#' computed automatically as \code{length(x)} and \code{k} corresponds to
#' \eqn{t}, a potential change point.
#'
#' @param x A numeric vector for the data set
#' @param k The potential change point at which the data set is split
#' @import stats
#' @return The estimated change-consistent variance
#' @examples
#' CPAT:::cpt_consistent_var(c(rnorm(500, mean = 0), rnorm(500, mean = 1)), k = 500)
cpt_consistent_var <- function(x, k) {
  n <- length(x)
  if (n < k | k < 0) {stop("k must be an integer between 1 and length(x)")}

  x1 <- x[1:k]
  x2 <- x[(k + 1):n]
  mu1 <- mean(x1)
  mu2 <- mean(x2)
  sse1 <- sum((x1 - mu1)^2)
  sse2 <- sum((x2 - mu2)^2)

  (sse1 + sse2)/n
}

#' Weights for Long-Run Variance
#'
#' Compute some weights for long-run variance. This code comes directly from the
#' source code of \pkg{cointReg}; see \code{\link[cointReg]{getLongRunWeights}}.
#'
#' @param n Length of weights' vector
#' @param bandwidth A number for the bandwidth
#' @param kernel The kernel function; see \code{\link[cointReg]{getLongRunVar}}
#'        for possible values
#' @return List with components \code{w} containing the vector of weights and
#'         \code{upper}, the index of the largest non-zero entry in \code{w}
#' @examples
#' CPAT:::getLongRunWeights(10, 1)
getLongRunWeights <- function(n, bandwidth, kernel = "ba") {
  w <- numeric(n - 1)
  bw<- bandwidth
  if (kernel == "tr") {
    w <- w + 1
    upper <- min(bw, n - 1)
  }
  else if (kernel == "ba") {
    upper <- ceiling(bw) - 1
    if (upper > 0) {
        j <- 1:upper
    }
    else {
        j <- 1
    }
    w[j] <- 1 - j/bw
  }
  else if (kernel == "pa") {
    upper1 <- floor(bw/2)
    if (upper1 > 0) {
        j <- 1:upper1
    }
    else {
        j <- 1
    }
    jj <- j/bw
    w[j] <- 1 - 6 * jj^2 + 6 * jj^3
    j2 <- (floor(bw/2) + 1):bw
    jj2 <- j2/bw
    w[j2] <- 2 * (1 - jj2)^3
    upper <- ceiling(bw) - 1
  }
  else if (kernel == "bo") {
    upper <- ceiling(bw) - 1
    if (upper > 0) {
        j <- 1:upper
    }
    else {
        j <- 1
    }
    jj <- j/bw
    w[j] <- (1 - jj) * cos(pi * jj) + sin(pi * jj)/pi
  }
  else if (kernel == "da") {
    upper <- n - 1
    j <- 1:upper
    w[j] <- sin(pi * j/bw)/(pi * j/bw)
  }
  else if (kernel == "qs") {
    sc <- 1.2 * pi
    upper <- n - 1
    j <- 1:upper
    jj <- j/bw
    w[j] <- 25/(12 * pi^2 * jj^2) * (sin(sc * jj)/(sc * jj) - 
        cos(sc * jj))
  }
  if (upper <= 0) 
    upper <- 1
  list(w = w, upper = upper)
}

#' Long-Run Variance Estimation With Possible Change Points
#'
#' Computes the estimates of the long-run variance in a change point context, as
#' described in \insertCite{horvathricemiller19}{CPAT}. By default it uses
#' kernel and bandwidth selection as used in the package \pkg{cointReg}, though
#' changing the parameters \code{kernel} and \code{bandwidth} can change this
#' behavior. If \pkg{cointReg} is not installed, the Bartlett internal (defined
#' internally) will be used and the bandwidth will be the square root of the
#' sample size.
#'
#' @param dat The data vector
#' @param kernel If character, the identifier of the kernel function as used in
#'               \pkg{cointReg} (see \code{\link[cointReg]{getLongRunVar}}); if
#'               function, the kernel function to be used for long-run variance
#'               estimation (default is the Bartlett kernel in \pkg{cointReg})
#' @param bandwidth If character, the identifier for how to compute the
#'               bandwidth as defined in \pkg{cointReg} (see
#'               \code{\link[cointReg]{getBandwidth}}); if function, a function
#'               to use for computing the bandwidth; if numeric, the bandwidth
#'               value to use (the default is to use Andrews' method, as used in
#'               \pkg{cointReg})
#' @return A vector of estimates of the long-run variance
#' @references
#'  \insertAllCited{}
#' @examples
#' x <- rnorm(1000)
#' CPAT:::get_lrv_vec(x)
#' CPAT:::get_lrv_vec(x, kernel = "pa", bandwidth = "nw")
get_lrv_vec <- function(dat, kernel = "ba", bandwidth = "and") {
  has_cointreg <- requireNamespace("cointReg", quietly = TRUE)

  n <- length(dat)

  if (is.character(bandwidth)) {
    if (!has_cointreg) {
      warning("cointReg is not installed! Defaulting to sqrt.")
      bandwidth <- sqrt
    } else {
      if (is.character(kernel) &&
          kernel %in%  c("ba", "pa", "qs", "th", "tr")) {
         kervar <- kernel
      } else {
         kervar = "ba"
      }
      h <- cointReg::getBandwidth(dat, bandwidth = bandwidth, kernel = kervar)
    }
  } else if (is.numeric(bandwidth)) {
    if (bandwidth <= 0) {
      stop("Bandwidth must be greater than zero.")
    } else {
      h <- bandwidth
    }
  } else if (!is.function(bandwidth)) {
    stop(paste("bandwidth must be a function, a valid character string, or a",
               "non-negative number."))
  } else {
    h <- bandwidth(n)
  }

  if (is.character(kernel)) {
    kern_vals <- c(1, getLongRunWeights(n, kernel = kernel,
                                        bandwidth = h)$w[-n])
  } else if (!is.function(kernel)) {
    stop("kernel must be a function or a valid character string.")
  } else {
    kern_vals <- sapply(1:(n - 1)/h, kernel)
  }

  # The maximum lag that needs to be checked
  max_l <- max(which(kern_vals != 0))

  # dat <- rnorm(100)
  #dat <- 1:10

  dat_mat <- matrix(rep(dat, times = n), byrow = FALSE, nrow = n)
  dat_lower <- dat_mat; dat_upper <- dat_mat;
  dat_lower[lower.tri(dat_mat, diag = FALSE)] <- NA
  dat_upper[!lower.tri(dat_mat, diag = FALSE)] <- NA

  x_bar <- colMeans(dat_lower, na.rm = TRUE)
  x_tilda <- colMeans(dat_upper, na.rm = TRUE)

  mean_mat <- matrix(rep(x_bar, times = n), byrow = TRUE, nrow = n)
  mean_mat[lower.tri(diag(n), diag = FALSE)] <- matrix(rep(x_tilda, times = n),
                                                       byrow = TRUE,
                                                       nrow = n)[lower.tri(
                                                         diag(n), diag = FALSE)]

  # l increases by row (starts at 0), t by column
  y <- dat_mat - mean_mat

  # Function implemented in Rcpp for speed
  # It accepts the matrix y and the evaluation of the kernel function stored
  # in kern_vals, and returns a vector containing estimated long-run variances
  # at points t
  sigma <- get_lrv_vec_cpp(y, kern_vals, max_l)

  # "Equivalent" R code (slower)
  # covs <- sapply(1:n, function(t) {
  #     sapply(0:(n - 1), function(l) {
  #         mean(y[1:(n - l),t]*y[(1 + l):n,t])
  #     })
  # })
  #
  # sigma <- sapply(2:(n - 2), function(t) {
  #     covs[0 + 1, t] + 2 * sum(vkernel(1:(n - 1)/h) * covs[1:(n - 1), t])
  # })

  if (any(sigma < 0)) {
    warning(paste("A negative variance was computed! This may be due to a bad",
                  "kernel being chosen."))
  }

  sigma
}

################################################################################
# TEST STATISTIC COMPUTATION FUNCTIONS
################################################################################

#' Compute the CUSUM Statistic
#'
#' This function computes the CUSUM statistic (and can compute weighted/trimmed
#' variants, depending on the values of \code{kn} and \code{tau}).
#'
#' The definition of the statistic is
#'
#' \deqn{T^{-1/2} \max_{1 \leq t \leq T} \hat{\sigma}_{t,T}^{-1} \left|
#'       \sum_{s = 1}^t X_s - \frac{t}{T}\sum_{s = 1}^T \right|}
#'
#' A more general version is
#'
#' \deqn{T^{-1/2} \max_{t_T \leq t \leq T - t_T} \hat{\sigma}_{t,T}^{-1}
#'       \left(\frac{t}{T}
#'       \left(\frac{T - t}{T}\right)\right)^{\tau} \left| \sum_{s = 1}^t X_s -
#'       \frac{t}{T}\sum_{s = 1}^T \right|}
#'
#' The parameter \code{kn} corresponds to the trimming parameter \eqn{t_T} and
#' the parameter \code{tau} corresponds to \eqn{\tau}.
#'
#' See \insertCite{horvathricemiller19}{CPAT} for more details.
#'
#' @param dat The data vector
#' @param kn A function corresponding to the trimming parameter \eqn{t_T} in the
#'           trimmed CUSUM variant; by default, is a function returning 1 (for
#'           no trimming)
#' @param tau The weighting parameter \eqn{\tau} for the weighted CUSUM
#'            statistic; by default, is 0 (for no weighting)
#' @param estimate Set to \code{TRUE} to return the estimated location of the
#'                 change point
#' @param use_kernel_var Set to \code{TRUE} to use kernel methods for long-run
#'                       variance estimation (typically used when the data is
#'                       believed to be correlated); if \code{FALSE}, then the
#'                       long-run variance is estimated using
#'                       \eqn{\hat{\sigma}^2_{T,t} = T^{-1}\left(
#'                       \sum_{s = 1}^t \left(X_s - \bar{X}_t\right)^2 +
#'                       \sum_{s = t + 1}^{T}\left(X_s -
#'                       \tilde{X}_{T - t}\right)^2\right)}, where
#'                       \eqn{\bar{X}_t = t^{-1}\sum_{s = 1}^t X_s} and
#'                       \eqn{\tilde{X}_{T - t} = (T - t)^{-1}
#'                       \sum_{s = t + 1}^{T} X_s}
#' @param kernel If character, the identifier of the kernel function as used in
#'               \pkg{cointReg} (see \code{\link[cointReg]{getLongRunVar}}); if
#'               function, the kernel function to be used for long-run variance
#'               estimation (default is the Bartlett kernel in \pkg{cointReg})
#' @param bandwidth If character, the identifier for how to compute the
#'               bandwidth as defined in \pkg{cointReg} (see
#'               \code{\link[cointReg]{getBandwidth}}); if function, a function
#'               to use for computing the bandwidth; if numeric, the bandwidth
#'               value to use (the default is to use Andrews' method, as used in
#'               \pkg{cointReg})
#' @param get_all_vals If \code{TRUE}, return all values for the statistic at
#'                     every tested point in the data set
#' @return If both \code{estimate} and \code{get_all_vals} are \code{FALSE}, the
#'         value of the test statistic; otherwise, a list that contains the test
#'         statistic and the other values requested (if both are \code{TRUE},
#'         the test statistic is in the first position and the estimated change
#'         point in the second)
#' @references
#'  \insertAllCited{}
#' @examples
#' CPAT:::stat_Vn(rnorm(1000))
#' CPAT:::stat_Vn(rnorm(1000), kn = function(n) {0.1 * n}, tau = 1/2)
#' CPAT:::stat_Vn(rnorm(1000), use_kernel_var = TRUE, bandwidth = "nw", kernel = "bo")
stat_Vn <- function(dat, kn = function(n) {1}, tau = 0, estimate = FALSE,
                   use_kernel_var = FALSE, kernel = "ba", bandwidth = "and",
                   get_all_vals = FALSE) {
  # Formerly named statVn()

  # Here is equivalent (slow) R code
  # n = length(dat)
  # return(n^(-1/2)*max(sapply(
  #     floor(max(kn(n),1)):min(n - floor(kn(n)),n-1), function(k)
  #         (k/n*(1 - k/n))^(-tau)/
  #         sqrt((sum((dat[1:k] -
  #            mean(dat[1:k]))^2)+sum((dat[(k+1):n] -
  #            mean(dat[(k+1):n]))^2))/n)*
  #         abs(sum(dat[1:k]) - k/n*sum(dat))
  # )))

  if (use_kernel_var) {
    lrv_est <- get_lrv_vec(dat, kernel = kernel, bandwidth = bandwidth)
  } else {
    # A vector must be passed to stat_Vn_cpp, so we'll pass one with an
    # unlikely value
    lrv_est <- c(-1)  
  }
  res <- stat_Vn_cpp(dat, kn(length(dat)), tau, use_kernel_var, lrv_est, 
                    get_all_vals)
  res[[2]] <- as.integer(res[[2]])
  if (!estimate & !get_all_vals) {
    return(res[[1]])
  } else {
    return(res[c(TRUE, estimate, get_all_vals)])
  }
}

#' Compute the Darling-Erdös Statistic
#'
#' This function computes the Darling-Erdös statistic.
#'
#' If \eqn{\bar{A}_T(\tau, t_T)} is the weighted and trimmed CUSUM statistic
#' with weighting parameter \eqn{\tau} and trimming parameter \eqn{t_T} (see
#' \code{\link{stat_Vn}}), then the Darling-Erdös statistic is
#'
#' \deqn{l(a_T) \bar{A}_T(1/2, 1) - u(b_T)}
#'
#' with \eqn{l(x) = \sqrt{2 \log x}} and \eqn{u(x) = 2 \log x + \frac{1}{2} \log
#' \log x - \frac{1}{2} \log \pi} (\eqn{\log x} is the natural logarithm of
#' \eqn{x}). The parameter \code{a} corresponds to \eqn{a_T} and \code{b} to
#' \eqn{b_T}; these are both \code{log} by default.
#'
#' See \insertCite{horvathricemiller19}{CPAT} to learn more.
#'
#' @param dat The data vector
#' @param a The function that will be composed with
#'          \eqn{l(x) = (2 \log x)^{1/2}}
#' @param b The function that will be composed with
#'          \eqn{u(x) = 2 \log x + \frac{1}{2} \log \log x - \frac{1}{2} \log
#'          \pi}
#' @param estimate Set to \code{TRUE} to return the estimated location of the
#'                 change point
#' @param use_kernel_var Set to \code{TRUE} to use kernel methods for long-run
#'                       variance estimation (typically used when the data is
#'                       believed to be correlated); if \code{FALSE}, then the
#'                       long-run variance is estimated using
#'                       \eqn{\hat{\sigma}^2_{T,t} = T^{-1}\left(
#'                       \sum_{s = 1}^t \left(X_s - \bar{X}_t\right)^2 +
#'                       \sum_{s = t + 1}^{T}\left(X_s -
#'                       \tilde{X}_{T - t}\right)^2\right)}, where
#'                       \eqn{\bar{X}_t = t^{-1}\sum_{s = 1}^t X_s} and
#'                       \eqn{\tilde{X}_{T - t} = (T - t)^{-1}
#'                       \sum_{s = t + 1}^{T} X_s}
#' @param kernel If character, the identifier of the kernel function as used in
#'               \pkg{cointReg} (see \code{\link[cointReg]{getLongRunVar}}); if
#'               function, the kernel function to be used for long-run variance
#'               estimation (default is the Bartlett kernel in \pkg{cointReg})
#' @param bandwidth If character, the identifier for how to compute the
#'               bandwidth as defined in \pkg{cointReg} (see
#'               \code{\link[cointReg]{getBandwidth}}); if function, a function
#'               to use for computing the bandwidth; if numeric, the bandwidth
#'               value to use (the default is to use Andrews' method, as used in
#'               \pkg{cointReg})
#' @param get_all_vals If \code{TRUE}, return all values for the statistic at
#'                     every tested point in the data set
#' @return If both \code{estimate} and \code{get_all_vals} are \code{FALSE}, the
#'         value of the test statistic; otherwise, a list that contains the test
#'         statistic and the other values requested (if both are \code{TRUE},
#'         the test statistic is in the first position and the estimated change
#'         point in the second)
#' @references
#'  \insertAllCited{}
#' @examples
#' CPAT:::stat_de(rnorm(1000))
#' CPAT:::stat_de(rnorm(1000), use_kernel_var = TRUE, bandwidth = "nw", kernel = "bo")
stat_de <- function(dat, a = log, b = log, estimate = FALSE,
                    use_kernel_var = FALSE, kernel = "ba", bandwidth = "and",
                    get_all_vals = FALSE) {

  # Formerly known as statDE()
  n <- length(dat)
  l <- function(x) {sqrt(2*log(x))}
  u <- function(x) {2*log(x) + 1/2*log(log(x)) - 1/2*log(pi)}

  res <- stat_Vn(dat, kn = function(n) {1}, tau = 1/2, estimate = TRUE,
                 use_kernel_var = use_kernel_var, kernel = kernel,
                 bandwidth = bandwidth, get_all_vals = get_all_vals)
  res[[2]] <- as.integer(res[[2]])
  res[[1]] <- l(a(n)) + res[[1]] - u(b(n))
  if (get_all_vals) {
    res[[3]] <- l(a(n)) + res[[3]] - u(b(n))
  }
  if (!estimate & !get_all_vals) {
    return(res[[1]])
  } else {
    return(res[c(TRUE, estimate, get_all_vals)])
  }
}

#' Compute the Hidalgo-Seo Statistic
#'
#' This function computes the Hidalgo-Seo statistic for a change in mean model.
#'
#' For a data set \eqn{x_t} with \eqn{n} observations, the test statistic is
#'
#' \deqn{\max_{1 \leq s \leq n - 1} 	(\mathcal{LM}(s) - B_n)/A_n}
#'
#' where \eqn{\hat{u}_t = x_t - \bar{x}} (\eqn{\bar{x}} is the sample mean),
#' \eqn{a_n = (2 \log \log n)^{1/2}}, \eqn{b_n = a_n^2 - \frac{1}{2} \log \log
#' \log n - \log \Gamma (1/2)}, \eqn{A_n = b_n / a_n^2}, \eqn{B_n =
#' b_n^2/a_n^2}, \eqn{\hat{\Delta} = \hat{\sigma}^2 = n^{-1} \sum_{t = 1}^{n}
#' \hat{u}_t^2}, and \eqn{\mathcal{LM}(s) = n (n - s)^{-1} s^{-1}
#' \hat{\Delta}^{-1} \left( \sum_{t = 1}^{s} \hat{u}_t\right)^2}.
#'
#' If \code{corr} is \code{FALSE}, then the residuals are assumed to be
#' uncorrelated. Otherwise, the residuals are assumed to be correlated and
#' \eqn{\hat{\Delta} = \hat{\gamma}(0) + 2 \sum_{j = 1}^{\lfloor \sqrt{n}
#' \rfloor} (1 - \frac{j}{\sqrt{n}}) \hat{\gamma}(j)} with \eqn{\hat{\gamma}(j)
#' = \frac{1}{n}\sum_{t = 1}^{n - j} \hat{u}_t \hat{u}_{t + j}}.
#'
#' This statistic was presented in \insertCite{hidalgoseo13}{CPAT}.
#'
#' @param dat The data vector
#' @param estimate Set to \code{TRUE} to return the estimated location of the
#'                 change point
#' @param corr If \code{TRUE}, the long-run variance will be computed under the
#'             assumption of correlated residuals; ignored if \code{custom_var}
#'             is not \code{NULL} or \code{use_kernel_var} is \code{TRUE}
#' @param get_all_vals If \code{TRUE}, return all values for the statistic at
#'                   every tested point in the data set
#' @param custom_var Can be a vector the same length as \code{dat} consisting of
#'                   variance-like numbers at each potential change point (so
#'                   each entry of the vector would be the "best estimate" of
#'                   the long-run variance if that location were where the
#'                   change point occured) or a function taking two parameters
#'                   \code{x} and \code{k} that can be used to generate this
#'                   vector, with \code{x} representing the data vector and
#'                   \code{k} the position of a potential change point; if
#'                   \code{NULL}, this argument is ignored
#' @param use_kernel_var Set to \code{TRUE} to use kernel methods for long-run
#'                       variance estimation (typically used when the data is
#'                       believed to be correlated); if \code{FALSE}, then the
#'                       long-run variance is estimated using
#'                       \eqn{\hat{sigma}^2_{T,t} = T^{-1}\left(
#'                       \sum_{s = 1}^t \left(X_s - \bar{X}_t\right)^2 +
#'                       \sum_{s = t + 1}^{T}\left(X_s -
#'                       \tilde{X}_{T - t}\right)^2\right)}, where
#'                       \eqn{\bar{X}_t = t^{-1}\sum_{s = 1}^t X_s} and
#'                       \eqn{\tilde{X}_{T - t} = (T - t)^{-1}
#'                       \sum_{s = t + 1}^{T} X_s}; if \code{custom_var} is not
#'                       \code{NULL}, this argument is ignored
#' @param kernel If character, the identifier of the kernel function as used in
#'               \pkg{cointReg} (see \code{\link[cointReg]{getLongRunVar}}); if
#'               function, the kernel function to be used for long-run variance
#'               estimation (default is the Bartlett kernel in \pkg{cointReg})
#' @param bandwidth If character, the identifier for how to compute the
#'               bandwidth as defined in \pkg{cointReg} (see
#'               \code{\link[cointReg]{getBandwidth}}); if function, a function
#'               to use for computing the bandwidth; if numeric, the bandwidth
#'               value to use (the default is to use Andrews' method, as used in
#'               \pkg{cointReg})
#' @return If both \code{estimate} and \code{get_all_vals} are \code{FALSE}, the
#'         value of the test statistic; otherwise, a list that contains the test
#'         statistic and the other values requested (if both are \code{TRUE},
#'         the test statistic is in the first position and the estimated change
#'         point in the second)
#' @references
#'  \insertAllCited{}
#' @examples
#' CPAT:::stat_hs(rnorm(1000))
#' CPAT:::stat_hs(rnorm(1000), corr = FALSE)
stat_hs <- function(dat, estimate = FALSE, corr = TRUE, get_all_vals = FALSE,
                    custom_var = NULL, use_kernel_var = FALSE, kernel = "ba",
                    bandwidth = "and") {
  # Formerly named statHS()

  n <- length(dat)
  mu <- mean(dat)
  u <- dat - mu

  if (corr) {
    m <- sqrt(n)
    Delta <- sum(u^2) / n + sum(sapply(1:m, function(j) {
      2 * (1 - j/m) * sum(u[(j + 1):n] * u[1:(n - j)])/n
    }))
  } else {
    # Delta will be the analog of lrv in stat_Zn, but this variable was already
    # defined and the naming is in alignment with Hidalgo's paper
    Delta <- ((n-1)/n) * var(u)
  }

  # To implement kernel LRV estimation that is consistent under H_A, I play
  # around with custom_var; this will be the avenue through which that
  # functionality is implemented.
  if (is.null(custom_var)) {
    if (use_kernel_var) {
      Delta <- get_lrv_vec(dat, kernel, bandwidth)
    } else {
      Delta <- rep(Delta, length(dat))
    }
  } else {
    if (is.function(custom_var)) {
      # To allow greater flexibility, Delta will be made a vector so that Delta
      # can take values at different possible change points; this allows for
      # greater experimentation on our part

      # This may seem silly, but this is so that error codes refer to
      # custom_var, and we don't want recursion either
      custom_var_temp <- custom_var
      custom_var <- purrr::partial(custom_var_temp, x = dat, .lazy = FALSE)
      custom_var_vec <- vapply(1:length(dat), custom_var,
                               FUN.VALUE = numeric(1))
    } else if (is.numeric(custom_var)) {
      if (length(custom_var) < length(dat)) stop("custom_var must have" %s%
                                                 "length at least" %s%
                                                 length(dat) %s0% ", the" %s%
                                                 "length of the data set")
      custom_var_vec <- custom_var
    } else {
      stop("Don't know how to handle custom_var of class" %s% class(custom_var))
    }
    if (any(custom_var_vec < 0, na.rm = TRUE)) {
      stop("custom_var suggests a negative variance, which is impossible.")
    }
    Delta <- custom_var_vec
  }

  la_mu <- sapply(1:(n - 1), function(s) {
    n/(n - s) * 1/s * sum(u[1:s])^2/Delta[s]
  })

  a_n <- sqrt(2 * log(log(n)))
  b_n <- 2 * log(log(n)) - log(log(log(n)))/2 - log(pi/4)/2

  A_n <- b_n/a_n^2  # Matching notation in Hidalgo and Seo (2013)
  B_n <- A_n * b_n

  stat <- (max(la_mu) - B_n)/A_n
  est <- which.max(la_mu)
  res <- list("statistic" = stat,
              "estimate" = est,
              "stat_vals" = (la_mu - B_n)/A_n)

  if (!estimate & !get_all_vals) {
    return(res[[1]])
  } else {
    return(res[c(TRUE, estimate, get_all_vals)])
  }
}

#' Univariate Andrews Test for End-of-Sample Structural Change
#'
#' This implements Andrews' test for end-of-sample change, as described by
#' \insertCite{andrews03;textual}{CPAT}. This test was derived for detecting a
#' change in univariate data. See \insertCite{andrews03}{CPAT} for
#' a description of the test.
#'
#' @param x Vector of the data to test
#' @param M Numeric index of the location of the first potential change point
#' @param pval If \code{TRUE}, return a p-value
#' @param stat If \code{TRUE}, return a test statistic
#' @return If both \code{pval} and \code{stat} are \code{TRUE}, a list
#'         containing both; otherwise, a number for one or the other, depending
#'         on which is \code{TRUE}
#' @references
#'  \insertAllCited{}
#' @examples
#' CPAT:::andrews_test(rnorm(1000), M = 900)
andrews_test <- function(x, M, pval = TRUE, stat = TRUE) {
  mu <- mean(x)
  u <- x - mu
  m <- length(x) - M  # Deriving m and n as described in Andrews (2003)
  n <- length(x) - m

  Sigma <- Reduce("+", lapply(1:(n + 1), function(j)
                              {u[j:(j + m - 1)] %*% t(u[j:(j + m - 1)])})) /
           (n + 1)

  S <- (t(u[(n + 1):(n + m)]) %*% solve(Sigma) %*% u[(n + 1):(n + m)])[1,1]

  submu <- sapply(1:(n - m + 1),
                  function(j) {
                    (sum(x[1:(j - 1)]) + sum(u[(j + ceiling(m/2)):n]))/
                  (n - ceiling(m/2))})

  Sj <- sapply(1:(n - m + 1), function(j) {
                  uj <- x - submu[j]
                  (t(uj[j:(j + m - 1)]) %*% solve(Sigma) %*%
                        uj[j:(j + m - 1)])[1,1]
  })

  res <- list("pval" = mean(S <= Sj), "stat" = S)[c(pval, stat)]
  if (length(res) == 1) {
    return(res[[1]])
  } else {
    return(res)
  }
}

#' Compute the Rényi-Type Statistic
#'
#' This function computes the Rényi-type statistic.
#'
#' The definition of the statistic is
#'
#' \deqn{\max_{t_T \leq t \leq T - t_T} \hat{\sigma}_{t,T}^{-1}
#'       \left|t^{-1}\sum_{s = 1}^{t}X_s - (T - t)^{-1}\sum_{s = t + 1}^{T}
#'       X_s \right|}
#'
#' The parameter \code{kn} corresponds to the trimming parameter \eqn{t_T}.
#'
#' @param dat The data vector
#' @param kn A function corresponding to the trimming parameter \eqn{t_T}; by
#'           default, the square root function
#' @param estimate Set to \code{TRUE} to return the estimated location of the
#'                 change point
#' @param use_kernel_var Set to \code{TRUE} to use kernel methods for long-run
#'                       variance estimation (typically used when the data is
#'                       believed to be correlated); if \code{FALSE}, then the
#'                       long-run variance is estimated using
#'                       \eqn{\hat{\sigma}^2_{T,t} = T^{-1}\left(
#'                       \sum_{s = 1}^t \left(X_s - \bar{X}_t\right)^2 +
#'                       \sum_{s = t + 1}^{T}\left(X_s -
#'                       \tilde{X}_{T - t}\right)^2\right)}, where
#'                       \eqn{\bar{X}_t = t^{-1}\sum_{s = 1}^t X_s} and
#'                       \eqn{\tilde{X}_{T - t} = (T - t)^{-1}
#'                       \sum_{s = t + 1}^{T} X_s}; if \code{custom_var} is not
#'                       \code{NULL}, this argument is ignored
#' @param custom_var Can be a vector the same length as \code{dat} consisting of
#'                   variance-like numbers at each potential change point (so
#'                   each entry of the vector would be the "best estimate" of
#'                   the long-run variance if that location were where the
#'                   change point occured) or a function taking two parameters
#'                   \code{x} and \code{k} that can be used to generate this
#'                   vector, with \code{x} representing the data vector and
#'                   \code{k} the position of a potential change point; if
#'                   \code{NULL}, this argument is ignored
#' @param kernel If character, the identifier of the kernel function as used in
#'               \pkg{cointReg} (see \code{\link[cointReg]{getLongRunVar}}); if
#'               function, the kernel function to be used for long-run variance
#'               estimation (default is the Bartlett kernel in \pkg{cointReg})
#' @param bandwidth If character, the identifier for how to compute the
#'               bandwidth as defined in \pkg{cointReg} (see
#'               \code{\link[cointReg]{getBandwidth}}); if function, a function
#'               to use for computing the bandwidth; if numeric, the bandwidth
#'               value to use (the default is to use Andrews' method, as used in
#'               \pkg{cointReg})
#' @param get_all_vals If \code{TRUE}, return all values for the statistic at
#'                     every tested point in the data set
#' @return If both \code{estimate} and \code{get_all_vals} are \code{FALSE}, the
#'         value of the test statistic; otherwise, a list that contains the test
#'         statistic and the other values requested (if both are \code{TRUE},
#'         the test statistic is in the first position and the estimated change
#'         point in the second)
#' @examples
#' CPAT:::stat_Zn(rnorm(1000))
#' CPAT:::stat_Zn(rnorm(1000), kn = function(n) {floor(log(n))})
#' CPAT:::stat_Zn(rnorm(1000), use_kernel_var = TRUE, bandwidth = "nw",
#'                kernel = "bo")
stat_Zn <- function(dat, kn = function(n) {floor(sqrt(n))}, estimate = FALSE,
                    use_kernel_var = FALSE, custom_var = NULL, kernel = "ba",
                    bandwidth = "and", get_all_vals = FALSE) {
  # Formerly known as statZn()
  if (use_kernel_var) {
    lrv <- get_lrv_vec(dat, kernel, bandwidth)
  } else if (!is.null(custom_var)) {
    use_kernel_var <- TRUE  # Otherwise stat_Zn_cpp() will ignore lrv
    if (is.function(custom_var)) {
      # This may seem silly, but this is so that error codes refer to
      # custom_var, and we don't want recursion either
      custom_var_temp <- custom_var
      custom_var <- purrr::partial(custom_var_temp, x = dat, .lazy = FALSE)
      custom_var_vec <- vapply(1:length(dat), custom_var,
                               FUN.VALUE = numeric(1))
    } else if (is.numeric(custom_var)) {
      if (length(custom_var) < length(dat)) stop("custom_var must have" %s%
                                                 "length at least" %s%
                                                 length(dat) %s0% ", the" %s%
                                                 "length of the data set")
      custom_var_vec <- custom_var
    } else {
      stop("Don't know how to handle custom_var of class" %s% class(custom_var))
    }
    if (any(custom_var_vec < 0)) stop("custom_var suggests a negative" %s%
                                      "variance, which is impossible.")
    lrv <- custom_var_vec
  } else {
    # A vector must be passed to stat_Zn_cpp, so we'll pass one with an
    # impossible value
    lrv <- c(-1)
  }

  res <- stat_Zn_cpp(dat, kn(length(dat)), use_kernel_var, lrv, get_all_vals)
  res[[2]] <- as.integer(res[[2]])
  if (!estimate & !get_all_vals) {
    return(res[[1]])
  } else {
    return(res[c(TRUE, estimate, get_all_vals)])
  }
  # Here is equivalen (slow) R code
  #n = length(dat)
  #return(sqrt(kn(n))*max(sapply(
  #    floor(kn(n)):(n - floor(kn(n))), function(k)
  #      abs(1/k*sum(dat[1:k]) -
  #            1/(n-k)*sum(dat[(k+1):n]))/sqrt((sum((dat[1:k] -
  #             mean(dat[1:k]))^2)+sum((dat[(k+1):n] -
  #             mean(dat[(k+1):n]))^2))/n))))
}

#' Multivariate Andrews' Test for End-of-Sample Structural Change
#'
#' This implements Andrews' test for end-of-sample change, as described by
#' \insertCite{andrews03;textual}{CPAT}. This test was derived for detecting a
#' change in multivarate data, aso originally described. See
#' \insertCite{andrews03}{CPAT} for a description of the test.
#'
#' @param formula The regression formula, which will be passed to
#'        \code{\link[stats]{lm}}
#' @param data \code{data.frame} containing the data
#' @inheritParams andrews_test
#' @return If both \code{pval} and \code{stat} are \code{TRUE}, a list
#'         containing both; otherwise, a number for one or the other, depending
#'         on which is \code{TRUE}
#' @references
#'  \insertAllCited{}
#' @examples
#' x <- rnorm(1000)
#' y <- 1 + 2 * x + rnorm(1000)
#' df <- data.frame(x, y)
#' CPAT:::andrews_test_reg(y ~ x, data = df, M = 900)
andrews_test_reg <- function(formula, data, M, pval = TRUE, stat = TRUE) {
  if (!is(formula, "formula")) stop("Bad formula passed to argument" %s%
                                    "\"formula\"")
  fit <- lm(formula = formula, data = data)
  beta <- coefficients(fit)
  d <- length(beta)
  X <- model.matrix(fit)
  u <- residuals(fit)
  y <- fit$model[[1]]
  m <- nrow(X) - M  # Deriving m and n as described in Andrews (2003)
  n <- nrow(X) - m

  Sigma <- Reduce("+", lapply(1:(n + 1), function(j)
                              {u[j:(j + m - 1)] %*% t(u[j:(j + m - 1)])})) /
           (n + 1)

  if (d <= m) {
    V <- t(X[(n + 1):(n + m),]) %*% solve(Sigma) %*% X[(n + 1):(n + m),]
    A <- t(X[(n + 1):(n + m),]) %*% solve(Sigma) %*% u[(n + 1):(n + m)]
    S <- (t(A) %*% solve(V) %*% A)[1,1]
  } else {
    S <- (t(u[(n + 1):(n + m)]) %*% solve(Sigma) %*% u[(n + 1):(n + m)])[1,1]
  }

  subbeta <- lapply(1:(n - m + 1),
                  function(j) {
                    coefficients(lm(formula = formula,
                                    data = data[c(1:(j - 1),
                                                  (j + ceiling(m/2)):n),]))
                  })

  Sj <- sapply(1:(n - m + 1), function(j) {
                 yj <- y[j:(j + m - 1)]
                 Xj <- X[j:(j + m - 1),]
                 uj <- yj - Xj %*% subbeta[[j]]

                 Vj <- t(Xj) %*%
                   solve(Sigma) %*% Xj
                 Aj <- t(Xj) %*%
                   solve(Sigma) %*% uj
                 (t(Aj) %*% solve(Vj) %*% A)[1,1]  # Sj
  })

  res <- list("pval" = mean(S <= Sj), "stat" = S)[c(pval, stat)]
  if (length(res) == 1) {
    return(res[[1]])
  } else {
    return(res)
  }
}

################################################################################
# STATISTICAL TEST INTERFACES
################################################################################

#' CUSUM Test
#'
#' Performs the (univariate) CUSUM test for change in mean, as described in
#' \insertCite{horvathricemiller19}{CPAT}. This is effectively an interface to
#' \code{\link{stat_Vn}}; see its documentation for more details. p-values are
#' computed using \code{\link{pkolmogorov}}, which represents the limiting
#' distribution of the statistic under the null hypothesis.
#'
#' @param x Data to test for change in mean
#' @param stat_plot Whether to create a plot of the values of the statistic at
#'        all potential change points
#' @inheritParams stat_Vn
#' @return A \code{htest}-class object containing the results of the test
#' @references
#'  \insertAllCited{}
#' @examples
#' CUSUM.test(rnorm(1000))
#' CUSUM.test(rnorm(1000), use_kernel_var = TRUE, kernel = "bo",
#'            bandwidth = "nw")
#' @export
CUSUM.test <- function(x, use_kernel_var = FALSE, stat_plot = FALSE,
                     kernel = "ba", bandwidth = "and") {
  testobj <- list()
  testobj$method <- "CUSUM Test for Change in Mean"
  testobj$data.name <- deparse(substitute(x))

  res <- stat_Vn(x,
                 estimate = TRUE,
                 use_kernel_var = use_kernel_var,
                 kernel = kernel,
                 bandwidth = bandwidth,
                 get_all_vals = stat_plot)
  stat <- res[[1]]
  est <- res[[2]]

  if (stat_plot) {
    plot.ts(res[[3]], main = "Value of Test Statistic", ylab = "Statistic")
  }
  attr(stat, "names") <- "A"
  attr(est, "names") <- "t*"
  testobj$p.value <- 1 - pkolmogorov(res[[1]])
  testobj$estimate <- est
  testobj$statistic <- stat

  class(testobj) <- "htest"
  testobj
}

#' Darling-Erdös Test
#'
#' Performs the (univariate) Darling-Erdös test for change in mean, as described
#' in \insertCite{horvathricemiller19}{CPAT}. This is effectively an interface
#' to \code{\link{stat_de}}; see its documentation for more details. p-values
#' are computed using \code{\link{pdarling_erdos}}, which represents the
#' limiting distribution of the test statistic under the null hypothesis when
#' \code{a} and \code{b} are chosen appropriately. (Change those parameters at
#' your own risk!)
#'
#' @param x Data to test for change in mean
#' @param stat_plot Whether to create a plot of the values of the statistic at
#'                  all potential change points
#' @inheritParams stat_de
#' @return A \code{htest}-class object containing the results of the test
#' @references
#'  \insertAllCited{}
#' @examples
#' DE.test(rnorm(1000))
#' DE.test(rnorm(1000), use_kernel_var = TRUE, kernel = "bo", bandwidth = "nw")
#' @export
DE.test <- function(x, a = log, b = log, use_kernel_var = FALSE,
                    stat_plot = FALSE, kernel = "ba", bandwidth = "and") {
  l <- function(x) {sqrt(2*log(x))}
  u <- function(x) {2*log(x) + 1/2*log(log(x)) - 1/2*log(pi)}

  testobj <- list()
  testobj$method <- "Darling-Erdos Test for Change in Mean"
  testobj$data.name <- deparse(substitute(x))
  params <- c(l(a(length(x))), u(b(length(x))))
  names(params) <- c("a(" %s0% deparse(substitute(a)) %s0% "(T))", "b(" %s0%
                     deparse(substitute(b)) %s0% "(T))")
  testobj$parameter <- params

  res <- stat_de(x,
                 estimate = TRUE,
                 a = a,
                 b = b,
                 use_kernel_var = use_kernel_var,
                 kernel = kernel,
                 bandwidth = bandwidth,
                 get_all_vals = stat_plot)
  stat <- res[[1]]
  est <- res[[2]]

  if (stat_plot) {
    plot.ts(res[[3]], main = "Value of Test Statistic", ylab = "Statistic")
  }
  attr(stat, "names") <- "A"
  attr(est, "names") <- "t*"
  testobj$p.value <- 1 - pdarling_erdos(res[[1]])
  testobj$estimate <- est
  testobj$statistic <- stat

  class(testobj) <- "htest"
  testobj
}

#' Rényi-Type Test
#'
#' Performs the (univariate) Rényi-type test for change in mean, as described in
#' \insertCite{horvathricemiller19}{CPAT}. This is effectively an interface to
#' \code{\link{stat_Zn}}; see its documentation for more details. p-values are
#' computed using \code{\link{pZn}}, which represents the limiting distribution
#' of the test statistic under the null hypothesis, which represents the
#' limiting distribution of the test statistic under the null hypothesis when
#' \code{kn} represents a sequence \eqn{t_T} satisfying \eqn{t_T \to \infty}
#' and \eqn{t_T/T \to 0} as \eqn{T \to \infty}. (\code{\link[base]{log}} and
#' \code{\link[base]{sqrt}} should be good choices.)
#'
#' @param x Data to test for change in mean
#' @param stat_plot Whether to create a plot of the values of the statistic at
#'        all potential change points
#' @inheritParams stat_Zn
#' @return A \code{htest}-class object containing the results of the test
#' @references
#'  \insertAllCited{}
#' @examples
#' HR.test(rnorm(1000))
#' HR.test(rnorm(1000), use_kernel_var = TRUE, kernel = "bo", bandwidth = "nw")
#' @export
HR.test <- function(x, kn = log, use_kernel_var = FALSE, stat_plot = FALSE,
                  kernel = "ba", bandwidth = "and") {
  testobj <- list()
  testobj$method <- "Horvath-Rice Test for Change in Mean"
  testobj$data.name <- deparse(substitute(x))

  res <- stat_Zn(x,
                 kn = kn,
                 estimate = TRUE,
                 use_kernel_var = use_kernel_var,
                 kernel = kernel,
                 bandwidth = bandwidth,
                 get_all_vals = stat_plot)
  stat <- res[[1]]
  est <- res[[2]]

  if (stat_plot) {
    series <- ts(res[[3]], start = ceiling(kn(length(x))))
    plot.ts(res[[3]], main = "Value of Test Statistic", ylab = "Statistic")
  }

  kn_val <- kn(length(x))
  attr(kn_val, "names") <- deparse(substitute(kn)) %s0% "(T)"

  attr(stat, "names") <- "D"
  attr(est, "names") <- "t*"
  testobj$parameters <- kn_val
  testobj$p.value <- 1 - pZn(res[[1]])
  testobj$estimate <- est
  testobj$statistic <- stat

  class(testobj) <- "htest"
  testobj
}

#' Hidalgo-Seo Test
#'
#' Performs the (univariate) Hidalgo-Seo test for change in mean, as described
#' in \insertCite{horvathricemiller19}{CPAT}. This is effectively an interface
#' to \code{\link{stat_hs}}; see its documentation for more details. p-values
#' are computed using \code{\link{phidalgo_seo}}, which represents the limiting
#' distribution of the test statistic when the null hypothesis is true.
#'
#' @param x Data to test for change in mean
#' @param stat_plot Whether to create a plot of the values of the statistic at
#'        all potential change points
#' @inheritParams stat_hs
#' @return A \code{htest}-class object containing the results of the test
#' @references
#'  \insertAllCited{}
#' @examples
#' HS.test(rnorm(1000))
#' HS.test(rnorm(1000), corr = FALSE)
#' @export
HS.test <- function(x, corr = TRUE, stat_plot = FALSE) {
  testobj <- list()
  testobj$method <- "Hidalgo-Seo Test for Change in Mean"
  testobj$data.name <- deparse(substitute(x))
  params <- c(corr)
  names(params) <- c("Correlated Residuals")
  testobj$parameter <- params

  res <- stat_hs(x, estimate = TRUE, corr = corr, get_all_vals = stat_plot)
  stat <- res[[1]]
  est <- res[[2]]

  if (stat_plot) {
    plot.ts(res[[3]], main = "Value of Test Statistic", ylab = "Statistic")
  }
  attr(stat, "names") <- "A"
  attr(est, "names") <- "t*"
  testobj$p.value <- 1 - phidalgo_seo(res[[1]])
  testobj$estimate <- est
  testobj$statistic <- stat

  class(testobj) <- "htest"
  testobj
}

#' Andrews' Test for End-of-Sample Structural Change
#'
#' Performs Andrews' test for end-of-sample structural change, as described in
#' \insertCite{andrews03}{CPAT}. This function works for both univariate and
#' multivariate data depending on the nature of \code{x} and whether
#' \code{formula} is specified. This function is thus an interface to
#' \code{\link{andrews_test}} and \code{\link{andrews_test_reg}}; see the
#' documentation of those functions for more details.
#'
#' @param x Data to test for change in mean (either a vector or
#'          \code{data.frame})
#' @inheritParams andrews_test_reg
#' @return A \code{htest}-class object containing the results of the test
#' @references
#'  \insertAllCited{}
#' @examples
#' Andrews.test(rnorm(1000), M = 900)
#' x <- rnorm(1000)
#' y <- 1 + 2 * x + rnorm(1000)
#' df <- data.frame(x, y)
#' Andrews.test(df, y ~ x, M = 900)
#' @export
Andrews.test <- function(x, M, formula = NULL) {
  testobj <- list()
  testobj$method <- "Andrews' Test for Structural Change"
  testobj$data.name <- deparse(substitute(x))

  if (is.numeric(x)) {
    mchange <- length(x) - M
    res <- andrews_test(x, M, pval = TRUE, stat = TRUE)
  } else if (is.data.frame(x)) {
    mchange <- nrow(x) - M
    res <- andrews_test_reg(formula, x, M, pval = TRUE, stat = TRUE)
  } else {
    stop("x must be vector-like or a data frame")
  }

  stat <- res[["stat"]]
  attr(mchange, "names") <- "m"

  attr(stat, "names") <- "S"
  testobj$p.value <- res[["pval"]]
  testobj$parameters <- mchange
  testobj$statistic <- stat

  class(testobj) <- "htest"
  testobj
}
