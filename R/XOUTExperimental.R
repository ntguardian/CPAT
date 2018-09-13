################################################################################
# XOUTExperimental.R
################################################################################
# 2018-08-27
# Curtis Miller
################################################################################
# Experimental functions that are unstable (perhaps they don't even work) and
# should not be included in any public distribution of CPAT.
################################################################################

################################################################################
# ESTIMATORS AND STATISTICAL FUNCTIONS
################################################################################

#' Get \eqn{\Sigma} and \eqn{\Omega} Matrices
#'
#' This function obtains the \eqn{\Sigma(k)} and \eqn{\Omega(k)} matrices as
#' described in \insertCite{ling07;textual}{CPAT}.
#'
#' @param grad1 A matrix consisting of all gradients at each \eqn{t} up to
#'              \eqn{k}
#' @param hess1 A 3D array consisting of all Hessian matrices at each \eqn{t} up
#'              to \eqn{k}
#' @param grad2 A matrix consisting of all gradients at each \eqn{t} after
#'              \eqn{k}; if \code{NULL}, only \code{grad1} will be used
#' @param hess2 A 3D array consisting of all Hessian matrices at each \eqn{t}
#'              after \eqn{k}; if \code{NULL}, only \code{hess1} will be used
#' @return A list containing a matrix \code{Omega} and a matrix \code{Sigma}
#'         corresponding to the respective matrices described in
#'         \insertCite{ling07;textual}{CPAT}
#' @references
#'  \insertAllCited{}
#' @examples
#' sig <- c(0.68, 18.60,  0.23,  2.10,  1.47, 2.92,  0.39, 40.31, 0.60,  0.18)
#' x   <- c(0.74, -0.97, -0.73, -1.02, -0.80, 0.79, -0.73, -3.06, 0.75, -0.48)
#' mats <- get_gradient_hessian(var = sig^2, eps = x^2, omega = 0.5,
#'                              alpha = 0.1, beta = 0.1)
#' get_sigma_omega(mats$gradient, mats$hessian)
get_sigma_omega <- function(grad1, hess1, grad2 = NULL, hess2 = NULL) {
  if (is.null(grad2) & is.null(hess2)) {
    # Compute Sigma(n) and Omega(n), that is, for no change point
    DDT <- array(apply(grad1, 2, function(col) {col %*% t(col)}),
                      dim = c(3, 3, ncol(grad1)))
    Omega <- apply(DDT, c(1, 2), sum)
    Sigma <- -apply(hess1, c(1, 2), sum)
  } else {
    # Compute Sigma(k) and Omega(k)
    DDT1 <- array(apply(grad1, 2, function(col) {col %*% t(col)}),
                      dim = c(3, 3, ncol(grad1)))
    DDT2 <- array(apply(grad2, 2, function(col) {col %*% t(col)}),
                      dim = c(3, 3, ncol(grad2)))
    Omega <- apply(DDT1, c(1, 2), sum) + apply(DDT2, c(1, 2), sum)
    Sigma <- -apply(hess1, c(1, 2), sum) - apply(hess2, c(1, 2), sum)
  }
  
  list("Omega" = Omega, "Sigma" = Sigma)
}

#' Theoretical Covariance Matrix for GARCH(1,1) Models
#'
#' This computes the covariance matrix appropriate to GARCH(1,1) models using
#' the  "theoretical" gradient and Hessian matrix (as opposed to numerical
#' approximations).
#'
#' @param x A series corresponding to the GARCH(1,1) process for which to
#'          obtain the matrix
#' @param k An integer designating where to split the data (a potential change
#'          point)
#' @param return_coef If \code{TRUE}, return a list that also include the fitted
#'                    coefficients
#' @param ... Arguments to pass to \code{\link[rugarch]{ugarchfit}}
#' @return A matrix corresponding to the covariance matrix (or a list including
#'         this as an entry if \code{return_coef} is \code{TRUE})
#' @examples
#' garch11_covmat_theo(rnorm(1000))
garch11_covmat_theo <- function(x, k = NULL, return_coef = FALSE, ...) {
  # TODO: curtis: IMPORT rugarch NAMESPACE -- Wed 29 Aug 2018 03:48:35 PM MDT

  arglist <- list(...)
  arglist$spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0),
                                               include.mean = FALSE))
  x <- as.vector(x)
  n <- length(x)
  
  if (is.null(k)) {
    arglist$data <- x
    
    fit <- do.call(ugarchfit, arglist)
    coef <- fit@fit$coef
    
    derivs <- get_gradient_hessian(var = fit@fit$sigma^2,
                                   eps = x^2,
                                   omega = coef["omega"],
                                   alpha = coef["alpha1"],
                                   beta = coef["beta1"])
    mats <- get_sigma_omega(grad1 = derivs$gradient, hess1 = derivs$hessian)
  } else {
    x1 <- x[1:k]
    x2 <- x[(k + 1):n]
    arglist1 <- arglist
    arglist2 <- arglist
    arglist1$data <- x1
    arglist2$data <- x2
    
    fit1 <- do.call(ugarchfit, arglist1)
    fit2 <- do.call(ugarchfit, arglist2)
    coef1 <- fit1@fit$coef
    coef2 <- fit2@fit$coef
    
    derivs1 <- get_gradient_hessian(var = fit1@fit$sigma^2,
                                    eps = x1^2,
                                    omega = coef1["omega"],
                                    alpha = coef1["alpha1"],
                                    beta = coef1["beta1"])
    derivs2 <- get_gradient_hessian(var = fit2@fit$sigma^2,
                                    eps = x2^2,
                                    omega = coef2["omega"],
                                    alpha = coef2["alpha1"],
                                    beta = coef2["beta1"])
    mats <- get_sigma_omega(grad1 = derivs1$gradient, hess1 = derivs1$hessian,
                            grad2 = derivs2$gradient, hess2 = derivs2$hessian)
  }
  Sigma = mats$Sigma
  Omega = mats$Omega
  
  C <- solve(Sigma) %*% Omega %*% solve(Sigma)
  rownames(C) <- c("omega", "alpha1", "beta1")
  colnames(C) <- rownames(C)
  
  if (return_coef) {
    if (is.null(k)) {
      return(list("coef" = coef, "C" = C))
    } else {
      return(list("coef1" = coef1, "coef2" = coef2, "C" = C))
    }
  } else {
    return(C)
  }
}

#' Compute the Gradient and Hessian Matrix For GARCH(1,1) Model
#'
#' Effectively an interface to a C++ function; gets the gradient and Hessian
#' matrix for a GARCH(1,1) process of the form:
#'
#' \deqn{\epsilon_t^2 = \omega + \alpha \epsilon_{t - 1}^2 +
#' \beta \sigma_{t - 1}^2}
#'
#' @param var A vector of conditional variances
#' @param eps A vector of the \strong{squared} GARCH(1,1) process
#' @param omega The \eqn{\omega} parameter in the GARCH(1,1) model
#' @param alpha The \eqn{\alpha} parameter in the GARCH(1,1) model
#' @param beta The \eqn{\beta} parameter in the GARCH(1,1) model
#' @param init_vals A vector of initial values to be used when computing the
#'                  gradient and Hessian matrix
#' @return A list with the \code{gradient} element containing the gradient and
#'         \code{hessian} containing the Hessian matrix; actually, this is a
#'         matrix and a 3D array (respectively) that contain gradient and
#'         Hessian matrices at each \eqn{t}, to be condensed later
#' @examples
#' sig <- c(0.68, 18.60,  0.23,  2.10,  1.47, 2.92,  0.39, 40.31, 0.60,  0.18)
#' x   <- c(0.74, -0.97, -0.73, -1.02, -0.80, 0.79, -0.73, -3.06, 0.75, -0.48)
#' get_gradient_hessian(var = sig^2, eps = x^2, omega = 0.5, alpha = 0.1,
#'                      beta = 0.1)
get_gradient_hessian <- function(var, eps, omega, alpha, beta, 
                                 init_vals = NULL) {
  if (is.null(init_vals)) {
    init_vals_cpp <- c(1/(1 - alpha - beta),
                       omega/(1 - alpha - beta)^2,
                       omega/(1 - alpha - beta)^2,
                       1/(1 - alpha - beta)^2,
                       1/(1 - alpha - beta)^2,
                       2 * omega/(1 - alpha - beta)^3,
                       2 * omega/(1 - alpha - beta)^3,
                       2 * omega/(1 - alpha - beta)^3)
  } else {
    init_vals_cpp <- c(1,
                       init_vals["eps"],
                       init_vals["var"],
                       0,
                       0,
                       0,
                       0,
                       0)
  }
  
  var <- as.vector(var)
  eps <- as.vector(eps)
  
  cond_var_derivatives <- cond_var_gradient_hessian_cpp(var, eps, omega, alpha,
                                                        beta, init_vals_cpp)
  
  lt_grad <- (var - eps)/var^2
  gradient <- cond_var_derivatives$gradient * rep(lt_grad, each = 3)
  
  lt_hess <- (2 * eps - var)/var^3
  outer_grad <- array(apply(cond_var_derivatives$gradient, 2,
                            function(col) {col %*% t(col)}),
                      dim = c(3, 3, length(lt_grad)))
  
  hessian <- rep(lt_hess, each = 9) * outer_grad #+
    rep(lt_grad, each = 9) * cond_var_derivatives$hessian
  
  list("gradient" = gradient, "hessian" = hessian)
}

################################################################################
# EXPERIMENTAL STATISTICS
################################################################################

#' Multivariate Rényi-Type Statistic
#'
#' Computes the Rényi-type statistic in the context of structural change in a
#' regression model, specifically a linear model.
#'
#' This statistic performs badly; I've almost never seen the statistic fail to
#' reject the null hypothesis when it should. Why this occurs is a mystery to me
#' (Curtis).
#'
#' @param dat A \code{timeSeries} or \code{data.frame} object containing the
#'            variables in the model; eligible to be passed to the argument
#'            \code{data} in \code{\link[stats]{lm}}
#' @param kn A function that corresponds to the sequence \eqn{t_T} in the
#'           definition of the Rényi-type statistc; by default, is \eqn{\lfloor
#'           \sqrt{T} \rfloor}
#' @param estimate If \code{TRUE}, return the estimated location of the change
#'                 point
#' @param verbose If \code{TRUE}, print progress reports to output
#' @param wholevcov If \code{TRUE}, use the covariance matrix estimated on the
#'                  whole sample
#' @param DEBUG If \code{TRUE}, turn on debug mode
#' @param ... Further arguments to pass to \code{\link[stats]{lm}} (such as the
#'            regression formula)
#' @return A vector containing the test statistic and the location of the
#'         change, if requested
#' @examples
#' x <- rnorm(1000)
#' y <- 1 + 2 * x + rnorm(1000)
#' df <- data.frame(x, y)
#' stat_Zn_model_reg(dat = df, formula = y ~ x)
stat_Zn_model_reg <- function(dat, kn = function(n) {floor(sqrt(n))},
                              estimate = FALSE, verbose = FALSE, 
                              wholevcov = TRUE, DEBUG = FALSE, ...) {
  # Formerly known as statZnModelReg()
  args <- list(...)

  # if (is.vector(dat) | is.ts(dat) | is.matrix(dat)) {
  #   st = deparse(substitute(dat))
  #   dat <- data.frame(dat)
  #   names(dat) <- st
  # }
  n <- nrow(dat)
  k <- kn(n)

  args$data <- dat
  # Get the covariance matrix
  temp <- do.call(lm, args = args)
  # TODO: curtis: IMPORT sandwich NAMESPACE -- Wed 29 Aug 2018 11:28:59 PM MDT
  C <- vcovHC(temp, type = "HC0") * n

  if (verbose) {
    cat("\n")
  }
  if (DEBUG) {
    DEBUG_list <- list()
  }
  # Get all values of the test statistic (will take max later)
  vals <- sapply(k:(n - k), function(t) {
    if (verbose) {
      cat(sprintf("\rChange points checked: %%%3.0f", (t - k)/(n - 2 * k) * 100))
    }
    # Split indices
    idx1 <- 1:t
    idx2 <- (t + 1):n

    # Parameter estimates around split
    args1 <- args
    args2 <- args
    args1$data <- dat[idx1,]
    args2$data <- dat[idx2,]
    mod1 <- do.call(lm, args = args1)
    mod2 <- do.call(lm, args = args2)
    coef1 <- coef(mod1)
    coef2 <- coef(mod2)

    if (wholevcov) {
      Ck <- C
    } else {
      H1 <- bread(mod1)
      H2 <- bread(mod2)

      GtG1 <- meatHC(mod1, type = "HC0")
      GtG2 <- meatHC(mod2, type = "HC0")

      H <- solve(solve(H1) + solve(H2))
      GtG <- GtG1 + GtG2

      Ck <- H %*% GtG %*% H * 2
    }

    if (DEBUG) {
      DEBUG_list[[t]] <<- list("coef1" = coef1, "coef2" = coef2, "Ck" = Ck)
    }
    return(t(coef1 - coef2) %*% solve(Ck) %*% (coef1 - coef2))
  })
  if (verbose) {
    cat("\n")
  }

  # Get maximal value
  tau <- (k:(n - k))[which(vals == max(vals, na.rm = TRUE))]
  stat <- sqrt(k * max(vals, na.rm = TRUE))

  res <- list("stat" = stat)
  if (estimate) {
    res$estimate <- round(tau)
  }
  if (DEBUG) {
    res$quad_term <- sqrt(k * vals)
    res$DEBUG_list <- DEBUG_list
  }

  res
}

#' GARCH Model Rényi-Type Statistic
#'
#' Computes the Rényi-type statistic in the context of structural change of an
#' ARFIMA-GARCH model.
#'
#' This statistic performs badly; I've almost never seen the statistic fail to
#' reject the null hypothesis when it should. Why this occurs is a mystery to me
#' (Curtis).
#'
#' This function should be rewritten to depend on \pkg{rugarch}, though that
#' implementation isn't likely to be stable either. (\pkg{fGarch} is no longer
#' being maintained.)
#'
#' TODO: curtis: DOCUMENTATION IS ALL LIES!!! -- Wed 29 Aug 2018 02:41:21 PM MDT
#'
#' @param dat A \code{timeSeries} or \code{data.frame} object containing the
#'            variables in the model; eligible to be passed to data in
#'            \code{\link[fGarch]{garchFit}}
#' @param kn A function that corresponds to the sequence \eqn{t_T} in the
#'           definition of the Rényi-type statistc; by default, is \eqn{\lfloor
#'           \sqrt{T} \rfloor}
#' @param estimate If \code{TRUE}, return the estimated location of the change
#'                 point
#' @param verbose If \code{TRUE}, print progress reports to output
#' @param wholevcov If \code{TRUE}, use the covariance matrix estimated on the
#'                  whole sample
#' @param DEBUG If \code{TRUE}, turn on debug mode
#' @param ... Further arguments to pass to \code{\link[stats]{lm}} (such as the
#'            regression formula)
#' @return A vector containing the test statistic and the location of the
#'         change, if requested
#' @examples
#' stat_Zn_model(rnorm(1000))
#' stat_Zn_model(rnorm(1000), formula = ~ garch(1, 1))
stat_Zn_model <- function(dat, kn = function(n) {floor(sqrt(n))},
                          estimate = FALSE, verbose = FALSE, wholevcov = TRUE,
                          use_theo = FALSE, DEBUG = FALSE, ...) {
  # Formerly called statZnModel()

  # TODO: curtis: THIS SHOULD DEPEND ON rugarch -- Wed 29 Aug 2018 12:34:16 PM
  args <- list(...)

  # if (is.vector(dat) | is.ts(dat) | is.matrix(dat)) {
  #   st = deparse(substitute(dat))
  #   dat <- data.frame(dat)
  #   names(dat) <- st
  # }
  n <- length(dat)
  k <- kn(n)

  args$data <- dat
  # Get the covariance matrix
  temp <- do.call(ugarchfit, args = args)
  if (use_theo) {
    C <- garch11_covmat_theo(dat)
  } else {
    C <- vcov(temp, robust = TRUE) * n  # THIS IS WHERE THE ERROR WAS!
  }

  if (verbose) {
    cat("\n")
  }
  if (DEBUG) {
    DEBUG_list <- list()
  }
  # Get all values of the test statistic (will take max later)
  vals <- sapply(k:(n - k), function(t) {
    if (verbose) {
      cat(sprintf("\rChange points checked: %%%3.0f", (t - k)/(n - 2 * k) * 100))
    }
    if (use_theo) {
      args1 <- args
      args1$x <- dat
      args1$k <- t
      args1$return_coef <- TRUE
      res <- do.call(garch11_covmat_theo, args1)

      Ck <- res$C
      coef1 <- res$coef1
      coef2 <- res$coef2
    } else {
      # Split indices
      idx1 <- 1:t
      idx2 <- (t + 1):n

      # Parameter estimates around split
      args1 <- args
      args2 <- args
      args1$data <- dat[idx1]
      args2$data <- dat[idx2]
      mod1 <- do.call(ugarchfit, args = args1)
      mod2 <- do.call(ugarchfit, args = args2)
      coef1 <- coef(mod1)
      coef2 <- coef(mod2)
    }

    if (wholevcov) {
      Ck <- C
    } else if (!use_theo) {
      H1 <- mod1@fit$hessian
      H2 <- mod2@fit$hessian

      V1 <- vcov(mod1, robust = TRUE) * t
      V2 <- vcov(mod2, robust = TRUE) * (n - t)

      GtG1 <- H1 %*% V1 %*% H1
      GtG2 <- H2 %*% V2 %*% H2

      H <- H1 + H2
      GtG <- GtG1 + GtG2

      Ck <- solve(H) %*% GtG %*% solve(H)
    }

    if (DEBUG) {
      DEBUG_list[[t]] <<- list("coef1" = coef1, "coef2" = coef2, "Ck" = Ck)
    }
    return(t(coef1 - coef2) %*% solve(Ck) %*% (coef1 - coef2))
  })
  if (verbose) {
    cat("\n")
  }

  # Get maximal value
  tau <- (k:(n - k))[which(vals == max(vals))]
  stat <- sqrt(k * max(vals))

  res <- list("stat" = stat)
  if (estimate) {
    res$estimate <- tau
  }
  if (DEBUG) {
    res$quad_term <- sqrt(k * vals)
    res$DEBUG_list <- DEBUG_list
  }

  res
}

#' Ling's Multivariate Darling-Erdös-Type Statistic
#'
#' Computes a statistic for testing for structural change in a regression model
#' as described in \insertCite{ling07;textual}{CPAT}.
#'
#' This statistic seems to perform okay. (Curtis)
#'
#' @param dat A \code{timeSeries} or \code{data.frame} object containing the
#'            variables in the model; eligible to be passed to the argument
#'            \code{data} in \code{\link[stats]{lm}}
#' @param kn A function that corresponds to the \eqn{k_n} function in the
#'           statistic; by default, is \eqn{\lfloor \sqrt{n} \rfloor}
#' @param estimate If \code{TRUE}, return the estimated location of the change
#'                 point
#' @param verbose If \code{TRUE}, print progress reports to output
#' @param wholevcov If \code{TRUE}, use the covariance matrix estimated on the
#'                  whole sample
#' @param DEBUG If \code{TRUE}, turn on debug mode
#' @param ... Further arguments to pass to \code{\link[stats]{lm}} (such as the
#'            regression formula)
#' @return A vector containing the test statistic and the location of the
#'         change, if requested
#' @references
#'  \insertAllCited{}
#' @examples
#' x <- rnorm(1000)
#' y <- 1 + 2 * x + rnorm(1000)
#' df <- data.frame(x, y)
#' stat_de_model_reg(dat = df, formula = y ~ x)
stat_de_model_reg <- function(dat, kn = function(n) {floor(sqrt(n))},
                              estimate = FALSE, verbose = FALSE, DEBUG = TRUE,
                              ...) {
  # Formerly called statDEModelReg()
  args <- list(...)
  n <- nrow(dat)
  if (is.function(kn)) {
    k <- kn(n)
  } else {
    k <- kn
  }

  # Constants involved in the statistic
  bn <- (2 * log(log(n)) + k * log(log(log(n)))/2 - log(gamma(k/2)))^2/(
    2*log(log(n)))
  an <- sqrt(bn/(2*log(log(n))))

  args$data <- dat

  if (verbose) {
    cat("\n")
  }
  if (DEBUG) {
    DEBUG_list <- list()
  }
  # Get all values of the test statistic (will take max later)
  vals <- sapply(k:(n - k), function(t) {
    if (verbose) {
      cat(sprintf("\rChange points checked: %%%3.0f", (t - k)/(n - 2 * k) * 100))
    }
    # Split indices
    idx1 <- 1:t
    idx2 <- (t + 1):n

    # Parameter estimates around split
    args1 <- args
    args2 <- args
    args1$data <- dat[idx1,]
    args2$data <- dat[idx2,]
    mod1 <- do.call(lm, args = args1)
    mod2 <- do.call(lm, args = args2)
    coef1 <- coef(mod1)
    coef2 <- coef(mod2)

    # TODO: curtis: IMPORT sandwich NAMESPACE -- Wed 29 Aug 2018 11:30:05 PM MDT
    H1 <- bread(mod1)
    H2 <- bread(mod2)

    GtG1 <- meatHC(mod1, type = "HC0")
    GtG2 <- meatHC(mod2, type = "HC0")

    H <- solve(solve(H1) + solve(H2))
    GtG <- GtG1 + GtG2

    Ck <- H %*% GtG %*% H * 2

    if (DEBUG) {
      DEBUG_list[[t]] <<- list("coef1" = coef1, "coef2" = coef2, "Ck" = Ck)
    }
    return(t(coef1 - coef2) %*% solve(Ck) %*% (coef1 - coef2))
  })
  if (verbose) {
    cat("\n")
  }

  # Get maximal value
  t <- k:(n - k)
  vals <- t * (n - t) * vals / n^2
  tau <- (k:(n - k))[which(vals == max(vals))]

  stat <- max((vals - bn)/an)

  res <- list("stat" = stat)
  if (estimate) {
    res$estimate <- tau
  }
  if (DEBUG) {
    res$quad_term <- (vals - bn)/an
    res$DEBUG_list <- DEBUG_list
  }

  res
}

#' Ling's GARCH Model Darling-Erdös-Type Statistic
#'
#' Computes \insertCite{ling07}{CPAT}'s statistic for structural change of an
#' ARFIMA-GARCH model.
#'
#' This statistic performs badly; I've almost never seen the statistic fail to
#' reject the null hypothesis when it should. Why this occurs is a mystery to me
#' (Curtis).
#'
#' This function should be rewritten to depend on \pkg{rugarch}, though that
#' implementation isn't likely to be stable either. (\pkg{fGarch} is no longer
#' being maintained.)
#'
#' TODO: curtis: DOCUMENTATION IS ALL LIES!!! -- Wed 29 Aug 2018 02:41:21 PM MDT
#'
#' @param dat A \code{timeSeries} or \code{data.frame} object containing the
#'            variables in the model; eligible to be passed to data in
#'            \code{\link[fGarch]{garchFit}}
#' @param kn A function that corresponds to the sequence \eqn{t_T} in the
#'           definition of the Rényi-type statistc; by default, is \eqn{\lfloor
#'           \sqrt{T} \rfloor}
#' @param estimate If \code{TRUE}, return the estimated location of the change
#'                 point
#' @param verbose If \code{TRUE}, print progress reports to output
#' @param use_theo If \code{TRUE}, use the theoretically derived gradient and
#'                 Hessian matrices
#' @param DEBUG If \code{TRUE}, turn on debug mode
#' @param ... Further arguments to pass to \code{\link[stats]{lm}} (such as the
#'            regression formula)
#' @return A vector containing the test statistic and the location of the
#'         change, if requested
#' @references
#'  \insertAllCited{}
#' @examples
#' stat_de_model(rnorm(1000))
#' stat_de_model(rnorm(1000), formula = ~ garch(1, 1))
stat_de_model <- function(dat, kn = function(n) {floor(sqrt(n))},
                          estimate = FALSE, verbose = FALSE, use_theo = TRUE,
                          DEBUG = FALSE, ...) {
  # Formerly called statDEModel()

  # TODO: curtis: IMPORT rugarch NAMESPACE -- Wed 29 Aug 2018 04:50:35 PM MDT
  args <- list(...)
  n <- length(dat)
  if (is.function(kn)) {
    k <- kn(n)
  } else {
    k <- kn
  }

  # Constants involved in the statistic
  bn <- (2 * log(log(n)) + k * log(log(log(n)))/2 - log(gamma(k/2)))^2/(
    2*log(log(n)))
  an <- sqrt(bn/(2*log(log(n))))

  args$data <- dat

  if (verbose) {
    cat("\n")
  }
  if (DEBUG) {
    DEBUG_list <- list()
  }
  # Get all values of the test statistic (will take max later)
  vals <- sapply(k:(n - k), function(t) {
    if (verbose) {
      cat(sprintf("\rChange points checked: %%%3.0f", (t - k)/(n - 2 * k) * 100))
    }
    if (use_theo) {
      args1 <- args
      args1$x <- dat
      args1$k <- t
      args1$return_coef <- TRUE
      res <- do.call(garch11_covmat_theo, args1)

      Ck <- res$C
      coef1 <- res$coef1
      coef2 <- res$coef2
    } else {
      # Split indices
      idx1 <- 1:t
      idx2 <- (t + 1):n

      # Parameter estimates around split
      args1 <- args
      args2 <- args
      args1$data <- dat[idx1]
      args2$data <- dat[idx2]
      mod1 <- do.call(ugarchfit, args = args1)
      mod2 <- do.call(ugarchfit, args = args2)
      coef1 <- coef(mod1)
      coef2 <- coef(mod2)
    }

    if (!use_theo) {
      H1 <- mod1@fit$hessian
      H2 <- mod2@fit$hessian

      V1 <- vcov(mod1, robust = TRUE) * t
      V2 <- vcov(mod2, robust = TRUE) * (n - t)

      GtG1 <- H1 %*% V1 %*% H1
      GtG2 <- H2 %*% V2 %*% H2

      H <- H1 + H2
      GtG <- GtG1 + GtG2

      Ck <- solve(H) %*% GtG %*% solve(H)
    }

    if (DEBUG) {
      DEBUG_list[[t]] <<- list("coef1" = coef1, "coef2" = coef2, "Ck" = Ck)
    }
    return(t(coef1 - coef2) %*% solve(Ck) %*% (coef1 - coef2))
  })
  if (verbose) {
    cat("\n")
  }

  # Get maximal value
  t <- k:(n - k)
  vals <- t * (n - t) * vals / n^2
  tau <- (k:(n - k))[which(vals == max(vals, na.rm = TRUE))]

  stat <- max((vals - bn)/an, na.rm = TRUE)

  res <- list("stat" = stat)
  if (estimate) {
    res$estimate <- tau
  }
  if (DEBUG) {
    res$quad_term <- (vals - bn)/an
    res$DEBUG_list <- DEBUG_list
  }

  res
}

################################################################################
# EXPERIMENTAL STATISTICAL TESTS
################################################################################

#' Rényi-Type Test for Structural Change for Regression Models
#'
#' Cleanly applies the implied test for the Rényi-type test for regression
#' models. This function acts as a clean interface to
#' \code{\link{stat_Zn_model_reg}}, and thus inherits the method's problems.
#'
#' @param x Data to test for structural change
#' @param formula The formula to represent the linear model
#' @param stat_plot Whether to create a plot of the values of the statistic at
#'        all potential change points
#' @inheritParams stat_Zn_model_reg
#' @return A \code{htest}-class object containing the results of the test
#' @examples
#' x <- rnorm(1000)
#' y <- 1 + 2 * x + rnorm(1000)
#' df <- data.frame(x, y)
#' HR.Reg.test(df, y ~ x)
HR.Reg.test <- function(x, formula, kn = log, stat_plot = FALSE,
                        wholevcov = TRUE, verbose = FALSE) {
  testobj <- list()
  testobj$method <-
    "Horvath-Rice Test for Change in Model Coefficents (Regression)"
  testobj$data.name <- deparse(substitute(x))

  res <- stat_Zn_model_reg(x,
                           kn = kn,
                           estimate = TRUE,
                           DEBUG = stat_plot,
                           formula = formula,
                           wholevcov = wholevcov,
                           verbose = verbose)
  stat <- res[[1]]
  est <- res[[2]]

  if (stat_plot) {
    series <- ts(res$quad_term, start = ceiling(kn(length(x))))
    plot.ts(res$quad_term, main = "Value of Test Statistic", ylab = "Statistic")
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

#' Rényi-Type Test for Structural Change for ARFIMA-GARCH Models
#'
#' Cleanly applies the implied test for the Rényi-type test for ARFIMA-GARCH
#' models. This function acts as a clean interface to
#' \code{\link{stat_Zn_model}}, and thus inherits the method's problems.
#'
#' TODO: curtis: THE DOCUMENTATION IS A LIE! -- Wed 29 Aug 2018 11:06:29 PM MDT
#'
#' @param x Data to test for structural change
#' @param spec Specification for the model; will be passed to
#'             \code{\link[fGarch]{garchFit}}
#' @inheritParams stat_Zn_model
#' @return A \code{htest}-class object containing the results of the test
#' @examples
#' HR.GARCH.test(rnorm(1000))
#' HR.GARCH.test(rnorm(1000), ~ garch(1, 1))
HR.GARCH.test <- function(x, spec, kn = log, stat_plot = FALSE,
                          wholevcov = TRUE, verbose = FALSE, ...) {
  testobj <- list()
  testobj$method <-
    "Horvath-Rice Test for Change in Model Coefficents (GARCH)"
  testobj$data.name <- deparse(substitute(x))

  args <- list(...)
  args$dat <- x
  args$spec <- spec
  args$kn <- kn
  args$estimate <- TRUE
  args$DEBUG <- stat_plot
  args$wholevcov <- wholevcov
  args$verbose <- verbose

  res <- do.call(stat_Zn_model, args)
  stat <- res$stat
  est <- res$estimate

  if (stat_plot) {
    series <- ts(res$quad_term, start = ceiling(kn(length(x))))
    plot.ts(res$quad_term, main = "Value of Test Statistic", ylab = "Statistic")
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

#' Ling's Darling-Erdös-Type Test for Structural Change for ARFIMA-GARCH Models
#'
#' Cleanly applies the implied test for \insertCite{ling07;textual}{CPAT}
#' Darling-Erdös-type test for ARFIMA-GARCH models. This function acts as a
#' clean interface to \code{\link{stat_de_model}}, and thus inherits the
#' method's problems.
#'
#' TODO: curtis: THE DOCUMENTATION IS A LIE! -- Wed 29 Aug 2018 11:06:29 PM MDT
#'
#' @param x Data to test for structural change
#' @param spec Specification for the model; will be passed to
#'             \code{\link[fGarch]{garchFit}}
#' @param stat_plot If \code{TRUE}, create a plot of the values of the statistic
#'                  at all potential change points
#' @inheritParams stat_de_model
#' @return A \code{htest}-class object containing the results of the test
#' @references
#'  \insertAllCited{}
#' @examples
#' DE.GARCH.test(rnorm(1000))
#' DE.GARCH.test(rnorm(1000), ~ garch(1, 1))
DE.GARCH.test <- function(x, spec, kn = log, stat_plot = FALSE,
                          wholevcov = TRUE, verbose = FALSE, ...) {
  testobj <- list()
  testobj$method <-
    "Ling Test for Change in Model Coefficents (GARCH)"
  testobj$data.name <- deparse(substitute(x))

  args <- list(...)
  args$dat <- x
  args$spec <- spec
  args$kn <- kn
  args$estimate <- TRUE
  args$DEBUG <- stat_plot
  args$wholevcov <- wholevcov
  args$verbose <- verbose

  res <- suppressWarnings(do.call(stat_de_model, args))
  stat <- res$stat
  est <- res$estimate

  if (stat_plot) {
    series <- ts(res$quad_term, start = ceiling(kn(length(x))))
    plot.ts(res$quad_term, main = "Value of Test Statistic", ylab = "Statistic")
  }

  kn_val <- kn(length(x))
  attr(kn_val, "names") <- deparse(substitute(kn)) %s0% "(T)"

  attr(stat, "names") <- "A"
  attr(est, "names") <- "t*"
  testobj$parameters <- kn_val
  testobj$p.value <- 1 - pdarling_erdos(res[[1]])
  testobj$estimate <- est
  testobj$statistic <- stat

  class(testobj) <- "htest"
  testobj
}

#' Ling's Darling-Erdös-Type Test for Structural Change for Regression Models
#'
#' Cleanly applies the implied test for \insertCite{ling07;textual}{CPAT}
#' Darling-Erdös-type test for regression models. This function acts as a clean
#' interface to \code{\link{stat_de_model_reg}}, and thus inherits the method's
#' problems.
#'
#' @param x Data to test for structural change
#' @param formula The formula to represent the linear model
#' @param stat_plot Whether to create a plot of the values of the statistic at
#'        all potential change points
#' @inheritParams stat_de_model_reg
#' @return A \code{htest}-class object containing the results of the test
#' @references
#'  \insertAllCited{}
#' @examples
#' x <- rnorm(1000)
#' y <- 1 + 2 * x + rnorm(1000)
#' df <- data.frame(x, y)
#' DE.Reg.test(df, y ~ x)
DE.Reg.test <- function(x, formula, kn = log, stat_plot = FALSE,
                        wholevcov = TRUE, verbose = FALSE) {
testobj <- list()
  testobj$method <-
    "Ling Test for Change in Model Coefficents (Regression)"
  testobj$data.name <- deparse(substitute(x))

  res <- suppressWarnings(stat_de_model_reg(x,
                                            kn = kn,
                                            estimate = TRUE,
                                            DEBUG = stat_plot,
                                            formula = formula,
                                            wholevcov = wholevcov,
                                            verbose = verbose))
  stat <- res[[1]]
  est <- res[[2]]

  if (stat_plot) {
    series <- ts(res$quad_term, start = ceiling(kn(length(x))))
    plot.ts(res$quad_term, main = "Value of Test Statistic", ylab = "Statistic")
  }

  kn_val <- kn(length(x))
  attr(kn_val, "names") <- deparse(substitute(kn)) %s0% "(T)"

  attr(stat, "names") <- "A"
  attr(est, "names") <- "t*"
  testobj$parameters <- kn_val
  testobj$p.value <- 1 - pdarling_erdos(res[[1]])
  testobj$estimate <- est
  testobj$statistic <- stat

  class(testobj) <- "htest"
  testobj
}
