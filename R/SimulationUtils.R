################################################################################
# SimulationUtils.R
################################################################################
# 2018-08-27
# Curtis Miller
################################################################################
# Utilities for performing and managing statistical tests and simulations used
# in our paper (Horváth et. al. 2019).
################################################################################

################################################################################
# POWER SIMULATION OBJECT AND FILE MANAGEMENT
################################################################################

#' Convert Rényi-Type Statistic Power Simulation Save List to Data Frame
#'
#' This function will convert the power simulation data generated in a list in
#' our simulation scripts to a \code{data.frame}. Given such a list and a
#' critical value to determine whether the null hypothesis should be rejected,
#' the function will return a \code{data.frame} with columns \code{power},
#' \code{stat}, \code{dist}, \code{kn}, \code{n}, \code{cpt}, and \code{delta},
#' which correspond to: the empirical power of the statistic; the identifier of
#' the statistic; the generating distribution of the statistic was computed on;
#' the \code{kn} parameter; the identifier of how change points were computed;
#' and the size of the change.
#'
#' @param obj A list containing simulated statistic values
#' @param crit The critical value determining whether a statistic should lead to
#'        the rejection of the null hypothesis
#' @return A \code{data.frame} summarizing the results of the data stored in
#'         \code{obj}
#' @examples
#' saveobj <- list("norm" = list(
#'                   "log" = list(
#'                     "n50" = list(
#'                       "c4rt" = list(
#'                         "d_0" = c(1.551, 1.276, 1.348, 1.982, 1.423)
#' )))))
#' power_sim_Zn_to_df(saveobj, qZn(.95))
power_sim_Zn_to_df <- function(obj, crit) {
  # Formerly known as powerSim_Zn_to_df()

  return_obj <- list(power = c(), stat = c(), dist = c(),
                     kn = c(), n = c(), cpt = c(),
                     delta = c())
  for (dist in names(obj)) {
    for (kn in names(obj[[dist]])) {
      for (n in names(obj[[dist]][[kn]])) {
        for (cpt in names(obj[[dist]][[kn]][[n]])) {
          newprops <- sapply(names(
            obj[[dist]][[kn]][[n]][[cpt]]),
            function(delta) {
              return(mean(
                obj[[dist]][[kn]][[n]][[cpt]][[
                  delta]] >= crit))})
          return_obj$power <- c(return_obj$power,
                                newprops)
          return_obj$dist <- c(return_obj$dist,
                               rep(dist,
                                   times = length(newprops)))
          return_obj$kn <- c(return_obj$kn,
                             rep(kn,times =
                                   length(newprops)))
          return_obj$n <- c(return_obj$n,
                            rep(n, times = length(newprops)))
          return_obj$cpt <- c(return_obj$cpt,
                              rep(cpt,
                                  times = length(newprops)))
          return_obj$delta <- c(return_obj$delta,
                                names(
                                  obj[[dist]][[kn]][[n]][[
                                    cpt]]))
        }
      }
    }
  }
  return_obj$stat <- rep("Zn",
                         times = length(return_obj$power))

  as.data.frame(return_obj)
}

#' Convert CUSUM-Type Statistic Power Simulation Save List to Data Frame
#'
#' This function will convert the power simulation data generated in a list in
#' our simulation scripts to a \code{data.frame}. Given such a list and a
#' critical value to determine whether the null hypothesis should be rejected,
#' the function will return a \code{data.frame} with columns \code{power},
#' \code{stat}, \code{dist}, \code{n}, \code{cpt}, and \code{delta},
#' which correspond to: the empirical power of the statistic; the identifier of
#' the statistic; the generating distribution of the statistic was computed on;
#' the identifier of how change points were computed; and the size of the
#' change.
#'
#' @inheritParams power_sim_Zn_to_df
#' @return A \code{data.frame} summarizing the results of the data stored in
#'         \code{obj}
#' @examples
#' saveobj <- list("norm" = list(
#'                   "n50" = list(
#'                     "c4rt" = list(
#'                       "d_0" = c(1.551, 1.276, 1.348, 1.982, 1.423)
#' ))))
#' power_sim_Vn_to_df(saveobj, qkolmogorov(.95))
power_sim_Vn_to_df <- function(obj, crit) {
  # Formerly known as powerSim_Vn_to_df()

  return_obj <- list(power = c(), stat = c(), dist = c(),
                     n = c(), cpt = c(), delta = c())
  for (dist in names(obj)) {
    for (n in names(obj[[dist]])) {
      for (cpt in names(obj[[dist]][[n]])) {
        newprops <- sapply(names(obj[[dist]][[n]][[cpt]]),
                           function(delta) {return(
                             mean(obj[[dist]][[n]][[cpt]][[
                               delta]] >= crit))})
        return_obj$power <- c(return_obj$power, newprops)
        return_obj$dist <- c(return_obj$dist,
                             rep(dist,
                                 times = length(newprops)))
        return_obj$n <- c(return_obj$n,
                          rep(n, times = length(newprops)))
        return_obj$cpt <- c(return_obj$cpt,
                            rep(cpt,
                                times = length(newprops)))
        return_obj$delta <- c(return_obj$delta,
                              names(obj[[dist]][[n]][[cpt]]))
      }
    }
  }
  return_obj$stat <- rep("Vn",
                         times = length(return_obj$power))

  as.data.frame(return_obj)
}

#' Create Power Simulation Results Data Frame
#'
#' Creates a \code{data.frame} that contains power simulation results from files
#' containing power simulations. This function should automate the use of
#' \code{\link{power_sim_Zn_to_df}} and \code{\link{power_sim_Vn_to_df}} for
#' collecting power simulation data. It takes two \code{CSV} files, one passed
#' (as a character string) to \code{file_meta} and the other to
#' \code{stat_meta}, describing how the files (named and described in
#' \code{file_meta}) should be handled.
#'
#' @param file_meta The location of a \code{CSV} file that contains file names
#'                  and the statistics that those files correspond to
#' @param stat_meta The location of a \code{CSV} file that contains statistic
#'                  (\code{stat}) labels (used in \code{file_meta}). the name
#'                  of the variable for the statistic, and the name of the
#'                  function that converts a file (mentioned in
#'                  \code{file_meta}) to a \code{data.frame} of power data
#' @param prefix Character string representing a prefix for file names mentioned
#'               in \code{file_meta}; could be used for adding path information
#'               to those names, in case the files are not in the working
#'               directory and there is no desire to edit \code{file_meta}'s
#'               data
#' @return A data frame containing the power simulation data
#' @import utils
#' @examples
#' \dontrun{
#' power_sim_stat_df_creator("FileStatMeta.csv", "StatMeta.csv")
#' }
power_sim_stat_df_creator <- function(file_meta, stat_meta, prefix = "",
                                      alpha = 0.05) {
  stat_meta <- read.csv(stat_meta, stringsAsFactors = FALSE)
  file_meta <- read.csv(file_meta, stringsAsFactors = FALSE)
  file_meta$file <- paste0(prefix, file_meta$file)

  stats <- unique(file_meta$statistic)

  perc_Zn_theo <- qZn(1 - alpha)
  perc_Vn_theo <- qkolmogorov(1 - alpha)
  perc_de_theo <- qdarling_erdos(1 - alpha)
  perc_hs_theo <- qhidalgo_seo(1 - alpha)

  Zn_df <- bind_power_sim_objs(with(file_meta, file[statistic == "Zn"]),
                               with(stat_meta, eval(as.name(
                                 crit_value[statistic == "Zn"]))),
                               with(stat_meta, eval(as.name(
                                 conv_func[statistic == "Zn"]))),
                               "Zn")

  df <- Zn_df

  kn <- unique(Zn_df$kn)
  for (s in stats[which(stats != "Zn")]) {
    temp <- bind_power_sim_objs(with(file_meta, file[statistic == s]),
                                with(stat_meta, eval(as.name(
                                  crit_value[statistic == s]))),
                                with(stat_meta, eval(as.name(
                                  conv_func[statistic == s]))),
                                s)
    for (k in kn) {
      temp$kn <- k
      df <- rbind(df, temp)
    }
  }

  df$delta <- as.character(df$delta)
  df$delta <- as.numeric(substr(df$delta, 3, nchar(df$delta)))

  df$n <- as.character(df$n)
  df$n <- as.integer(substr(df$n, 2, nchar(df$n)))

  df[with(df, order(stat, dist, kn, n, cpt, delta)),]
}

#' Power Result Data Frame Creation
#'
#' Creates a \code{data.frame} containing power simulation results. Effectively
#' a better, higher-level interface to \code{\link{power_sim_Zn_to_df}} and
#' \code{\link{power_sim_Vn_to_df}}.
#'
#' @param files A character vector of file names
#' @param crit_value The critical value against which to compare a test
#'                   statistic
#' @param conv_func The function responsible for converting a list containing
#'                  simulated statistic values under different conditions to a
#'                  \code{data.frame}
#' @param stat_name The label of the statistic
#' @return A \code{data.frame} containing power levels
#' @examples
#' \dontrun{
#' filenames <- c("powerSimulations_sdest_norm_DE.rda",
#'                "powerSimulations_sdest_ar1_0.5_DE.rda")
#' bind_power_sim_objs(filenames, crit_value = qdarling_erdos(.95),
#'                     conv_func = power_sim_Vn_to_df, stat_name = "de")
#' }
bind_power_sim_objs <- function(files, crit_value, conv_func, stat_name) {
  # Formerly known as bind_powerSim_objs()

  load(files[1])
  df <- conv_func(saveobj, crit_value)
  df$stat <- stat_name
  # print(head(df))
  # cat("\n\n")

  for (f in files[-1]) {
    load(f)
    temp <- conv_func(saveobj, crit_value)
    temp$stat <- stat_name
    # print(head(temp))
    # cat("\n\n")
    df <- rbind(df, temp)
    # print(tail(df))
    # cat("\n\n")
  }

  df
}

################################################################################
# CHANGING WINDOWN P-VALUE COMPUTATION
################################################################################

#' Expanding Window p-Values
#'
#' Gets p-values for the CUSUM, Darling-Erdös, Hidalgo-Seo, Andrews, and
#' Rényi-type tests when applied to an expanding window of data.
#'
#' @param dat The dataset for which to test for change in mean
#' @param m The location of the first potential change point for Andrews' test
#' @return A matrix containing p-values for an expanding sample size, with each
#'         row corresponding to one observation larger; columns are labeled for
#'         each statistic
#' @examples
#' if (require("foreach") & require("doParallel")) {
#'   get_expanding_window_pvals(rnorm(1000), m = 900)
#' }
get_expanding_window_pvals <- function(dat, m = Inf) {
  has_parallel <- requireNamespace("foreach", quietly = TRUE) &&
                  requireNamespace("doParallel", quietly = TRUE)
  if (!has_parallel) {stop("The foreach and doParallel packages are needed" %s%
                           "to use this function")}
  `%dopar%` <- foreach::`%dopar%`
  foreach <- foreach::foreach
  test_p_vals <- foreach(n = 3:length(dat), .combine = rbind) %dopar% {
      cusum <- CUSUM.test(dat[1:n])
      de <- DE.test(dat[1:n])
      hr <- HR.test(dat[1:n])
      hs <- HS.test(dat[1:n], corr = FALSE)
  
      cusum_ker <- CUSUM.test(dat[1:n], use_kernel_var = TRUE)
      de_ker <- DE.test(dat[1:n], use_kernel_var = TRUE)
      hr_ker <- HR.test(dat[1:n], use_kernel_var = TRUE)
      hs_ker <- HS.test(dat[1:n], corr = TRUE)

      if (n > m) {
        andrew <- Andrews.test(dat[1:n], m = m)$p.value[[1]]
      } else {
        andrew <- NA
      }
      
  
      c("n" = n,
        "cusum" = cusum$p.value[[1]],
        "de" = de$p.value[[1]],
        "hr" = hr$p.value[[1]],
        "hs" = hs$p.value[[1]],
        "cusum_ker" = cusum_ker$p.value[[1]],
        "de_ker" = de_ker$p.value[[1]],
        "hr_ker" = hr_ker$p.value[[1]],
        "hs_ker" = hs_ker$p.value[[1]],
        "andrew" = andrew)
  }
  
  rownames(test_p_vals) <- NULL
  test_p_vals
}

#' Expanding Window p-Values for Regression Models
#'
#' Gets p-values for the CUSUM, Darling-Erdös, Hidalgo-Seo, Andrews, and
#' Rényi-type tests when applied to an expanding window of data for a regression
#' model.
#'
#' @param data A \code{data.frame}, the dataset for which to test for
#'                structural change
#' @param formula The regression model formula, which will be passed to
#'                \code{\link[stats]{lm}}
#' @param min_n An integer; the minimum sample size
#' @param m The location of the first potential change point for Andrews' test
#' @param verbose If \code{TRUE}, send messages to output
#' @return A matrix containing p-values for an expanding sample size, with each
#'         row corresponding to one observation larger; columns are labeled for
#'         each statistic
#' @examples
#' x <- rnorm(1000)
#' y <- 1 + 2 * x + rnorm(1000)
#' df <- data.frame(x, y)
#' if (require("foreach") & require("doParallel")) {
#'   get_expanding_window_pvals_reg(y ~ x, data = df, min_n = 4, m = 900)
#' }
get_expanding_window_pvals_reg <- function(formula, data, min_n = 3, m = Inf,
                                           verbose = FALSE) {
  has_parallel <- requireNamespace("foreach", quietly = TRUE) &&
                  requireNamespace("doParallel", quietly = TRUE)
  if (!has_parallel) {stop("The foreach and doParallel packages are needed" %s%
                           "to use this function")}
  `%dopar%` <- foreach::`%dopar%`
  foreach <- foreach::foreach

  if (requireNamespace("cointReg", quietly = TRUE)) {
    pkgs <- c("stats", "cointReg")
  } else {
    pkgs <- c("stats")
  }

  test_p_vals <- foreach(n = min_n:nrow(data), .combine = rbind,
                         .packages = pkgs) %dopar% {
    if (verbose) {
      cat("Working on n = ", n, "out of ", nrow(data), "\n")
    }
    vec <- residuals(lm(formula = formula, data = data))
    cusum <- CUSUM.test(vec[1:n])
    de <- DE.test(vec[1:n])
    hr <- HR.test(vec[1:n])
    hs <- HS.test(vec[1:n], corr = FALSE)
    
    cusum_ker <- CUSUM.test(vec[1:n], use_kernel_var = TRUE)
    de_ker <- DE.test(vec[1:n], use_kernel_var = TRUE)
    hr_ker <- HR.test(vec[1:n], use_kernel_var = TRUE)
    hs_ker <- HS.test(vec[1:n], corr = TRUE)
    
    if (n > m) {
      andrew <- tryCatch(Andrews.test(x = data[1:n,], formula = formula,
                                      m = m)$p.value[[1]],
                         error = function(e) {return(NA)})
    } else {
      andrew <- NA
    }
    
    c("n" = n,
      "cusum" = cusum$p.value[[1]],
      "de" = de$p.value[[1]],
      "hr" = hr$p.value[[1]],
      "hs" = hs$p.value[[1]],
      "cusum_ker" = cusum_ker$p.value[[1]],
      "de_ker" = de_ker$p.value[[1]],
      "hr_ker" = hr_ker$p.value[[1]],
      "hs_ker" = hs_ker$p.value[[1]],
      "andrew" = andrew)
  }
  
  rownames(test_p_vals) <- rownames(data[test_p_vals[, "n"],])
  test_p_vals
}
