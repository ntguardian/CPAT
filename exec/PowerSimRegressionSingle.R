#!/usr/bin/Rscript
################################################################################
# PowerSimRegressionSingle.R
################################################################################
# 2019-03-04
# Curtis Miller
################################################################################
# Computes simulated change point test statistics in the OLS linear regression
# context depending on how data should be simulated and what change point
# contexts are requested for a single scenario from contexts.
################################################################################

# optparse: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse")
  require("optparse")
}

################################################################################
# FUNCTIONS
################################################################################

`%s%` <- CPAT:::`%s%`
`%s0%` <- CPAT:::`%s0%`
stop_with_message <- CPAT:::stop_with_message
check_envir_has_objects <- CPAT:::check_envir_has_objects
base_file_name <- CPAT:::base_file_name

################################################################################
# MAIN FUNCTION DEFINITION
################################################################################

main <- function(SIMINPUT, CONTEXTINPUT, TESTINPUT, SIMINPUTPOST = "",
                 heterobreak = 1, seed = 20190219, seedless = FALSE,
                 replications = 5000, cores = NULL, stat = NULL, cpt = NULL,
                 n = NULL, case = NULL, alpha = 0.05, conflevel = .95,
                 help = FALSE) {
  # This function will be executed when the script is called from the command
  # line; the help parameter does nothing, but is needed for do.call() to work

  suppressPackageStartupMessages(library(doParallel))
  suppressPackageStartupMessages(library(doRNG))
  suppressPackageStartupMessages(library(purrr))
  suppressPackageStartupMessages(library(CPAT))

  if (is.null(cores)) {
    cores = max(1, detectCores() - 1)
  }
  registerDoParallel(cores = cores)
  if (seedless) {
    seed <- sample(1:99999999, 1)
  }
  set.seed(seed)

  if (SIMINPUTPOST == "") {
    SIMINPUTPOST <- SIMINPUT
  }

  simulation_tools <- new.env()
  simulation_post_tools <- new.env()
  context_tools <- new.env()
  test_tools <- new.env()
  load(SIMINPUT, envir = simulation_tools)
  load(SIMINPUTPOST, envir = simulation_post_tools)
  load(CONTEXTINPUT, envir = context_tools)
  load(TESTINPUT, envir = test_tools)
  simulation_tools_expected_objects <- c("df_generator", "eps_generator")
  context_tools_expected_objects <- c("n_values", "kstar_functions",
                                      "struc_models")
  test_tools_expected_objects <- c("stat_functions")
  check_envir_has_objects(simulation_tools_expected_objects,
                          envir = simulation_tools, blame_string = SIMINPUT)
  check_envir_has_objects(simulation_tools_expected_objects,
                          envir = simulation_post_tools,
                          blame_string = SIMINPUTPOST)
  check_envir_has_objects(context_tools_expected_objects, envir = context_tools,
                          blame_string = CONTEXTINPUT)
  check_envir_has_objects(test_tools_expected_objects, envir = test_tools,
                          blame_string = TESTINPUT)

  # It should be safe now to pull what we want from the files we loaded
  # But always check assumptions
  n_values <- context_tools$n_values
  stop_with_message(is.integer(n_values) & all(n_values > 0),
                    "Invalid n_values from" %s% CONTEXTINPUT %s0% "; must" %s%
                    "be positive integers")
  kstar_functions <- context_tools$kstar_functions
  stop_with_message(is.vector(kstar_functions) &
                    is.function(kstar_functions[[1]]) &
                      !any(is.null(names(kstar_functions))),
                    "Invalid kstar_functions from" %s% CONTEXTINPUT %s0% ";" %s%
                    "must be named vector of functions")
  struc_models <- context_tools$struc_models
  stop_with_message(is.list(struc_models) &
                    !any(is.null(names(struc_models))) &
                    all(sapply(struc_models, is.matrix)) &
                    all(sapply(struc_models, function(m) {ncol(m) == 2})) &
                    is.numeric(as.vector(sapply(struc_models, as.vector))),
                    "Invalid struc_models from" %s% CONTEXTINPUT %s0% "; it" %s%
                    "must be a named list of numeric matrices that have two" %s%
                    "columns")

  stat_functions <- test_tools$stat_functions
  stop_with_message(is.vector(stat_functions) &
                    is.function(stat_functions[[1]]) &
                    !any(is.null(names(stat_functions))), "Invalid" %s%
                    "stat_functions from" %s% TESTINPUT %s0% "; must be" %s%
                    "a named vector of functions")
  pval_functions <- test_tools$pval_functions
  stop_with_message(is.vector(pval_functions) &
                    is.function(pval_functions[[1]]) &
                    !any(is.null(names(pval_functions))) &
                    all(sort(names(pval_functions)) ==
                      sort(names(stat_functions))),
                    "Invalid pval_functions from" %s% TESTINPUT %s0% ";" %s%
                    "must be a named vector of functions with names" %s%
                    "agreeing with those in stat_functions (which could" %s%
                    "also be the source of this error)")
  
  df_generator <- simulation_tools$df_generator
  stop_with_message(is.function(df_generator) & all(c("n", "beta", "eps") %in%
                    names(formals(df_generator))), "Invalid df_generator" %s%
                    "from" %s% SIMINPUT %s0% "; must be a function with" %s%
                    "arguments n, beta, and eps")
  eps_generator <- simulation_tools$eps_generator
  stop_with_message(is.function(df_generator) & "n" %in%
                    names(formals(eps_generator)), "Invalid eps_generator" %s%
                    "from" %s% SIMINPUT %s0% "; must be a function with" %s%
                    "argument n")
  
  df_generator_post <- simulation_post_tools$df_generator
  stop_with_message(is.function(df_generator_post) & all(c("n", "beta",
                      "eps") %in% names(formals(df_generator_post))),
                    "Invalid df_generator from" %s% SIMINPUTPOST %s0%
                    "; must be a function with" %s%
                    "arguments n, beta, and eps")
  eps_generator_post <- simulation_post_tools$eps_generator
  stop_with_message(is.function(eps_generator_post) & "n" %in%
                    names(formals(eps_generator_post)), "Invalid" %s%
                    "eps_generator from" %s% SIMINPUTPOST %s0% "; must be a" %s%
                    "function with argument n")

  stop_with_message(heterobreak <= 1 & heterobreak >= 0, "heterobreak must" %s%
                    "be a number between 0 and 1")

  if (is.null(stat)) {
    stat <- names(stat_functions)[[1]]
  }
  if (is.null(cpt)) {
    cpt <- names(kstar_functions)[[1]]
  }
  if (is.null(n)) {
    n <- n_values[[1]]
  }
  if (is.null(case)) {
    case <- names(struc_models)[[1]]
  }

  stop_with_message(stat %in% names(stat_functions), "stat must identify" %s%
                    "one of the functions provided in stat_functions from" %s%
                    TESTINPUT)
  stop_with_message(cpt %in% names(kstar_functions), "cpt must identify" %s%
                    "one of the functions provided in kstar_functions from" %s%
                    CONTEXTINPUT)
  stop_with_message(is.integer(n) & n > 0, "n must be a positive integer")
  stop_with_message(case %in% names(struc_models), "case must identify" %s%
                    "one of the models in struc_models from" %s% CONTEXTINPUT)

  # TODO: curtis: THIS IS A CRUDE SOLUTION TO THE PROBLEM OF COMMUNICATING WHAT
  #               d IS TO THE P-VALUE FUNCTIONS; A MORE FLEXIBL, ROBUST, AND
  #               GENERAL SOLUTION IS NEEDED IN MORE COMPLEX SIMULATIONS.
  #               EFFECTIVELY ASSUMING THAT d IS CONSTANT, WHICH MAY NOT BE
  #               ALWAYS THE CASE -- Mon 15 Apr 2019 05:11:35 PM MDT
  d <- nrow(struc_models[[1]])

  # Collecting parameters, now that condition checking is done
  n_theta <- ceiling(n * heterobreak)
  kstar <- kstar_functions[[cpt]]
  k <- floor(kstar(n))
  stop_with_message(is.numeric(k) & k >= 1 & k <= n, "The function in" %s%
                    "kstar_functions called" %s% cpt %s% "did" %s%
                    "not return a number between 1 and" %s% n %s% "when" %s%
                    "expected; this should happen when" %s% cpt %s0%
                    "(n) is called")
  stat_name <- stat
  stat <- stat_functions[[stat]]
  stop_with_message(all(c("data", "formula") %in% names(formals(stat))),
                    "The function in stat_functions called" %s% 
                    stat_name %s% "does not have all of the" %s%
                    "expected arguments (data, formula)")
  r <- struc_models[[case]]
  beta1 <- r[, 1]
  beta2 <- r[, 2]
  d <- nrow(r)

  # Now the simulations can commence
  old_seed <- seed
  seed <- seed + which(names(struc_models) == case)[[1]] - 1  # [[1]] is safe
  simulation <- foreach(throwaway = 1:replications,
                        .options.RNG = seed, .combine = c) %dorng% {
    eps <- eps_generator(n = n)
    if (n_theta >= 1 & k >= 1) {
      df1 <- df_generator(n = min(n_theta, k), beta = beta1,
        eps = eps[1:min(n_theta, k)])
    } else {
      df1 <- data.frame()
    }

    if (n_theta == n) {
      df3 <- data.frame()
    } else {
      df3 <- df_generator_post(n = n - max(n_theta, k), beta = beta2,
                               eps = eps[(max(n_theta, k) + 1):n])
    }

    if (n_theta < k) {
      df2 <- df_generator_post(n = k - n_theta, beta = beta1,
                               eps = eps[(n_theta + 1):k])
    } else if (n_theta > k){
      df2 <- df_generator(n = n_theta - k, beta = beta2,
                          eps = eps[(k + 1):n_theta])
    } else {  # Equality
      df2 <- data.frame()
    }

    df <- rbind(df1, df2, df3)
    stop_with_message("y" %in% names(df), "y must be one of the" %s%
                      "columns generated by the data frame" %s%
                      "returned by df_generator(); it is the" %s%
                      "dependent variable")
    names(df[which(names(df) != "y")]) <- paste0("X", 1:(d - 1))
    f <- "y" %s% "~" %s% paste0("X", 1:(d - 1), collapse = " + ")
    f <- as.formula(f)
    
    stat(formula = f, data = df)
  }
  pvals <- sapply(simulation, function(z) pval_functions[[stat_name]](z, d = d))
  rejections <- pvals <= alpha

  # Reporting
  cat("RESULTS\n-------\n\n")
  cat("Simulation Definition File:", SIMINPUT %s0% "\n")
  cat("Test Statistic:", stat_name %s0% "\n")
  cat("Change Point Generator:", cpt %s0% "\n")
  cat("Change Point:", k %s0% "\n")
  cat("Change Case:", case %s0% "\n")
  cat("Sample Size:", n %s0% "\n")
  cat("Significance Level:", alpha %s0% "\n")
  cat("Simulations:", replications %s0% "\n")
  cat("Seed:", old_seed %s0% "\n")
  cat("Empirical Rejection Rate:", mean(rejections) %s0% "\n")
  ci <- binom.test(sum(rejections), length(rejections), conf.level = conflevel)
  ci <- ci$conf.int
  cat((100 * conflevel) %s0% "%", "Confidence Interval for Rejection Rate:",
       "(" %s0% ci[[1]] %s0% ", " %s0% ci[[2]] %s0% ")\n")
}

################################################################################
# INTERFACE SETUP
################################################################################

if (sys.nframe() == 0) {
  cl_args <- parse_args(OptionParser(
        description = "Computes simulated change point test statistics in" %s%
                      "OLS linear regression context for a single case.",
        option_list = list(
          make_option(c("--SIMINPUT", "-S"), type = "character",
                      help = "Name of .Rda file containing data defining" %s%
                             "how data simulations are to be performed"),
          make_option(c("--CONTEXTINPUT", "-C"), type = "character",
                      help = "Name of .Rda file defining what contexts" %s%
                             "in which to perform simulations"),
          make_option(c("--TESTINPUT", "-T"), type = "character",
                      help = "Name of .Rda file defining statistical tests" %s%
                             "that are simulated"),
          make_option(c("--SIMINPUTPOST", "-I"), type = "character",
                      default = "", help = "Name of .Rda file defining how" %s%
                                           "data simulations are to be" %s%
                                           "performed post-change in second" %s%
                                           "moment of data generating" %s%
                                           "process; if not set, ignored"),
          make_option(c("--heterobreak", "-e"), type = "numeric", default = 1,
                      help = "Location (as proportion of sample) at which" %s%
                             "to break from using SIMINPUT to generate data" %s%
                             "and use SIMINPUTPOST instead; if 1," %s%
                             "effectively ignored (always rounds up)"),
          make_option(c("--seed", "-s"), type = "integer", default = 20190219,
                      help = "The seed for simulations"),
          make_option(c("--seedless", "-R"), action = "store_true",
                      default = FALSE,
                      help = "Don't set a seed (causes --seed to be ignored)"),
          make_option(c("--replications", "-N"), type = "integer",
                      default = 5000,
                      help = "Number of replications per context"),
          make_option(c("--cores", "-c"), type = "integer", default = NULL,
                      help = "Number of cores to use in simulation; set to" %s%
                             "1 to effectively disable parallelization, or" %s%
                             "leave unset to use default (all but one core)"),
          make_option(c("--stat", "-t"), type = "character", default = NULL,
                      help = "The test statistic (from TESTINPUT) to" %s%
                             "simulate"),
          make_option(c("--cpt", "-p"), type = "character", default = NULL,
                      help = "The function to determine the change point" %s%
                             "(if any) in the data (from CONTEXTINPUT)"),
          make_option(c("--n", "-n"), type = "integer", default = NULL,
                      help = "The sample size to simulate"),
          make_option(c("--case", "-d"), type = "character", default = NULL,
                      help = "The structural change situation to simulate" %s%
                             "(from CONTEXTINPUT)"),
          make_option(c("--alpha", "-a"), type = "double", default = NULL,
                      help = "The level of significance at which to" %s%
                             "determine whether to reject the null hypothesis"),
          make_option(c("--conflevel", "-o"), type = "double", default = 0.95,
                      help = "The confidence level of the reported interval")
         )
      ))

  do.call(main, cl_args)
}

