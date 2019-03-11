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

#' Command-Line Utility for Printing Through Loop
#'
#' This is a tool for moving through an iterable object in the R command line.
#' \code{obj} is the object being iterated through, and \code{accessor} is a
#' function that can take this object and some index and return the item in that
#' (integer) index of the object. This is then printed to the screen until
#' stopped.
#'
#' When the utility is called, it enters a loop. The user can type the index of
#' the object she wishes to see, or \code{j} for the previous entry, or \code{k}
#' for the next entry. Pressing \code{q} ends the loop.
#'
#' @param obj The object through which to iterate
#' @param accessor A function with parameters \code{index} and \code{obj} that
#'                 will access the entry in position \code{index} of \code{obj}
#'                 and return the result
#' @export
#' @examples
#' \dontrun{
#' loop_print_iterator(1:10, function(index, obj) {obj[[index]]})
#' }
loop_print_iterator <- function(obj, accessor) {
  stopifnot(is.function(accessor))
  stopifnot(all(c("index", "obj") %in% names(formals(accessor))))
  index <- 1
  input <- ""
  while (TRUE) {
    system("clear")
    m <- 1
    if (input != "") {
      if (nchar(input) > 1 & suppressWarnings(is.na(as.integer(input)))) {
        n <- nchar(input)
        prefix <- substr(input, 1, n - 1)
        input <- substr(input, n, n)
        prefix <- suppressWarnings(as.integer(prefix))
        if (is.na(prefix) | input == "q") {
          input <- "z"  # A value meant to cause a later error
        } else {
          m <- prefix
        }
      }
      if (input == "j") {
        index <- max(index - m, 1)
      } else if (input == "k") {
        index <- index + m
	    } else if (input == "q") {
        break
      } else {
        index <- suppressWarnings(as.integer(input))
      }

      cat("INDEX:", index, "\n\n")
      x <- try(accessor(index = index, obj = obj), silent = TRUE)
      if (is(x, "try-error") | is.na(index)) {
        cat("Bad index! Could not access element. Index reset to 1.\n")
        index <- 1
      } else {
        print(x)
      }
    }
    input <- readline("(index, j, k, q)> ")
  }
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

