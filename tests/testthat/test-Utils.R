################################################################################
# test-Utils.R
################################################################################
# 2019-02-19
# Curtis Miller
################################################################################
# Check that certain utility functions work correctly.
################################################################################

################################################################################
# SCAFFOLDING
################################################################################

context("Utility functions")
library(CPAT)

################################################################################
# USEFUL OBJECTS
################################################################################

`%s0%` <- CPAT:::`%s0%`
`%s%` <- CPAT:::`%s%`

x <- 1

test_env <- new.env()
test_env$y <- 1

################################################################################
# OPERATORS TESTS
################################################################################

test_that("Concatenation functions work", {
  expect_equal("Hello" %s% "world!", "Hello world!")
  expect_equal("Hello" %s0% "world!", "Helloworld!")
})

################################################################################
# CLASS CHECKING FUNCTIONS TEST
################################################################################

test_that("is-type functions work", {
  expect_true(CPAT:::is.formula(y ~ x))
  expect_false(CPAT:::is.formula("Hello"))
})

################################################################################
# ERROR CHECKING FUNCTIONS TESTS
################################################################################

test_that("Environment checking function works", {
  expect_silent(CPAT:::check_envir_has_objects(c("y"), envir = test_env))
  expect_error(CPAT:::check_envir_has_objects(c("x", "y"), envir = test_env),
               "test_env does not have all expected objects; must have x, y")
  expect_error(CPAT:::check_envir_has_objects(c("x"), envir = test_env,
                                              blame_string = "TEST"),
               "TEST does not have all expected objects; must have x")
})

test_that("stop_with_message() works", {
  expect_silent(CPAT:::stop_with_message(x == 1))
  expect_silent(CPAT:::stop_with_message(c(x == 1, TRUE)))
  expect_error(CPAT:::stop_with_message(c(x == 2, TRUE)))
  expect_error(CPAT:::stop_with_message(x == 2, message = "Hello!"), "Hello!")
})

################################################################################
# FILESYSTEM UTILITIES TESTS
################################################################################

test_that("base_file_name() works properly", {
  expect_equal(CPAT:::base_file_name("./some/tree/test.txt"), "test")
})
