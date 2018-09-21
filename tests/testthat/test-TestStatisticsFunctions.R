################################################################################
# TestStatisticsFunctions.R
################################################################################
# 2018-09-20
# Curtis Miller
################################################################################
# Test statistics functions for proper behavior.
################################################################################

################################################################################
# SCAFFOLDING
################################################################################

context("Statistics functions")
library(CPAT)

has_cointReg <- requireNamespace("cointReg", quietly = TRUE)

#' Check for cointReg
#'
#' Check for the cointReg package; if not present, skip the test.
#' 
#' @examples
#' check_cointReg()
check_cointReg <- function() {
  if (!has_cointReg) {
    skip("cointReg not available")
  }
}

#' Some Kernel Function
#'
#' A function to be used for testing kernel-based methods
#'
#' @param j A number to pass to the kernel function
#' @return Value of the kernel function
#' @examples
#' ker_func(0.5)
ker_func <- function(j) {
  max(1 - abs(j), 0)
}

# Dataset for testing
dat <- c(-0.64, -1.17, 0.76,  1.77, 1.77, -0.28, -0.20, 0.83, -0.08,  2.57, 
         -0.13, -1.22, 0.47, -0.32, 0.18, -1.38,  0.21, 0.38,  0.15, -0.31)

################################################################################
# UNDERLYING FUNCTION TESTING
################################################################################

test_that("stat_Zn() functions properly", {
  expect_equal(CPAT:::stat_Zn(dat), 1.58902205341715)

  res1 <- CPAT:::stat_Zn(dat, estimate = TRUE, get_all_vals = TRUE)
  expect_equal(res1$statistic,  1.58902205341715)
  expect_equal(res1$estimate, 10)
  expect_equal(res1$stat_vals,
               c(0.0303450116854491, 0.907107559783405, 0.584115126572313,
                 0.3722609420561530, 0.638146467227364, 0.514186903623801,
                 1.5890220534171500, 1.457098869814150, 0.845285623368565,
                 1.0378836334000500, 0.873366776528059, 0.987448725093092,
                 0.1530582795440690))

  expect_equal(CPAT:::stat_Zn(dat, use_kernel_var = TRUE, kernel = ker_func,
                              bandwidth = sqrt), 2.27610096950933)
  expect_equal(CPAT:::stat_Zn(dat, custom_var = function(x, k) {
                                var(x[1:(min(k + 1, length(x)))])
                              }), 1.25197565798885)
  expect_equal(CPAT:::stat_Zn(dat, custom_var = function(x, k) var(x)),
               1.439370888561)
})

test_that("get_lrv_vec() functions properly, interfaces where it should", {
  expect_equal(CPAT:::get_lrv_vec(dat, kernel = ker_func, bandwidth = sqrt),
               c(0.581238887175265, 0.632229813007764, 0.681431014190319,
                 0.646220314031592, 0.611266882865371, 0.578769643980172,
                 0.582735929700783, 0.533450929616323, 0.522333689063721,
                 0.411455302969274, 0.382884504978510, 0.438867732426083,
                 0.456809375948772, 0.518036697090424, 0.559787452269136,
                 0.630613678319110, 0.633485003925464, 0.633403884907476,
                 0.636668741508810))

  check_cointReg()

  expect_equal(CPAT:::get_lrv_vec(dat),
               c(0.944000055543033, 0.854345524470940, 0.938887787065566,
                 0.988543921777035, 0.950674807065429, 0.966696492122214,
                 0.977280759919950, 0.962346545793554, 0.963543535868417,
                 0.844971457475420, 0.858397192435424, 0.934364474783002,
                 0.924364070223716, 0.941827034624808, 0.944039246554919,
                 0.986652076557760, 0.986822029414594, 0.981427647381865,
                 0.976271063827923))
})
