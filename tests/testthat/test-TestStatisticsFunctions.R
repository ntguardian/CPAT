################################################################################
# test-TestStatisticsFunctions.R
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

`%s%` <- CPAT:::`%s%`

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

dat2 <-  c(-0.73, -0.87,  0.5, 2.34, 1.41,  0.42,  0.06, 0.86,  0.89, 1.04, 
            1.48,  0.65, 1.84, 0.95, 0.94, -0.02, -0.08, 0.15, -0.07, 0.81)

df <- data.frame(x = 1:20, y = 1 + 2 * 1:20 + dat)
df2 <- data.frame(x = dat2, y = 1 + 2 * dat2 + dat)
X <- model.matrix(y ~ x, data = df2)
eX <- X * residuals(lm(y ~ x, data = df2))

################################################################################
# UNDERLYING FUNCTION TESTING
################################################################################

test_that("stat_Zn() functions properly", {
  expect_equal(CPAT:::stat_Zn(dat), 1.58902205341715)

  res <- CPAT:::stat_Zn(dat, estimate = TRUE, get_all_vals = TRUE)
  expect_equal(res$statistic,  1.58902205341715)
  expect_equal(res$estimate, 10)
  expect_equal(res$stat_vals,
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

  check_cointReg()

  expect_equal(CPAT:::stat_Zn(dat, use_kernel_var = TRUE), 1.58829744129333)
})

test_that("stat_Zn_reg() functions properly", {
  expect_error(CPAT:::stat_Zn_reg(dat), "Bad formula passed")
  expect_equal(CPAT:::stat_Zn_reg(y ~ x, data = df), 5.93426034417689)
  expect_equal(CPAT:::stat_Zn_reg(y ~ x, data = df,
                                  custom_var = function(x, k) {diag(2)}),
               49.1565680600069)
  expect_equal(CPAT:::stat_Zn_reg(y ~ x, data = df,
                                  custom_var = function(x, k) {diag(2)},
                                  estimate = TRUE)$estimate, 11L)
  expect_equal(CPAT:::stat_Zn_reg(y ~ x, data = df,
                                  custom_var = function(x, k) {diag(2)},
                                  get_all_vals = TRUE)$stat_vals,
               c(1.49379676274333, 3.84317040943000, 0.727071394520901, 
                 3.14210359828251, 1.10739365339191, 4.648741703797920,
                 15.3837927214878, 49.1565680600069, 16.00624104013680,
                 18.8677456895055, 10.1824122337828, 10.99171939307290,
                 10.5272602747702))
  expect_is(CPAT:::stat_Zn_reg(y ~ x, data = df,
                               custom_var = function(x, k) {diag(2)},
                               get_all_vals = TRUE, estimate = TRUE), "list")
  expect_equal(names(CPAT:::stat_Zn_reg(y ~ x, data = df,
                                        custom_var = function(x, k) {diag(2)},
                                        get_all_vals = TRUE, estimate = TRUE)),
               c("statistic", "estimate", "stat_vals"))

  check_cointReg()

  expect_equal(CPAT:::stat_Zn_reg(y ~ x, data = df, use_kernel_var = TRUE),
               6.08134617119088)
})

test_that("stat_Vn() functions properly", {
  expect_equal(CPAT:::stat_Vn(dat), 0.888290332296761)

  res <- CPAT:::stat_Vn(dat, estimate = TRUE, get_all_vals = TRUE)
  expect_equal(res$statistic, 0.888290332296761)
  expect_equal(res$estimate, 10)
  expect_equal(res$stat_vals, 
               c(0.1860482588333920, 0.520634431435331, 0.360310650616115,
                 0.0108565614250703, 0.380316406233684, 0.274285437237333,
                 0.1893711755857890, 0.342465331277215, 0.284564825236631,
                 0.8882903322967610, 0.806397600403535, 0.453627867421295,
                 0.5279770762213650, 0.410110571137388, 0.414000463800686,
                 0.0547597948127525, 0.064285500675461, 0.112576638143350,
                 0.1087825244876340))

  expect_equal(CPAT:::stat_Vn(dat, use_kernel_var = TRUE, kernel = ker_func,
                              bandwidth = sqrt), 1.27237912286901)
  expect_equal(CPAT:::stat_Vn(dat, custom_var = function(x, k) {
      var(x[1:(min(k + 1, length(x)))])
                              }), 0.699875669359526)
  expect_equal(CPAT:::stat_Vn(dat, custom_var = function(x, k) var(x)),
               0.804632787914169)

  check_cointReg()

  expect_equal(CPAT:::stat_Vn(dat, use_kernel_var = TRUE), 0.887885261805216)
})

test_that("stat_de() functions properly", {
  expect_equal(CPAT:::stat_de(dat), 1.58953572310406)

  res <- CPAT:::stat_de(dat, estimate = TRUE, get_all_vals = TRUE)
  expect_equal(res$statistic, 1.58953572310406)
  expect_equal(res$estimate, 10)
  expect_equal(res$stat_vals, 
               c( 0.6666030152225960,  1.54840316329497000, 0.822026714752378,
                 -0.1599035379267860,  0.69125817657553300, 0.411494952874413,
                  0.2099849222681120,  0.51200948869585500, 0.384951879021424,
                  1.5895357231040600,  1.43387522359022000, 0.738919065418096,
                  0.9198957847897980,  0.70788969448132500, 0.769048175390362,
                 -0.0501454544575806, -0.00700952891422735, 0.188210518988372,
                  0.3120835122629970))

  expect_equal(CPAT:::stat_de(dat, use_kernel_var = TRUE, kernel = ker_func,
                              bandwidth = sqrt), 2.35771330424856)
  expect_equal(CPAT:::stat_de(dat, custom_var = function(x, k) {
      var(x[1:(min(k + 1, length(x)))])
                              }), 2.02497320978937)
  expect_equal(CPAT:::stat_de(dat, custom_var = function(x, k) var(x)),
               1.42222063433888)

  check_cointReg()

  expect_equal(CPAT:::stat_de(dat, use_kernel_var = TRUE), 1.58872558212097)
})

test_that("stat_hs() functions properly", {
  expect_equal(CPAT:::stat_hs(dat), 1.21370945972696)

  res <- CPAT:::stat_hs(dat, estimate = TRUE, get_all_vals = TRUE)
  expect_equal(res$statistic, 1.21370945972696)
  expect_equal(res$estimate, 10)
  expect_equal(res$stat_vals, 
               c(-1.370582690156750,  1.07517939544294, -1.03102406018458, 
                 -2.267843009815840, -1.31990422541547, -1.81917583478302,
                 -2.068984134345630, -1.65939233959489, -1.85754100603907,
                  1.213709459726960,  0.69794177537658, -1.21847690153309,
                 -0.793814665362957, -1.28502483213330, -1.15205095731958,
                 -2.244864368803760, -2.22744398880951, -2.09014950652634,
                 -1.954438886999530))

  expect_equal(CPAT:::stat_hs(dat, use_kernel_var = TRUE, kernel = ker_func,
      bandwidth = sqrt), 3.99463145080931)
  expect_equal(CPAT:::stat_hs(dat, custom_var = function(x, k) {
                 var(x[1:(min(k + 1, length(x)))])
               }), 2.46376948017636)
  expect_equal(CPAT:::stat_hs(dat, custom_var = function(x, k) var(x)),
               0.236019124621259)
  expect_equal(CPAT:::stat_hs(dat, corr = TRUE), 1.21370945972696)

  check_cointReg()

  expect_equal(CPAT:::stat_hs(dat, use_kernel_var = TRUE), 0.78115982765937)
})

test_that("stat_hs_reg() functions properly", {
  expect_error(CPAT:::stat_hs_reg(dat), "Bad formula passed")
  expect_equal(CPAT:::stat_hs_reg(y ~ x, data = df), 4.31295336184296)
  expect_equal(CPAT:::stat_hs_reg(y ~ x, data = df, m = 6), 1.90718155854372)
  expect_equal(CPAT:::stat_hs_reg(y ~ x, data = df, m = function(n) {n/2}),
               3.8253699534694)
  expect_equal(CPAT:::stat_hs_reg(y ~ x, data = df, estimate = TRUE,
                                  get_all_vals = TRUE),
               list(statistic = 4.31295336184296, estimate = 1L,
                    stat_vals = c( 4.312953361842960, -0.464210620180501,
                                  -1.209761372351950, -1.511710762774970,
                                  -1.636033980821910, -1.423330875070440,
                                  -1.664290280254540,  3.863475451748580,
                                   2.729950006095740, -1.453799593143130,
                                  -0.422480788176892, -1.476981256289980,
                                  -0.937368463192843, -0.636798519922383,
                                  -1.379690943378280, -2.284485815966370)))
})

test_that("andrews_test() functions properly", {
  expect_error(CPAT:::andrews_test(dat), "argument \"M\" is missing")
  expect_equal(CPAT:::andrews_test(dat, 15),
               list(pval = 1, stat = 2.41061543748607))
  expect_equal(CPAT:::andrews_test(c(dat, c(1, 3, 2, 3, 2)), 20),
               list(pval = 0.3125, stat = 7.94772816422427))
  expect_equal(CPAT:::andrews_test(c(dat, c(1, 3, 2, 3, 2)), 20, pval = FALSE),
               7.94772816422427)
  expect_equal(CPAT:::andrews_test(c(dat, c(1, 3, 2, 3, 2)), 20, stat = FALSE),
               0.3125)
})

test_that("andrews_test_reg() functions properly", {
  expect_error(CPAT:::andrews_test_reg(dat), "Bad formula passed")
  expect_error(CPAT:::andrews_test_reg(y ~ x, data = df), "argument .* missing")
  expect_equal(CPAT:::andrews_test_reg(y ~ x, data = df, M = 15),
               list(pval = 0.545454545454545, stat = 0.578318267262971))
  expect_equal(CPAT:::andrews_test_reg(y ~ x, data = df, M = 15, pval = FALSE),
               0.578318267262971)
  expect_equal(CPAT:::andrews_test_reg(y ~ x, data = df, M = 15, stat = FALSE),
               0.545454545454545)
})

################################################################################
# LRV TEST
################################################################################

test_that("get_lrv_vec_old1() functions properly, interfaces where it should", {
  expect_equal(CPAT:::get_lrv_vec_old1(dat, kernel = ker_func, bandwidth = sqrt),
               c(0.658500239630696, 0.805352909354794, 0.807664511365906,
                 0.697010305338316, 0.653237233734026, 0.593204166349428,
                 0.582917880096948, 0.487181301225369, 0.486163002545009,
                 0.327311984831554, 0.286361056170444, 0.365052316130534,
                 0.379029670077905, 0.467345014456849, 0.530268017910291,
                 0.668028324699801, 0.670453151630940, 0.664687710663636,
                 0.674880875114663))


  check_cointReg()

  expect_equal(CPAT:::get_lrv_vec_old1(dat),
               c(0.944000055543033, 0.854345524470940, 0.938887787065566,
                 0.988543921777035, 0.950674807065429, 0.966696492122214,
                 0.977280759919950, 0.962346545793554, 0.963543535868417,
                 0.844971457475420, 0.858397192435424, 0.934364474783002,
                 0.924364070223716, 0.941827034624808, 0.944039246554919,
                 0.986652076557760, 0.986822029414594, 0.981427647381865,
                 0.976271063827923))
})

test_that("get_lrv_vec() functions properly, interfaces where it should", {
  expect_equal(CPAT:::get_lrv_vec(dat, kernel = ker_func, bandwidth = sqrt),
               c(0.658500239630696, 0.805352909354794, 0.807664511365906,
                 0.697010305338316, 0.653237233734026, 0.593204166349428,
                 0.582917880096948, 0.487181301225369, 0.486163002545009,
                 0.327311984831554, 0.286361056170444, 0.365052316130534,
                 0.379029670077905, 0.467345014456849, 0.530268017910291,
                 0.668028324699801, 0.670453151630940, 0.664687710663636,
                 0.674880875114663))


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

test_that("get_lrv_arr() functions properly, interfaces where it should", {
  expect_equal(CPAT:::get_lrv_arr(eX, kernel = ker_func,
                                  bandwidth = (1:20)^0.5), 
                                  structure(
                                    c( 0.0377436292675294, -0.02755284936529650, 
                                      -0.0275528493652965,  0.02011358003666640,
                                       0.0194323731416807, -0.01611575463985900, 
                                      -0.0161157546398590,  0.01370293740923840,
                                       0.1012676956781630,  0.08151555513683590, 
                                       0.0815155551368359,  0.07687196074256860, 
                                       0.1790562778887310,  0.23984964946534200, 
                                       0.2398496494653420,  0.30677989063403600,
                                       0.5033867114419030,  0.71473583768765000,
                                       0.7147358376876500,  1.02665503604525000,
                                       0.3808678937603070,  0.55022173504226600,
                                       0.5502217350422660,  0.83625566582327500,
                                       0.3416908809010060,  0.48688054784433200,
                                       0.4868805478443320,  0.73255448906337700,
                                       0.3285801403106890,  0.45369928419248600,
                                       0.4536992841924860,  0.66891511235547800,
                                       0.2804228428696770,  0.39529498440455300,
                                       0.3952949844045530,  0.59106418309249200,
                                       0.6555683276776510,  0.78089315165869200,
                                       0.7808931516586920,  0.97963713553377100,
                                       1.1743994982232900,  1.09218887468317000,
                                       1.0921888746831700,  1.17246740859850000,
                                       0.9201474736423480,  0.71345132366320000,
                                       0.7134513236632000,  0.70210108917842500,
                                       0.5214482068419940,  0.35658375351207900,
                                       0.3565837535120790,  0.42789508317684000,
                                       0.4408478548952370,  0.16721284592264600,
                                       0.1672128459226460,  0.13645033701335600,
                                       0.3180525170107110,  0.04056079098580750,
                                       0.0405607909858075,  0.04826595734921470,
                                       0.3484832690734660,  0.02458514149442420,
                                       0.0245851414944242,  0.05122158251914640,
                                       0.4135751093189110,  0.04702349028691040,
                                       0.0470234902869104,  0.06751698458219490,
                                       0.3038248227581000,  0.06844779724119610,
                                       0.0684477972411961,  0.09073231649025490,
                                       0.1326189635924890,  0.07702722946203870,
                                       0.0770272294620387,  0.13354594669119600,
                                       0.3744307637848250,  0.30328891866570800,
                                       0.3032889186657080,  0.2456640241192240),
                                    .Dim = c(2L, 2L, 20L)))
  expect_error(CPAT:::get_lrv_arr(eX, kernel = "sqrt"))
  expect_error(CPAT:::get_lrv_arr(eX, kernel = sqrt))
  expect_error(CPAT:::get_lrv_arr(eX, kernel = ker_func, bandwidth = 1))

  check_cointReg()

  expect_error(CPAT:::get_lrv_arr(eX, kernel = "bar"))
  expect_error(CPAT:::get_lrv_arr(eX, kernel = "ba", bandwidth = "are"))
  expect_equal(CPAT:::get_lrv_arr(eX[1:3, ], kernel = "ba", bandwidth = "and"),
    structure(c(0.0377436292675294, -0.0275528493652965, -0.0275528493652965, 
                0.0201135800366664,  0.1738058309806920,  0.1108999179634800,
                0.1108999179634800,  0.1181484400166350,  0.4717011356410840,
                0.2358505678205420,  0.2358505678205420,  0.117925283910271),
              .Dim = c(2L, 2L, 3L)))
  expect_equal(CPAT:::get_lrv_arr(eX[1:3, ], kernel = "ba", bandwidth = "nw"),
    structure(c(0.0377436292675294, -0.0275528493652965, -0.0275528493652965,
                0.0201135800366664,  0.2628631842616570,  0.0944243076065014,
                0.0944243076065014,  0.0794084913394151,  0.4717011356410840,
                0.2358505678205420,  0.2358505678205420,  0.1179252839102710),
              .Dim = c(2L, 2L, 3L)))
  expect_equal(CPAT:::get_lrv_arr(eX[1:3, ], kernel = "pa", bandwidth = "and"),
    structure(c(0.0377436292675294, -0.0275528493652965, -0.0275528493652965,
                0.0201135800366664,  0.2286052227710890,  0.1007620304822570,
                0.1007620304822570,  0.0943107045878122,  0.4717011356410840,
                0.2358505678205420,  0.2358505678205420,  0.117925283910271),
              .Dim = c(2L, 2L, 3L)))
  expect_equal(CPAT:::get_lrv_arr(eX[1:3, ], kernel = "pa", bandwidth = "nw"),
    structure(c(0.0377436292675294, -0.0275528493652965, -0.0275528493652965,
                0.0201135800366664,  0.2170750690692540,  0.1028951089170960,
                0.1028951089170960,  0.0993263214481104,  0.4717011356410840,
                0.2358505678205420,  0.2358505678205420,  0.1179252839102710),
              .Dim = c(2L, 2L, 3L)))
  expect_equal(CPAT:::get_lrv_arr(eX[1:3, ], kernel = "tr"),
    structure(c(0.0377436292675294, -0.0275528493652965, -0.0275528493652965,
                0.0201135800366664, -0.0564095123713094,  0.1534897564836000,
                0.1534897564836000,  0.2182921143747550,  0.4717011356410840,
                0.2358505678205420,  0.2358505678205420,  0.1179252839102710),
              .Dim = c(2L, 2L, 3L)))
  expect_equal(CPAT:::get_lrv_arr(eX[1:3, ], kernel = "th"),
    structure(c(0.0377436292675294, -0.0275528493652965, -0.0275528493652965,
                0.0201135800366664,  0.1633100765625930,  0.1128416325308280,
                0.1128416325308280,  0.1227140931885080,  0.4717011356410840,
                0.2358505678205420,  0.2358505678205420,  0.1179252839102710),
              .Dim = c(2L, 2L, 3L)))
  expect_equal(CPAT:::get_lrv_arr(eX[1:3, ], kernel = "qs"),
    structure(c(0.0377436292675294, -0.0275528493652965, -0.0275528493652965,
                0.0201135800366664,  0.2046668652970270,  0.1051906266149580,
                0.1051906266149580,  0.1047238900890290,  0.4717011356410840,
                0.2358505678205420,  0.2358505678205420,  0.1179252839102710),
              .Dim = c(2L, 2L, 3L)))
  expect_equal(CPAT:::get_lrv_arr(eX[1:3, ], bandwidth = 1),
    structure(c(0.0377436292675294, -0.0275528493652965, -0.0275528493652965,
                0.0201135800366664,  0.2224799850369420,  0.1018951994630740,
                0.1018951994630740,  0.0969751830021661,  0.4717011356410840,
                0.2358505678205420,  0.2358505678205420,  0.1179252839102710),
              .Dim = c(2L, 2L, 3L)))
})

################################################################################
# INTERFACE TESTS
################################################################################

test_that("CUSUM.test() functions properly", {
  expect_s3_class(CUSUM.test(dat), "htest")
  expect_equal(CUSUM.test(dat)$data.name, "dat")
  expect_equal(CUSUM.test(dat)$p.value, 0.409099915842597)
  expect_equal(CUSUM.test(dat)$estimate, c(`t*` = 10))
  expect_equal(CUSUM.test(dat)$statistic, c(A = 0.888290332296761))
  expect_equal(CUSUM.test(dat, use_kernel_var = TRUE)$p.value,
               0.409673316747916)
  expect_error(CUSUM.test(c("a", "b")), "Don't know how to handle x of type" %s%
                                        "character")
  expect_error(CUSUM.test(df), "Formula needed")
  expect_equal(CUSUM.test(df, y ~ x)$statistic, c(A = 0.677882334483995))
})

test_that("HR.test() functions properly", {
  expect_s3_class(HR.test(dat), "htest")
  expect_equal(HR.test(dat)$data.name, "dat")
  expect_equal(HR.test(dat)$p.value, 0.00194727161751529)
  expect_equal(HR.test(dat)$estimate, c(`t*` = 2))
  expect_equal(HR.test(dat)$statistic, c(D = 3.48777644674448))
  expect_equal(HR.test(dat)$parameters, c(`log(T)` = 2.99573227355399))
  expect_equal(HR.test(dat, use_kernel_var = TRUE)$p.value,
               0.00349647304622447)
  expect_equal(HR.test(df2, y ~ x)$p.value, 0.0457699073347345)
  expect_equal(HR.test(df2, y ~ x, use_kernel_var = FALSE, kn = sqrt)$p.value,
               4.81098723914553e-06)
  expect_equal(HR.test(df2, y ~ x, use_kernel_var = FALSE,
                       kn = sqrt)$parameters,
               c(`sqrt(T)` = 4.47213595499958))
})

test_that("DE.test() functions properly", {
  expect_s3_class(DE.test(dat), "htest")
  expect_equal(DE.test(dat)$data.name, "dat")
  expect_equal(DE.test(dat)$p.value, 0.335048134286715)
  expect_equal(DE.test(dat)$estimate, c(`t*` = 10))
  expect_equal(DE.test(dat)$statistic, c(A = 1.58953572310406))
  expect_equal(DE.test(dat)$parameter, c(`a(log(T))` = 1.4813431070248, 
                                         `b(log(T))` = 1.66838804851426))
  expect_equal(DE.test(dat, use_kernel_var = TRUE)$p.value,
               0.335268000418675)
  expect_error(DE.test(c("a", "b")), "Don't know how to handle x of type" %s%
                                     "character")
  expect_error(DE.test(df), "Formula needed")
  expect_equal(DE.test(df, y ~ x)$statistic, c(A = 2.07256284012386))
})

test_that("HS.test() functions properly", {
  expect_s3_class(HS.test(dat), "htest")
  expect_equal(HS.test(dat)$data.name, "dat")
  expect_equal(HS.test(dat, m = 6)$statistic, c(A = 1.24993636457263))
  expect_equal(HS.test(dat)$p.value, 0.663825560827292)
  expect_equal(HS.test(dat)$statistic, c(A = 1.21370945972696))
  expect_equal(HS.test(dat, corr = FALSE)$p.value, 0.81061923831157)
  expect_error(HS.test(c("a", "b")), "Don't know how to handle x of type" %s%
                                      "character")
  expect_error(HS.test(df), "Formula needed")
  expect_equal(HS.test(df, y ~ x)$statistic, c(A = 4.31295336184296))
  expect_equal(HS.test(df, y ~ x, m = 6)$statistic, c(A = 1.90718155854372))
})

test_that("Andrews.test() functions properly", {
  expect_error(Andrews.test(dat), "argument \"M\" is missing")
  expect_s3_class(Andrews.test(dat, M = 15), "htest")
  expect_equal(Andrews.test(dat, M = 15)$data.name, "dat")
  expect_equal(Andrews.test(dat, M = 15)$p.value, 1)
  expect_equal(Andrews.test(dat, M = 15)$statistic, c(S = 2.41061543748607))
  expect_equal(Andrews.test(dat, M = 15)$parameters, c(m = 5))
  expect_equal(Andrews.test(df, 15, y ~ x)$statistic, c(S = 0.578318267262971))
})
