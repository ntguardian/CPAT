################################################################################
# Data.R
################################################################################
# 2018-09-22
# Curtis Miller
################################################################################
# Documentation for package datasets
################################################################################

#' Fama-French Five Factors
#'
#' Data set containing the five factors described by
#' \insertCite{famafrench15;textual}{CPAT}, from the data library maintained by
#' Kenneth French. Data ranges from July 1, 1963 to October
#' 31, 2017.
#'
#' @format A data frame with 13679 rows and 6 variables:
#' \describe{
#'   \item{Mkt.RF}{Market excess returns}
#'   \item{RF}{The risk-free rate of return}
#'   \item{SMB}{The return on a diversified portfolio of small stocks minus
#'              return on a diversified portfolio of big stocks}
#'   \item{HML}{The return of a portfolio of stocks with a high
#'              book-to-market (B/M) ratio minus the return of a portfolio
#'              of stocks with a low B/M ratio}
#'   \item{RMW}{The return of a portfolio of stocks with robust profitability
#'              minus a portfolio of stocks with weak profitability}
#'   \item{CMA}{The return of a portfolio of stocks with conservative
#'              investment minus the return of a portfolio of stocks with
#'              aggressive investment}
#' }
#'
#' Row names are dates in YYYYMMDD format.
#'
#' @source \url{http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}
"ff"

#' Bank Portfolio Returns
#'
#' Data set representing the returns of an industry portfolio representing the
#' banking industry based on company four-digit SIC codes, obtained from the
#' data library maintained by Kenneth French. Data ranges from July 1, 1926
#' to October 31, 2017.
#'
#' @format A data frame with 24099 rows and 1 variable:
#' \describe{
#'   \item{Banks}{The return of a portfolio representing the banking industry}
#' }
#'
#' Row names are dates in YYYY-MM-DD format.
#'
#' @source \url{http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}
"banks"
