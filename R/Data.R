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

#' United States and United Kingdom CPI and Exchange Rates
#'
#' Data set containing the CPI of the United Kingdom and the United States and
#' the exchange rate between the countries' currencies. Data is at monthly
#' frequency starting January 1971 and ending February 2019.
#'
#' CPI data is from the OECD database (downloaded April 5th, 2019) and exchange
#' rate data is from the FRED database (downloaded April 6th, 2019). 
#'
#' @format A data frame with 578 rows and 3 variables:
#' \describe{
#'   \item{USA_CPI}{United States CPI, with 2015 having an index of 100.}
#'   \item{GBR_CPI}{United Kingdom CPI, with 2015 having an index of 100.}
#'   \item{GBPUSD}{The exchange rate between the British pound sterling and the
#'                 U.S. dollar, in pounds per dollar.}
#' }
#'
#' Row names are month and year.
#'
#' @source \href{https://fred.stlouisfed.org/series/EXUSUK}{FRED},
#'         \href{https://data.oecd.org/price/inflation-cpi.htm}{OECD}
"usa_gbr_price"

#' Productivity and Compensation Indices for the U.S. Economy
#'
#' Data containing indices for (hourly) productivity and compensation of U.S.
#' workers from the U.S. Bureau of Labor Statistics. (2009 is 100.)
#'
#' @format A time series object with 279 rows and 2 variables:
#' \describe{
#'   \item{productivity}{Labor productivity (output per hour).}
#'   \item{compensation}{Real hourly compensation in business sector}
#' }
#'
#' @source \href{https://www.bls.gov/}{U.S. Bureau of Labor Statistics}
"prod_comp"

#' U.S. Quarterly Energy Consumption, Prices, and GDP
#'
#' Data set containing U.S. total consumption of energy and average price per
#' quarter, along with U.S. GDP. The period spans from Q1 1973 to Q4 2018.
#'
#' Raw data was obtained from the Federal Reserve Economic Database and from
#' Data.gov.
#'
#' @format A time series object with 184 rows and 3 variables:
#' \describe{
#'   \item{gdp}{U.S. GDP in billions of dollars, from FRED}
#'   \item{price}{Index of urban energy prices, from FRED, with 1982-1984 being
#'                100; made quarterly by averaging the monthly index values}
#'   \item{consumption}{Total U.S. primary energy consumption, from Data.gov, in
#'                      trillions of Btu}
#' }
#'
#' @source \href{https://fred.stlouisfed.org/series/GDP}{FRED: GDP}
#'         \href{https://fred.stlouisfed.org/series/CPIENGSL}{FRED: Energy CPI}
#'         \href{https://catalog.data.gov/dataset/monthly-energy-consumption-by-sector}{Data.gov}
"energy_demand"

#' U.S. Natural Gas Consumption, Prices, and GDP
#'
#' Data set containing U.S. total consumption of natural gas and average price
#' per quarter, along with U.S. GDP. The period spans from Q1 1976 to Q4 2012.
#'
#' Raw data was obtained from the Federal Reserve Economic Database and from the
#' U.S. Energy Information Administration.
#'
#' @format A time series object with 148 rows and 3 variables:
#' \describe{
#'   \item{gdp}{U.S. GDP in millions of dollars, from FRED}
#'   \item{price}{Average natural gas price in dollars per thousand cubic feet,
#'                from the EIA}
#'   \item{consumption}{Natural gas consumption in billions of cubic feet, from
#'                      the EIA}
#'   \item{population}{U.S. quarterly population in thousands, from FRED}
#'   \item{pcgdp}{U.S. quarterly per capita GDP}
#' }
#' 
#' @source \href{https://fred.stlouisfed.org/series/NA000334Q}{FRED: GDP}
#'         \href{https://fred.stlouisfed.org/series/B230RC0Q173SBEA}{
#'               FRED: Population}
#'         \href{https://www.eia.gov/totalenergy/data/browser/}{EIA}
"natural_gas_demand"

#' Daily Log Returns of Corrections Corporation of America/CoreCivic Stock (CXW)
#'
#' Data set containing daily log returns (as a percentage) for the stock of
#' Corrections Corporation of America/CoreCivic, a private prison company, from
#' July 16, 1997 to October 31, 2017. Stock prices are adjusted for stock splits
#' and dividends and was obtained from Quandl, via the command
#' \code{Quandl::Quandl("WIKI/CXW")} (see \code{\link[Quandl]{Quandl}}).
#'
#' @format A \code{link[zoo]{zoo}} time series object with 5109 rows and 1
#'         variable:
#' \describe{
#'   \item{CXW}{The daily (natural) log return of CXW, as a percentage}
#' }
#' @source \href{https://www.quandl.com/}
"CXW"
