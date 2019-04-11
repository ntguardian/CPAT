# CPAT News

## Version 0.2.0
---

- Updated `DE.test()` to work with linear regression models (computes residuals
    of model then applies test to them); added `formula` argument
- Updated `HS.test()` to work with linear regression models (implements Hidalgo
    and Seo's test when derived for linear regression models); added `formula`
    argument
- Added `m` argument to `HS.test()` that adds control for how many terms are
    summed when estimating the long-run variance
- Updated `CUSUM.test()` to work with linear regression models (computes
    residuals of model then applies test to them); added `formula` argument
- `Andrews.test()` now gives a different error message when data is not numeric
    or if a regression setup was incorrectly specified
- `DE.test()`, `HS.test()`, `CUSUM.test()`, and `Andrews.test()` now need
    argument `x` to be either `numeric` or a `data.frame`, and will complain if
    this is not the case
- Internal function `get_lrv_vec()` was rewritten; long-run variance estimation
    for affected functions (`DE.test()`, `HS.test()`, `CUSUM.test()`) should now
    be much faster, especially when `kernel = "qs"`
- Added data set `usa_gbr_price` containing CPI and exchange rates of the United
    States and Great Britain from January 1971 to February 2019.
- Added data set `prod_comp` containing U.S. worker hourly productivity and
    compensation indices from 1947 to 2016.
- Added data set `energy_demand` containing U.S. GDP, energy consumption, and
    energy pricing from 1973 to 2018.

## Version 0.1.0
---

- Created `Andrews.test()`, a function that implements D. W. K. Andrews' test
    for end-of-sample instability, both for univariate and regression data
- Created `DE.test()`, a function that implements the Darling-Erdös test for
    structural change
- Created `HS.test()`, a function that implements the Hidalgo-Seo test for
    structural change
- Created `HR.test()`, a function that implements the Horváth et. al. test for
    structural change
- Created `CUSUM.test()`, a function that implements the CUSUM test for
    structural change
- Created `ff`, a data set containing Fama-French five factor data
- Created `banks`, a data set containing the returns of a portfolio of banking
    sector stocks, as compiled by K. French
- Files and functions related to simulations done for the accompanying paper
    were removed

## Version 0.0.0.9000
---

### Setup

- Package **CPAT** created.
