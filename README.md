---
bibliography: 'inst/REFERENCES.bib'
output:
  md_document:
    variant: 'markdown\_github'
---

<!-- README.md is generated from README.Rmd. Please edit that file -->
CPAT
====

**CPAT** is a package implementing some statistical tests for detecting
structural change in a series of data. The tests made publicly available
are:

-   The CUSUM test, via `CUSUM.test()` [see @auehorvath13]
-   The Darling-Erdös test, via `DE.test()`
-   The Hidalgo-Seo test, via `HS.test()`, as described in
    [@hidalgoseo13]
-   Andrews test, via `Andrews.test()`, as described in [@andrews03]
-   The Rényi-type test as described by Horváth, Rice and Miller in
    [-@horvathricemiller19], via `HR.test()`

This package was written to facilitate the simulations performed in
[@horvathricemiller19] and thus is geared to change point tests capable
of detecting early/late changes in a sample. That said, it is general
purpose.

Change Point Testing
--------------------

Change point testing is performed on sequential data (such as time
series) to determine whether the data shares a common structure. In
particular, let $X_t = \mu_t + \epsilon_t$ where $\epsilon_t$ is a
random noise process, with $1 \leq t \leq T$. With the exception of
Andrews' test, the tests mentioned above can be used to decide between
the hypotheses:

$$H_0 : \mu_1 = \ldots \mu_T = \mu$$

$$H_A : \mu_1 = \ldots \mu_{t^{*} - 1} \neq \mu_{t^{*}} = \ldots = \mu_T$$

$t^{*}$ is an integer satisfying $1 \leq t^{*} \leq T$. Critically,
$t^{*}$ is not assumed to be known, so the alternative hypothesis states
that the change occurs at an unknown location in the sample.

Andrews' test is an exception; his test assumes some information about
where the change occured. Suppose $\mu_t = \mu$ for
$1 \leq t \leq T^{'}$ with $T^{'} < T$ and $T^{'}$ known. His test
decides between the hypotheses:

$$H_0: \mu_{T^{'} + 1} = \ldots = \mu_{T}$$

$$H_A: \mu_{T^{'} + 1} = \ldots = \mu_{t^{*} - 1} \neq \mu_{t^{*}} = \ldots =
\mu_{T}$$

In this case, $T^{'} \leq t^{*} \leq T$.

Change point testing traces its roots to quality control procedures; if
one imagines $X_t$ being some measurement of a part produced by an
machine, $H_0$ states that the machine was always calibrated while the
alternative hypothesis claimes that the machine became uncalibrated at
some unknown porint in time. Another view of a change point test is that
it's yet another test to check if a time series is stationary, which is
a critical assumption made in most analyses involving time series.

While the above formulation describes a test primarily concerned about
the mean of a process, these tests, even in this form, can make
statements about structural change other than the mean. For instance, if
$X_t$ represents the residuals of a regression model, the test checks
whether the parameters of the model are stable over time or not.

All of the tests included in this package are asymptotic tests; the test
performs better for large $T$ and one should be cautious when using
these tests for small $T$. Simulation studies comparing these tests
(except for Andrews' test) are presented in [@horvathricemiller19]. In
short, the Rényi-type test seems to perform best when a change occurs
near the ends of a sample. In unreported simulations studies, the CUSUM
test seemed to perform best when the change occured mid-sample.

Example
-------

Each function can accept a univariate dataset as input, and will return
an `htest`-class object with the results of the test.

``` {.r}
library(CPAT)
#>       ________ _________ _________ __________
#>      /       //  ___   //  ___   //         /
#>     /   ____//  /  /  //  /  /  //___   ___/
#>    /   /    /  /__/  //  /__/  /    /  /
#>   /   /___ /  ______//  ___   /    /  /
#>  /       //  /      /  /  /  /    /  /
#> /_______//__/      /__/  /__/    /__/        v. 0.0.0.9000
#> Type citation("CPAT") for citing this R package in publications

set.seed(20180924)
x <- c(rnorm(250), rnorm(50) + 1)

CUSUM.test(x)
#> 
#>  CUSUM Test for Change in Mean
#> 
#> data:  x
#> A = 1.8875, p-value = 0.001609
#> sample estimates:
#>  t* 
#> 250
DE.test(x)
#> 
#>  Darling-Erdos Test for Change in Mean
#> 
#> data:  x
#> A = 3.7528, a(log(T)) = 1.8661, b(log(T)) = 3.1872, p-value =
#> 0.04582
#> sample estimates:
#>  t* 
#> 251
HS.test(x)
#> 
#>  Hidalgo-Seo Test for Change in Mean
#> 
#> data:  x
#> A = 11.783, Correlated Residuals = 1, p-value = 0.00551
#> sample estimates:
#>  t* 
#> 251
HR.test(x)
#> 
#>  Horvath-Rice Test for Change in Mean
#> 
#> data:  x
#> D = 2.2145, log(T) = 5.7038, p-value = 0.1043
#> sample estimates:
#>  t* 
#> 284
Andrews.test(x, M = 250)
#> 
#>  Andrews' Test for Structural Change
#> 
#> data:  x
#> S = 52.63, m = 50, p-value = 0.2388
```

Often these tests return a test statistic, a $p$-value, and the
estimated location of the change. This latter quantity is the arg-max of
the terms the test statistics maximize. (There is theory supporting this
estimator for the CUSUM test, but no theory for the Rényi-type test.)

`Andrews.test()` also allows for testing for structural change in a
linear model directly.

``` {.r}
df <- data.frame(x = x, y = 1 + 2 * x + c(rnorm(250), rnorm(50) + 1))
Andrews.test(df, M = 250, formula = y ~ x)
#> 
#>  Andrews' Test for Structural Change
#> 
#> data:  df
#> S = 6.3494, m = 50, p-value < 2.2e-16
```

References {#references .unnumbered}
----------
