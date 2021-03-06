% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChangePointTests.R
\name{HR.test}
\alias{HR.test}
\title{Rényi-Type Test}
\usage{
HR.test(x, kn = log, use_kernel_var = FALSE, stat_plot = FALSE,
  kernel = "ba", bandwidth = "and")
}
\arguments{
\item{x}{Data to test for change in mean}

\item{kn}{A function corresponding to the trimming parameter \eqn{t_T}; by
default, the square root function}

\item{use_kernel_var}{Set to \code{TRUE} to use kernel methods for long-run
variance estimation (typically used when the data is
believed to be correlated); if \code{FALSE}, then the
long-run variance is estimated using
\eqn{\hat{\sigma}^2_{T,t} = T^{-1}\left(
\sum_{s = 1}^t \left(X_s - \bar{X}_t\right)^2 +
\sum_{s = t + 1}^{T}\left(X_s -
\tilde{X}_{T - t}\right)^2\right)}, where
\eqn{\bar{X}_t = t^{-1}\sum_{s = 1}^t X_s} and
\eqn{\tilde{X}_{T - t} = (T - t)^{-1}
\sum_{s = t + 1}^{T} X_s}; if \code{custom_var} is not
\code{NULL}, this argument is ignored}

\item{stat_plot}{Whether to create a plot of the values of the statistic at
all potential change points}

\item{kernel}{If character, the identifier of the kernel function as used in
\pkg{cointReg} (see \code{\link[cointReg]{getLongRunVar}}); if
function, the kernel function to be used for long-run variance
estimation (default is the Bartlett kernel in \pkg{cointReg})}

\item{bandwidth}{If character, the identifier for how to compute the
bandwidth as defined in \pkg{cointReg} (see
\code{\link[cointReg]{getBandwidth}}); if function, a function
to use for computing the bandwidth; if numeric, the bandwidth
value to use (the default is to use Andrews' method, as used in
\pkg{cointReg})}
}
\value{
A \code{htest}-class object containing the results of the test
}
\description{
Performs the (univariate) Rényi-type test for change in mean, as described in
\insertCite{horvathricemiller19}{CPAT}. This is effectively an interface to
\code{\link{stat_Zn}}; see its documentation for more details. p-values are
computed using \code{\link{pZn}}, which represents the limiting distribution
of the test statistic under the null hypothesis, which represents the
limiting distribution of the test statistic under the null hypothesis when
\code{kn} represents a sequence \eqn{t_T} satisfying \eqn{t_T \to \infty}
and \eqn{t_T/T \to 0} as \eqn{T \to \infty}. (\code{\link[base]{log}} and
\code{\link[base]{sqrt}} should be good choices.)
}
\examples{
HR.test(rnorm(1000))
HR.test(rnorm(1000), use_kernel_var = TRUE, kernel = "bo", bandwidth = "nw")
}
\references{
\insertAllCited{}
}
