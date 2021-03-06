% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProbabilityFunctions.R
\name{pdarling_erdos}
\alias{pdarling_erdos}
\title{Darling-Erdös Statistic CDF}
\usage{
pdarling_erdos(q)
}
\arguments{
\item{q}{Quantile input to CDF}
}
\value{
If \eqn{Z} is the random variable with this distribution, the
        quantity \eqn{P(Z \leq q)}
}
\description{
CDF for the limiting distribution of the Darling-Erdös statistic.
}
\examples{
CPAT:::pdarling_erdos(0.1)
}
