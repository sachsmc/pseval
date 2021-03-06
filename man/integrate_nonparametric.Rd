% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integrate.R
\name{integrate_nonparametric}
\alias{integrate_nonparametric}
\title{Nonparametric integration model for the missing S(1)}
\usage{
integrate_nonparametric(formula, ...)
}
\arguments{
\item{formula}{Formula specifying the integration model for the surrogate
under treatment. Generally the candidate surrogate will be on the left side
in the formula, and the BIP or BIPs will be on the right side. In this case
the BIP and the S(1) must be categorical.}

\item{...}{Not currently used}
}
\description{
Both S(1) and the BIP or set of BIPs must be categorical. This model integrates over the estimated distribution of S(1) | BIP
}
