% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/design.R
\name{verify_trt}
\alias{verify_trt}
\title{Check that a variable is suitable for using as binary treatment indicator}
\usage{
verify_trt(D)
}
\arguments{
\item{D}{Vector that will be checked for 2-class labels}
}
\description{
Checks for two classes and gives a warning message indicating which level is assumed to be 0/1
}
\keyword{Internal}
