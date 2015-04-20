#' Specify a design for a principal surrogate evaluation
#'

psdesign <- function(data = NULL, treatment, outcome,
                     BIP = NULL, CPV = NULL, BSM = NULL, weights = NULL, ...){


  trt <- eval(substitute(treatment), envir = data)
  oot <- eval(substitute(outcome), envir = data)

  if(is.null(weights)) weights <- 1
  rval <- NULL

  rval$augdata <- data.frame(Z = trt, cdfweights = weights, cdfind = trt == 1, data)
  class(rval) <- "psdesign"

  rval

}