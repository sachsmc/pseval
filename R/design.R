#' Specify a design for a principal surrogate evaluation
#'
#' @param data Data frame containing data to be analyzed
#' @param treatment Expression defining the treatment variable
#' @param outcome Expression defining the outcome variable
#' @param BIP Optional expression defining the baseline irrelevant predictor
#' @param CPV Optional expression defining the closeout placebo vaccination measurement
#' @param BSM Optional expression defining the baseline surrogate measurement
#' @param weights optional expression defining weights to accomodate nonrandom subsampling, such as case control or two phase
#' @param ... Not currently used
#'
#' @export

psdesign <- function(data = NULL, treatment, outcome,
                     BIP = NULL, CPV = NULL, BSM = NULL, weights = NULL, ...){


  trt <- eval(substitute(treatment), envir = data)
  oot <- eval(substitute(outcome), envir = data)

  if(is.null(weights)) weights <- 1
  rval <- NULL

  rval$augdata <- data.frame(Z = trt, cdfweights = weights, cdfind = trt == 1, data)
  class(rval) <- c("ps", "psdesign")

  rval

}