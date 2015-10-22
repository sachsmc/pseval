#' Add imputation model to a psdesign object
#'
#' @param psdesign A psdesign object
#' @param imputation An imputation object, such as \link{impute_parametric}
#'
#' @export

add_imputation <- function(psdesign, imputation){

  stopifnot(inherits(psdesign, "psdesign"))
  stopifnot(inherits(imputation, "imputation"))

  psdesign <- imputation(psdesign)
  psdesign

}

#' Add risk model to a psdesign object
#'
#' @param psdesign A psdesign object
#' @param risk A risk model object, such as \link{risk_logistic}
#'
#' @export

add_riskmodel <- function(psdesign, riskmodel){

  stopifnot(inherits(psdesign, "psdesign"))
  stopifnot(inherits(riskmodel, "riskmodel"))

  psdesign <- riskmodel(psdesign)
  psdesign

}


#' Modify a psdesign object by adding on new components.
#'
#' This operator allows you to add objects to a psdesign object, such as imputation models and risk models
#'
#' If the first object is an object of class \code{psdesign}, you can add
#' the following types of objects, and it will return a modified psdesign
#' object.
#'
#' \itemize{
#'   \item \code{imputation}: Add or replace imputation model
#'   \item \code{riskmodel}: Add or replace risk model
#' }
#'
#' @export

"+.ps" <- function(p1, p2){

  if(inherits(p2, "imputation")){
    add_imputation(p1, p2)
    } else if(inherits(p2, "riskmodel")){
      add_riskmodel(p1, p2)
    }

}
