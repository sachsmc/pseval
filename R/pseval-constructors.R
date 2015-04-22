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
#' }
#'
#' @export

"+.ps" <- function(p1, p2){

  if(inherits(p2, "imputation")) add_imputation(p1, p2)

}
