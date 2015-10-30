#' Imputation models
#'
#' Add imputation model to a psdesign object
#'
#' @details This is a list of the available imputation models. The fundamental problem in surrogate evaluation is that there are unobserved values of the counterfactual surrogate reponses S(1). In the estimated maximum likelihood framework, for subjects missing the S(1) values, we use an auxiliary pre-treatment variable or set of variables W that is observed for every subject to estimate the distribution of S(1) | W. Typically, this W is a BIP. Then for each missing S(1), we impute likelihood contributions for each non-missing S(1) given their value of W, and average over the contributions.
#'
#' \itemize{
#' \item \link{impute_parametric} This is a parametric imputation model that fits a linear model for the mean of S(1) | W and assumes a Gaussian distribution.
#' \item \link{impute_bivnorm} This is another parametric imputation model that assumes that S(1) and W are jointly normally distributed. The user must specify their mean, variances and correlation.
#' \item \link{impute_nonparametric} This is a non-parametric imputation model that is only valid for categorical S(1) and W. It uses the observed proportions to estimate the joint distribution of S(1), W.
#' \item \link{impute_semiparametric} This is a semi-parametric model that uses the semi-parametric location scale model of Heagerty and Pepe (1999). Models are specified for the location of S(1) | W and the scale of S(1) | W. Then imputations are drawn from the empirical distribution of the residuals from that model, which are then transformed to the appropriate location and scale.
#' }
#'
#'
#' @param psdesign A psdesign object
#' @param imputation An imputation object
#'
#' @export
#'
#' @examples
#'
#' test <- psdesign(generate_example_data(n = 100), Z = Z, Y = Y.obs, S = S.1, BIP = X)
#'
#' add_imputation(test, impute_parametric())
#' test + impute_parametric()  # same as above
#'

add_imputation <- function(psdesign, imputation){

  stopifnot(inherits(psdesign, "psdesign"))
  stopifnot(inherits(imputation, "imputation"))

  psdesign <- imputation(psdesign)
  psdesign

}

#' Add risk model to a psdesign object
#'
#' @details The risk model component specifies the likelihood for the data. This involves specifying the distribution of the outcome variable, whether it is binary or time-to-event, and specifying how the surrogate S(1) and the treatment Z interact and affect the outcome. We use the formula notation to be consistent with other regression type models in R. Below is a list of available risk models.
#'
#' \itemize{
#' \item \link{risk_binary} This is a generic risk model for binary outcomes. The user can specify the formula, and link function using either \link{risk.expit} for the logistic link, or \link{risk.probit} for the probit link. Custom link functions may also be specified, which take a single numeric vector argument, and returns a vector of corresponding probabilities.
#' \item \link{risk_weibull} This is a parameterization of the Weibull model for time-to-event outcomes that is consistent with that of \link{rweibull}. The user specifies the formula for the linear predictor of the scale parameter.
#' \item \link{risk_exponential} This is a simple exponential model for a time-to-event outcome.
#' }
#'
#' @param psdesign A psdesign object
#' @param risk A risk model object
#'
#' @export
#'
#' @examples
#' test <- psdesign(generate_example_data(n = 100), Z = Z, Y = Y.obs, S = S.1, BIP = X)
#' add_riskmodel(test, risk_binary())
#' test + risk_binary() # same as above

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
