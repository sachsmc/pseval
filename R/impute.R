#' Imputation models for the missing S(1)
#'
#' @param formula Formula specifying the imputation model for the surrogate under treatment
#' @param distribution Assumed distribution for the imputation model. Must be compatible with the \code{family} argument of \link{glm}
#' @param ... Arguments passed to \link{glm}
#'
#' @export

impute_parametric <- function(formula, distribution = gaussian, ...){

  rval <- function(psdesign){
    bdata <- psdesign$augdata
    bdata <- bdata[bdata$cdfind, ]

    fit <- glm(formula, data = bdata, family = distribution, weights = cdfweights, ...)

    psdesign$impute_model <- "parametric"

    mindelta <- subset(psdesign$augdata, Z != 1)

    psdesign$cdf_sbarw <-
      function(S.1){

        mu <- predict(fit, newdata = mindelta, type = "response")
        sd <- sd(fit$residuals)

        sapply(S.1, function(s) pnorm(s, mean = mu, sd = sd))

      }
    psdesign$icdf_sbarw <-
      function(U.1){

        mu <- predict(fit, newdata = mindelta, type = "response")
        sd <- sd(fit$residuals)

        sapply(U.1, function(u) qnorm(u, mean = mu, sd = sd))

      }
    ## return function of S and W that returns a probability, i.e., a cdf
    psdesign
  }

  class(rval) <- c("ps", "imputation")
  rval

}
