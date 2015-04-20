#' Imputation models for the missing S(1)
#'

impute_parametric <- function(psdesign, formula, distribution = gaussian, ...){

  stopifnot(inherits(psdesign, what = "psdesign"))

  bdata <- psdesign$augdata
  bdata <- bdata[bdata$cdfind, ]

  fit <- glm(formula, data = bdata, family = distribution, weights = cdfweights, ...)

  psdesign$impute_model <- "parametric"
  psdesign$cdf_sbarw <-
    function(S.1){

      mu <- predict(fit, newdata = psdesign$augdata, type = "response")
      sd <- sd(fit$residuals)

      sapply(S.1, function(s) pnorm(s, mean = mu, sd = sd))

    }
  ## return function of S and W that returns a probability, i.e., a cdf

  psdesign

}
