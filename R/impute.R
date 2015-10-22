#' Parametric imputation model for the missing S(1)
#'
#' @param formula Formula specifying the imputation model for the surrogate
#'   under treatment. Generally the candidate surrogate will be on the left side
#'   in the formula, and the BIP or BIPs will be on the right side
#' @param distribution Assumed distribution for the imputation model. Must be
#'   compatible with the \code{family} argument of \link{glm}. Currenly only
#'   Gaussian models are supported
#' @param ... Arguments passed to \link{glm}
#'
#' @export

impute_parametric <- function(formula, distribution = gaussian, ...){

  stopifnot(identical(distribution, gaussian))

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    if(!"imputation.models" %in% names(psdesign)) psdesign$imputation.models <- NULL
    outname <- paste(formula[[2]])

    missdex <- !is.na(get(paste(formula[[2]]), psdesign$augdata))

    fit <- glm(formula, data = psdesign$augdata[missdex, ], family = distribution, weights = cdfweights, ...)

    psdesign$imputation.models[[outname]]$model <- list(model = "parametric", args = arglist)

    mindelta <- subset(psdesign$augdata, !missdex)

    psdesign$imputation.models[[outname]]$cdf_sbarw <-
      function(S.1){

        mu <- predict(fit, newdata = mindelta, type = "response")
        sd <- sd(fit$residuals)

        sapply(S.1, function(s) pnorm(s, mean = mu, sd = sd))

      }
    psdesign$imputation.models[[outname]]$icdf_sbarw <-
      function(U.1){

        mu <- predict(fit, newdata = mindelta, type = "response")
        sd <- sd(fit$residuals)

        sapply(U.1, function(u) qnorm(u, mean = mu, sd = sd))

      }

    psdesign
  }

  class(rval) <- c("ps", "imputation")
  rval

}


#' Imputation models for the missing S(1)
#'
#' This model assumes that the pair [S(1), W] are bivariate normal, where W is
#' the BIP. The means, standard deviations, and correlation are estimated or
#' fixed before calling this function. Then the conditional normal formula is
#' applied in order to get the distribution of S(1) | W. That distribution is
#' used to impute the missing S(1) values. This method requires a BIP in the
#' design.
#'
#' @param x, expression identifying the variable to be imputed. Typically this is S.1 or S.0
#' @param mu, means of the pair of surrogates, missing one first
#' @param sd, standard deviations of the pair, missing one first
#' @param rho, the correlation between X1 and X2
#'
#' @export

impute_bivnorm <- function(x = S.1, mu = c(0, 0), sd = c(1, 1), rho = .2){

  outname <- as.character(substitute(x))

  rval <- function(psdesign){


    if(!"imputation.models" %in% names(psdesign)) psdesign$imputation.models <- NULL

    psdesign$imputation.models[[outname]]$impute.model <- "bivnorm"

    missdex <- !is.na(get(outname, psdesign$augdata))
    mindelta <- subset(psdesign$augdata, !missdex)

    vmu <- mu[1] + sd[1]/sd[2] * rho * (mindelta$BIP - mu[2])
    vsd <- (1 - rho^2) * sd[1]^2


    psdesign$imputation.models[[outname]]$cdf_sbarw <-
      function(S.1){

        sapply(S.1, function(s) pnorm(s, mean = vmu, sd = vsd))

      }
    psdesign$imputation.models[[outname]]$icdf_sbarw <-
      function(U.1){

        sapply(U.1, function(s) qnorm(s, mean = vmu, sd = vsd))

      }

    psdesign
  }

  class(rval) <- c("ps", "imputation")
  rval

}


#' Semiparametric imputation model using the location-scale model
#'
#' @param formula.location Formula specifying the imputation model for the location component of the surrogate
#'   under treatment. Generally the candidate surrogate will be on the left side
#'   in the formula, and the BIP or BIPs will be on the right side
#
#' @param formula.scale Formula specifying the imputation model for the scale component of the surrogate
#'   under treatment. Generally the candidate surrogate will be on the left side
#'   in the formula, and the BIP or BIPs will be on the right side

impute_semiparametric <- function(formula.location, formula.scale, ...){
  arglist <- as.list(match.call())

  rval <- function(psdesign){

    if(!"imputation.models" %in% names(psdesign)) psdesign$imputation.models <- NULL

    missdex <- !is.na(get(paste(formula.location[[2]]), psdesign$augdata))

    fit <- sp_locscale(formula.location = formula.location, formula.scale = formula.scale,
                       data = psdesign$augdata[missdex, ], weights = cdfweights, ...)

    psdesign$imputation.models[[paste(formula[[2]])]]$model <- list(model = "semiparametric", args = arglist)

    mindelta <- subset(psdesign$augdata, !missdex)

    psdesign$imputation.models[[paste(formula[[2]])]]$cdf_sbarw <-
      function(S.1){

        mu <- predict(fit, newdata = mindelta, type = "response")
        sd <- sd(fit$residuals)

        sapply(S.1, function(s) pnorm(s, mean = mu, sd = sd))

      }
    psdesign$imputation.models[[paste(formula[[2]])]]$icdf_sbarw <-
      function(U.1){

        mu <- predict(fit, newdata = mindelta, type = "response")
        sd <- sd(fit$residuals)

        sapply(U.1, function(u) qnorm(u, mean = mu, sd = sd))

      }

    psdesign
  }

  class(rval) <- c("ps", "imputation")
  rval

}


sp_locscale <- function(formula.location, formula.scale, data, weights, ...){
  # do stuff
}

