#' Imputation models for the missing S(1)
#'
#' @param formula Formula specifying the imputation model for the surrogate
#'   under treatment
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

    missdex <- !is.na(get(paste(formula[[2]]), psdesign$augdata))

    fit <- glm(formula, data = psdesign$augdata[missdex, ], family = distribution, weights = cdfweights, ...)

    psdesign$impute.model <- list(model = "parametric", args = arglist)

    mindelta <- subset(psdesign$augdata, !missdex)

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


#' Imputation models for the missing S(1)
#'
#' This model assumes that the pair [S(1), W] are bivariate normal, where W is
#' the BIP. The means, standard deviations, and correlation are estimated or
#' fixed before calling this function. Then the conditional normal formula is
#' applied in order to get the distribution of S(1) | S(0). That distribution is
#' used to impute the missing S(1) values. This method requires a BIP in the
#' design.
#'
#' @param mu, means of the pair of surrogates, missing one first
#' @param sd, standard deviations of the pair, missing one first
#' @param rho, the correlation between X1 and X2
#'
#' @export

impute_bivnorm <- function(mu = c(0, 0), sd = c(1, 1), rho = .2){

  rval <- function(psdesign){

    psdesign$impute.model <- "bivnorm"

    mindelta <- subset(psdesign$augdata, Z != 1)

    vmu <- mu[1] + sd[1]/sd[2] * rho * (mindelta$W - mu[2])
    vsd <- (1 - rho^2) * sd[1]^2


    psdesign$cdf_sbarw <-
      function(S.1){

        sapply(S.1, function(s) pnorm(s, mean = vmu, sd = vsd))

      }
    psdesign$icdf_sbarw <-
      function(U.1){

        sapply(U.1, function(s) qnorm(s, mean = vmu, sd = vsd))

      }
    ## return function of S and W that returns a probability, i.e., a cdf
    psdesign
  }

  class(rval) <- c("ps", "imputation")
  rval

}

