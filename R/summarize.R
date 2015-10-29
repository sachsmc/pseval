#' Calculate the vaccine efficacy
#'
#' Computes the vaccince efficacy (VE) over the range of surrogate values
#' observed in the data. VE(s) is defined as 1 - risk(s, z = 1)/risk(s, z = 0),
#' where z is the treatment indicator. If any other variables are present in the
#' risk model, then the risk is computed at their median value.
#'
#' @return A data frame containing columns for the S values, the VE at those S
#'   values, and optionally standard errors computed using bootstrapped
#'   estimates.
#'
#' @param psdesign A psdesign object. It must contain a risk model, an
#'   imputation model, and estimated parameters. Bootstrapped parameters are
#'   optional
#' @param t For time to event outcomes, a fixed time \code{t} may be provided to
#'   compute the cumulative distribution function. If not, the restricted mean
#'   survival time is used. Omit for binary outcomes.
#'
#' @export
#'
VE <- function(psdesign, t){

  stopifnot("estimates" %in% names(psdesign) | "bootstraps" %in% names(psdesign))

  Sin <- psdesign$augdata$S.1
  Sin <- Sin[!is.na(Sin)]


  Splot <- seq(min(Sin), max(Sin), length.out = 1000)
  dat1 <- data.frame(S.1 = Splot, Z = rep(1, 1000))
  dat0 <- data.frame(S.1 = Splot, Z = rep(0, 1000))


  if(is.null(psdesign$risk.model$args$model)){

    others <- NULL

  } else {

    mod <- eval(psdesign$risk.model$args$model)
    mterms <- rownames(attr(terms(mod), "factors"))[-1]

    others <- mterms[!mterms %in% c("S.1", "Z")]
    for(j in others){

      tempval <- median(psdesign$augdata[, j], na.rm = TRUE)
      dat1[, j] <- tempval
      dat0[, j] <- tempval

    }

  }

  if(inherits(psdesign$augdata$Y, "Surv") && missing(t)){

    ttt <- summary(survfit(psdesign$augdata$Y ~ 1), rmean = TRUE)$table[["*rmean"]]
    R1 <- psdesign$risk.function(dat1, psdesign$estimates$par, t = ttt)
    R0 <- psdesign$risk.function(dat0, psdesign$estimates$par, t = ttt)


  } else if(inherits(psdesign$augdata$Y, "Surv")) {

    R1 <- psdesign$risk.function(dat1, psdesign$estimates$par, t)
    R0 <- psdesign$risk.function(dat0, psdesign$estimates$par, t)

  } else {

    R1 <- psdesign$risk.function(dat1, psdesign$estimates$par)
    R0 <- psdesign$risk.function(dat0, psdesign$estimates$par)

  }

  data.frame(S.1 = Splot, VE = 1 - R1/R0)

}


#' Plot summary statistics for a psdesign object
#'
#' Plot the vaccine efficacy versus S.1 for an estimated psdesign object
#'
#' @param psdesign A psdesign object that contains a risk model, imputation model, and valid estimates
#' @param t For time to event outcomes, a fixed time \code{t} may be provided to
#'   compute the cumulative distribution function. If not, the restricted mean
#'   survival time is used. Omit for binary outcomes.
#'   @param summary Summary statistic to plot. Currently only \link{VE} is supported.
#'  @param ... Other arguments passes to \link{plot}
#'
#'  @export

plot.psdesign <- function(psdesign, t, summary = "VE", ...){

  plot(VE ~ S.1, data = VE(psdesign, t), type = 'l', ...)

}





