#' Calculate the vaccine efficacy
#'
#' Computes the vaccince efficacy (VE) and risk in each treatment arm over the
#' range of surrogate values observed in the data. VE(s) is defined as 1 -
#' risk(s, z = 1)/risk(s, z = 0), where z is the treatment indicator. If any
#' other variables are present in the risk model, then the risk is computed at
#' their median value.
#'
#' @return A data frame containing columns for the S values, the VE, R0, and R1
#'   at those S values, and optionally standard errors and confidence intervals
#'   computed using bootstrapped estimates.
#'
#' @param psdesign A psdesign object. It must contain a risk model, an
#'   integration model, and estimated parameters. Bootstrapped parameters are
#'   optional
#' @param t For time to event outcomes, a fixed time \code{t} may be provided to
#'   compute the cumulative distribution function. If not, the restricted mean
#'   survival time is used. Omit for binary outcomes.
#' @param sig.level Significance level for bootstrap confidence intervals
#' @param bootstraps If true, will calculate bootstrap standard errors and
#'   confidence bands.
#'
#' @export
#'
VE <- function(psdesign, t, sig.level = .05, n.samps = 5000, bootstraps = TRUE){

  stopifnot("estimates" %in% names(psdesign))

  impped <- psdesign$integration.models$S.1$icdf_sbarw(runif(n.samps))
  randrows <- sample(1:nrow(impped), ncol(impped), replace = TRUE)
  integrated <- impped[cbind(randrows, 1:n.samps)]
  obss <- psdesign$augdata$S.1

  trueobs <- sample(obss[!is.na(obss)],
                    floor(n.samps * mean(!is.na(obss))),
                    prob = psdesign$augdata$cdfweights[!is.na(obss)], replace = TRUE)


  if(is.factor(obss)){
    integrated <- factor(integrated, levels = levels(obss))
  }
  Splot <- sort(unlist(list(integrated, trueobs)))

  dat1 <- data.frame(S.1 = Splot, Z = 1)
  dat0 <- data.frame(S.1 = Splot, Z = 0)


  if(is.null(psdesign$risk.model$args$model)){

    others <- NULL

  } else {

    mod <- eval(psdesign$risk.model$args$model)
    mterms0 <- rownames(attr(terms(mod), "factors"))[-1]
    mterms <- unlist(lapply(colnames(psdesign$augdata), function(x){
      if(length(grep(x, mterms0, fixed = TRUE) > 0)){
        return(x)
      } else return(NULL)
    }))

    others <- mterms[!mterms %in% c("S.1", "Z")]
    for(j in others){

      tempval <- median(psdesign$augdata[, j], na.rm = TRUE)
      dat1[, j] <- tempval
      dat0[, j] <- tempval

    }

  }

  if(inherits(psdesign$augdata$Y, "Surv") && missing(t)){

    ttt <- summary(survfit(psdesign$augdata$Y ~ 1), rmean = TRUE)$table[["*rmean"]]

    warning(sprintf("No time given for time to event outcome, using restricted mean survival: %.1f", ttt))
    R1 <- psdesign$risk.function(dat1, psdesign$estimates$par, t = ttt)
    R0 <- psdesign$risk.function(dat0, psdesign$estimates$par, t = ttt)


  } else if(inherits(psdesign$augdata$Y, "Surv")) {

    R1 <- psdesign$risk.function(dat1, psdesign$estimates$par, t)
    R0 <- psdesign$risk.function(dat0, psdesign$estimates$par, t)

  } else {

    R1 <- psdesign$risk.function(dat1, psdesign$estimates$par)
    R0 <- psdesign$risk.function(dat0, psdesign$estimates$par)

  }

  obsVE <- data.frame(S.1 = Splot, VE = 1 - R1/R0, R1 = R1, R0 = R0)

  if("bootstraps" %in% names(psdesign) & bootstraps){

    bsests <- psdesign$bootstraps
    bootVEs <- matrix(NA, ncol = length(Splot) + 1, nrow = nrow(bsests))
    bootR1 <- matrix(NA, ncol = length(Splot) + 1, nrow = nrow(bsests))
    bootR0 <- matrix(NA, ncol = length(Splot) + 1, nrow = nrow(bsests))
    for(i in 1:nrow(bsests)){

      thispar <- as.numeric(bsests[i, -ncol(bsests)])
      if(inherits(psdesign$augdata$Y, "Surv") && missing(t)){

        ttt <- summary(survfit(psdesign$augdata$Y ~ 1), rmean = TRUE)$table[["*rmean"]]
        R1 <- psdesign$risk.function(dat1, thispar, t = ttt)
        R0 <- psdesign$risk.function(dat0, thispar, t = ttt)


      } else if(inherits(psdesign$augdata$Y, "Surv")) {

        R1 <- psdesign$risk.function(dat1, thispar, t)
        R0 <- psdesign$risk.function(dat0, thispar, t)

      } else {

        R1 <- psdesign$risk.function(dat1, thispar)
        R0 <- psdesign$risk.function(dat0, thispar)

      }

      bootVEs[i, ] <- c(1 - R1/R0, bsests[i, "convergence"])
      bootR1[i, ] <- c(R1, bsests[i, "convergence"])
      bootR0[i, ] <- c(R0, bsests[i, "convergence"])

    }

    bootVEs <- as.data.frame(bootVEs)
    bootR1 <- as.data.frame(bootR1)
    bootR0 <- as.data.frame(bootR0)

    colnames(bootVEs)[ncol(bootVEs)] <- colnames(bootR0)[ncol(bootR0)] <- colnames(bootR1)[ncol(bootR1)] <- "convergence"
    A1 <- as.data.frame(summarize_bs(bootR1, sig.level = sig.level)$table)
    A2 <- as.data.frame(summarize_bs(bootR0, sig.level = sig.level)$table)
    A3 <- as.data.frame(summarize_bs(bootVEs, sig.level = sig.level)$table)

    colnames(A1) <- gsub("%", "", paste("R1", colnames(A1), sep = "."), fixed = TRUE)
    colnames(A2) <- gsub("%", "", paste("R0", colnames(A2), sep = "."), fixed = TRUE)
    colnames(A3) <- gsub("%", "", paste("VE", colnames(A3), sep = "."), fixed = TRUE)

    obsVE <- cbind(obsVE, A3, A2, A1)

  }

  obsVE

}



#' Summarize bootstrap samples
#'
#' @param bootdf Data frame containing bootstrapped estimates, with a column containing a convergence indicator
#' @param sig.level Significance level to use for confidence intervals

summarize_bs <- function(bootdf, sig.level = .05) {

  bs <- subset(bootdf, convergence == 0)[, -which(colnames(bootdf) == "convergence")]

  mary <- function(x){
    funlist <- list(median = stats::median, mean = base::mean, boot.se = stats::sd,
                  lower.CL = function(x) quantile(x, sig.level/2),
                  upper.CL = function(x) quantile(x, 1 - sig.level/2))

    sapply(funlist, function(f) f(x))
  }

  table <- t(sapply(bs, mary))
  conv <- c(nboot = nrow(bootdf), ncov = sum(bootdf$convergence == 0))

  list(table = table, conv = conv)

}


#' Compute the empirical Vaccine Efficacy
#'
#' @param psdesign
#' @param t Fixed time for time to event outcomes to compute VE. If missing, uses restricted mean survival.
#'
#' @export
empirical_VE <- function(psdesign, t){

  pd <- psdesign$augdata
  if(inherits(pd$Y, "Surv")){

    if(missing(t)){

      ttt <- summary(survival::survfit(pd$Y ~ 1), rmean = TRUE)$table[["*rmean"]]

    } else ttt <- t

    sf <- survival::survfit(pd$Y ~ pd$Z)
    sfs <- 1 - summary(sf, times = ttt, extend = TRUE)$surv

    1 - sfs[2]/sfs[1]

  } else {

    1 - mean(pd$Y[pd$Z == 1])/mean(pd$Y[pd$Z == 0])

  }

}


