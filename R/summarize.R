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
#' @param CI.type Character string, "pointwise" for pointwise confidence intervals, and "band" for simultaneous confidence band.
#' @param n.samps The number of samples to take over the range of S.1 at which the VE is calculated
#' @param bootstraps If true, will calculate bootstrap standard errors and
#'   confidence bands.
#'
#' @export
#'
VE <- function(psdesign, t, sig.level = .05, CI.type = "band", n.samps = 5000, bootstraps = TRUE){

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

    ttt <- summary(survival::survfit(psdesign$augdata$Y ~ 1), rmean = TRUE)$table[["*rmean"]]

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

  if(bootstraps && "bootstraps" %in% names(psdesign)){

    bsests <- psdesign$bootstraps
    bootVEs <- matrix(NA, ncol = length(Splot) + 1, nrow = nrow(bsests))
    bootDiff <- matrix(NA, ncol = length(Splot) + 1, nrow = nrow(bsests))
    bootR1 <- matrix(NA, ncol = length(Splot) + 1, nrow = nrow(bsests))
    bootR0 <- matrix(NA, ncol = length(Splot) + 1, nrow = nrow(bsests))
    for(i in 1:nrow(bsests)){

      thispar <- as.numeric(bsests[i, -ncol(bsests)])
      if(inherits(psdesign$augdata$Y, "Surv") && missing(t)){

        ttt <- summary(survival::survfit(psdesign$augdata$Y ~ 1), rmean = TRUE)$table[["*rmean"]]
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
      bootDiff[i, ] <- c(R1 - R0, bsests[i, "convergence"])
      bootR1[i, ] <- c(R1, bsests[i, "convergence"])
      bootR0[i, ] <- c(R0, bsests[i, "convergence"])

    }

    bootVEs <- as.data.frame(bootVEs)
    bootR1 <- as.data.frame(bootR1)
    bootR0 <- as.data.frame(bootR0)
    bootDiff <- as.data.frame(bootDiff)

    colnames(bootVEs)[ncol(bootVEs)] <- colnames(bootDiff)[ncol(bootDiff)] <- colnames(bootR0)[ncol(bootR0)] <- colnames(bootR1)[ncol(bootR1)] <- "convergence"
    A1 <- as.data.frame(summarize_bs(bootR1, obsVE$R1, sig.level = sig.level, CI.type = CI.type)$table)
    A2 <- as.data.frame(summarize_bs(bootR0, obsVE$R0, sig.level = sig.level, CI.type = CI.type)$table)
    A3 <- as.data.frame(summarize_bs(bootVEs, obsVE$VE, sig.level = sig.level, CI.type = CI.type)$table)
    A4 <- as.data.frame(summarize_bs(bootDiff, obsVE$R1 - obsVE$R0, sig.level = sig.level, CI.type = CI.type)$table)

    colnames(A1) <- gsub("%", "", paste("R1", colnames(A1), sep = "."), fixed = TRUE)
    colnames(A2) <- gsub("%", "", paste("R0", colnames(A2), sep = "."), fixed = TRUE)
    colnames(A3) <- gsub("%", "", paste("VE", colnames(A3), sep = "."), fixed = TRUE)
    colnames(A4) <- gsub("%", "", paste("Rdiff", colnames(A4), sep = "."), fixed = TRUE)

    obsVE <- cbind(obsVE, A3, A2, A1, A4)

  }

  obsVE

}



#' Summarize bootstrap samples
#'
#' @param bootdf Data frame containing bootstrapped estimates, with a column containing a convergence indicator
#' @param estdf Data frame containing full sample estimate
#' @param sig.level Significance level to use for confidence intervals
#' @param CI.type Character string, "pointwise" for pointwise confidence intervals, and "band" for simultaneous confidence band.

summarize_bs <- function(bootdf, estdf = NULL, sig.level = .05, CI.type = "band") {

  bs <- bootdf[bootdf$convergence == 0, -which(colnames(bootdf) == "convergence")]

  if(CI.type == "pointwise"){
    mary <- function(x){
      funlist <- list(median = stats::median, mean = base::mean, boot.se = stats::sd,
                    lower.CL = function(x) quantile(x, sig.level/2),
                    upper.CL = function(x) quantile(x, 1 - sig.level/2))

      sapply(funlist, function(f) f(x))
    }

    table <- t(sapply(bs, mary))

  } else if(CI.type == "band"){

    estmat <- matrix(rep(estdf, nrow(bs)), nrow = nrow(bs), byrow = TRUE)
    maxdiff <- apply(abs(as.matrix(bs) - estmat), MARGIN = 1, FUN = max)
    inCI <- as.matrix(bs)[which(maxdiff < quantile(maxdiff, 1 - sig.level)), ]
    upper.CL <- apply(inCI, MARGIN = 2, FUN = max)
    lower.CL <- apply(inCI, MARGIN = 2, FUN = min)

    table <- data.frame(upper.CL = upper.CL, lower.CL = lower.CL)
    colnames(table) <- paste(colnames(table), 1-sig.level, sep = ".")

  }

  conv <- c(nboot = nrow(bootdf), ncov = sum(bootdf$convergence == 0))

  list(table = table, conv = conv)

}


#' Compute the empirical Vaccine Efficacy
#'
#' @param psdesign An object of class \link{psdesign}
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


