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


  Splot <- seq(min(Sin), max(Sin), length.out = 1000) ## redo sampling from imputation model

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

  if(summary == "VE"){
    plot(VE ~ S.1, data = VE(psdesign, t), type = 'l', ...)
  } else {
    stop(paste("Plots of type", summary, "are not supported"))
  }

}


#' Concisely print information about a psdesign object
#'
#' @param psdesign A \link{psdesign} object
#' @param digits Number of significant digits to display
#' @param sig.level Significance level to use for computing bootstrapped confidence intervals
#'
#' @export
#'
print.psdesign <- function(psdesign, digits = 3, sig.level = .05){

  objs <- names(psdesign)

  cat("Augmented data frame: ", nrow(psdesign$augdata), " obs. by ", ncol(psdesign$augdata), " variables. \n")
  print(head(psdesign$augdata), digits = digits)

  cat("\nMapped variables: \n")
  temp <- lapply(names(psdesign$mapping), function(ln){

    lnm <- psdesign$mapping[[ln]]

    if(length(lnm) == 1){
      cat("\t", ln, " -> ", lnm, "\n")
      } else if(lnm[1] == "Surv") {

       cat("\t", ln, " -> ", paste0(lnm[1], "(", lnm[2], ", ", lnm[3], ")"), "\n")
      } else cat("\t", ln, " -> ", paste0(lnm, collapse = ", "), "\n")

  })

  cat("\nImputation models: \n")
  if(!"imputation.models" %in% objs) {

    cat("\tNone present, see ?impute for information on adding imputation models.\n")

  } else {

    for(j in names(psdesign$imputation.models)){
      cat("\t Imputation model for ", j, ":\n")
      tyj <- psdesign$imputation.models[[j]]$model
      cat("\t\t", paste0("impute_", tyj$model, "("))
      cat(paste(sapply(names(tyj$args[-1]), function(nj) paste0(nj, " = ", tyj$args[-1][nj])), collapse = ", "), ")\n")
    }

  }

  cat("\nRisk models: \n")
  if(!"risk.model" %in% objs) {
    cat("\tNone present, see ?risk for information on adding a risk model.\n")
  } else {

    tyj <- psdesign$risk.model
    cat("\t", paste0("risk_", tyj$model, "("))
    cat(paste(sapply(names(tyj$args[-1]), function(nj) paste0(nj, " = ", tyj$args[-1][nj])), collapse = ", "), ")\n\n")

  }

  if(!"estimates" %in% objs){

    cat("\tNo estimates present, see ?ps_estimate.\n")
  } else {
    cat("Estimated parameters:\n")
    print(psdesign$estimates$par, digits = digits)
    cat("\tConvergence: ", psdesign$estimates$convergence == 0, "\n\n")

  }

  if(!"bootstraps" %in% objs){

    cat("\tNo bootstraps present, see ?ps_bootstrap.\n")
  } else {

    cat("Bootstrap replicates:\n")
    sbs <- summarize_bs(psdesign$bootstraps, sig.level = sig.level)
    print(sbs$table, digits = digits)
    cat("\n\t Out of", sbs$conv[1], "bootstraps, ", sbs$conv[2], "converged (", round(100 * sbs$conv[2]/sbs$conv[1], 1), "%)\n")

  }

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
