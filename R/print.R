
#' Plot summary statistics for a psdesign object
#'
#' Plot the vaccine efficacy versus S.1 for an estimated psdesign object
#'
#' @param psdesign A psdesign object that contains a risk model, integration
#'   model, and valid estimates
#' @param t For time to event outcomes, a fixed time \code{t} may be provided to
#'   compute the cumulative distribution function. If not, the restricted mean
#'   survival time is used. Omit for binary outcomes.
#' @param summary Summary statistic to plot. Currently only \link{VE} is
#'   supported.
#' @param sig.level Significance level used for confidence bands on the VE
#'   curve. This is only used if bootstrapped estimates are available.
#'   @param n.samps Number of samples to use over the range of S.1 for plotting the curve
#' @param ... Other arguments passes to \link{plot}
#'
#' @export

plot.psdesign <- function(psdesign, t, summary = "VE", sig.level = .05, n.samps = 500, ...){

  if(summary == "VE"){

    VE.me <- VE(psdesign, t, sig.level = sig.level, n.samps = n.samps)
    plot(VE ~ S.1, data = VE.me, type = 'l', ...)

    if("VE.boot.se" %in% colnames(VE.me)){

      lnme <- grep("VE.lower.", colnames(VE.me), fixed = TRUE)
      unme <- grep("VE.upper.", colnames(VE.me), fixed = TRUE)

      if(is.factor(VE.me[, 1])){

        subVE <- unique(VE.me)
        segments(rep.int(as.integer(subVE[, "S.1"]), 2) - .4, rep.int(subVE[, lnme], 2),
                 x1 = rep.int(as.integer(subVE[, "S.1"]), 2) + .4, lty = 3, ...)

        segments(rep.int(as.integer(subVE[, "S.1"]), 2) - .4, rep.int(subVE[, unme], 2),
                 x1 = rep.int(as.integer(subVE[, "S.1"]), 2) + .4, lty = 3, ...)

      } else {
        lines(VE.me[, lnme] ~ VE.me$S.1, lty = 3, ...)
        lines(VE.me[, unme] ~ VE.me$S.1, lty = 3, ...)
        }
    }

  } else if(summary == "RR"){

    VE.me <- VE(psdesign, t, sig.level = sig.level, n.samps = n.samps)
    plot(1 - VE ~ S.1, data = VE.me, type = 'l', ...)

    if("VE.boot.se" %in% colnames(VE.me)){

      lnme <- grep("VE.lower.", colnames(VE.me), fixed = TRUE)
      unme <- grep("VE.upper.", colnames(VE.me), fixed = TRUE)

      if(is.factor(VE.me[, 1])){

        subVE <- unique(VE.me)
        segments(rep.int(as.integer(subVE[, "S.1"]), 2) - .4, rep.int(1 - subVE[, lnme], 2),
                 x1 = rep.int(as.integer(subVE[, "S.1"]), 2) + .4, lty = 3, ...)

        segments(rep.int(as.integer(subVE[, "S.1"]), 2) - .4, rep.int(1 - subVE[, unme], 2),
                 x1 = rep.int(as.integer(subVE[, "S.1"]), 2) + .4, lty = 3, ...)

      } else {
        lines(1 - VE.me[, lnme] ~ VE.me$S.1, lty = 3, ...)
        lines(1 - VE.me[, unme] ~ VE.me$S.1, lty = 3, ...)
        }
    }

  } else if(summary == "logRR"){

    VE.me <- VE(psdesign, t, sig.level = sig.level, n.samps = n.samps)
    plot(log(1 - VE) ~ S.1, data = VE.me, type = 'l', ...)

    if("VE.boot.se" %in% colnames(VE.me)){

      lnme <- grep("VE.lower.", colnames(VE.me), fixed = TRUE)
      unme <- grep("VE.upper.", colnames(VE.me), fixed = TRUE)

      if(is.factor(VE.me[, 1])){

        subVE <- unique(VE.me)
        segments(rep.int(as.integer(subVE[, "S.1"]), 2) - .4, rep.int(1 - subVE[, lnme], 2),
                 x1 = rep.int(as.integer(subVE[, "S.1"]), 2) + .4, lty = 3, ...)

        segments(rep.int(as.integer(subVE[, "S.1"]), 2) - .4, rep.int(1 - subVE[, unme], 2),
                 x1 = rep.int(as.integer(subVE[, "S.1"]), 2) + .4, lty = 3, ...)

      } else {
        lines(log(1 - VE.me[, lnme]) ~ VE.me$S.1, lty = 3, ...)
        lines(log(1 - VE.me[, unme]) ~ VE.me$S.1, lty = 3, ...)
        }
    }

  } else if(summary == "risk") {

      VE.me <- VE(psdesign, t, sig.level = sig.level, n.samps = n.samps)



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
  pout <- NULL

  cat("Augmented data frame: ", nrow(psdesign$augdata), " obs. by ", ncol(psdesign$augdata), " variables. \n")
  print(head(psdesign$augdata), digits = digits)

  cat("\nEmpirical VE: ", round(empirical_VE(psdesign), digits), "\n")

  cat("\nMapped variables: \n")
  temp <- lapply(names(psdesign$mapping), function(ln){

    lnm <- psdesign$mapping[[ln]]

    if(length(lnm) == 1){
      cat("\t", ln, " -> ", lnm, "\n")
    } else if(lnm[1] == "Surv") {

      cat("\t", ln, " -> ", paste0(lnm[1], "(", lnm[2], ", ", lnm[3], ")"), "\n")
    } else cat("\t", ln, " -> ", paste0(lnm, collapse = ", "), "\n")

  })

  cat("\nIntegration models: \n")
  if(!"integration.models" %in% objs) {

    cat("\tNone present, see ?add_integration for information on integration models.\n")

  } else {

    for(j in names(psdesign$integration.models)){
      cat("\t integration model for ", j, ":\n")
      tyj <- psdesign$integration.models[[j]]$model
      cat("\t\t", paste0("integrate_", tyj$model, "("))
      cat(paste(sapply(names(tyj$args[-1]), function(nj) paste0(nj, " = ", tyj$args[-1][nj])), collapse = ", "), ")\n")
    }

  }

  cat("\nRisk models: \n")
  if(!"risk.model" %in% objs) {
    cat("\tNone present, see ?add_riskmodel for information on risk models.\n")
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

    pout$boot.table <- sbs

  }

  invisible(pout)


}

#' Summary method for psdesign objects
#'
#' @param psdesign A \link{psdesign} object
#' @param digits Number of significant digits to display
#' @param sig.level Significance level to use for computing bootstrapped confidence intervals
#'
#' @export
#'
summary.psdesign <- function(psdesign, digits = 3, sig.level = .05){

  pout <- print(psdesign, digits = digits, sig.level = sig.level)

  ## compute marginal model and summarize VE

  if("risk.model" %in% names(psdesign)){
    pdat <- psdesign$risk.model$args
    pdat$model <- Y ~ Z

    psdesign2 <- psdesign + do.call(as.character(pdat[[1]]), pdat[-1])
    marg.est <- psdesign2 + ps_estimate()
    marg.VE <- colMeans(VE(marg.est, bootstraps = FALSE))[2]
    emp.VE <- empirical_VE(psdesign)
    cond.VE <- VE(psdesign, bootstraps = FALSE)
    cond.VE.est <- 1 - mean(cond.VE$R1)/mean(cond.VE$R0)
    VEtab <- c(empirical = emp.VE, marginal = marg.VE, model = cond.VE.est)
    print(VEtab, digits = digits)

    empdiff <- 100 * VEtab[3]/VEtab[1] - 100
    mardiff <- 100 * VEtab[3]/VEtab[2] - 100

    cat(sprintf("Model-based average VE is %.1f %% different from the empirical and %.1f %% different from the marginal.\n", empdiff, mardiff))
    if(abs(mardiff) > 100) warning("Check model and results carefully!")

    invisible(list(print = pout, VE.estimates = VEtab))

  }
}
