
#' Plot summary statistics for a psdesign object
#'
#' Plot the vaccine efficacy versus S.1 for an estimated psdesign object
#'
#' @param psdesign A psdesign object that contains a risk model, integration
#'   model, and valid estimates
#' @param t For time to event outcomes, a fixed time \code{t} may be provided to
#'   compute the cumulative distribution function. If not, the restricted mean
#'   survival time is used. Omit for binary outcomes.
#' @param summary Summary statistic to plot. \code{"VE"} for vaccine efficacy =
#'   1 - risk_1(s)/risk_0(s), \code{"RR"} for relative risk =
#'   risk_1(s)/risk_0(s), \code{"logRR"} for log of the relative risk,
#'   \code{"risk"} for the risk in each treatment arm, and \code{"RD"} for the
#'   risk difference = risk_1(s) - risk_0(s).
#' @param sig.level Significance level used for confidence bands on the VE
#'   curve. This is only used if bootstrapped estimates are available.
#' @param n.samps Number of samples to use over the range of S.1 for plotting
#'   the curve
#'   @param col Vector of integers specifying colors for each curve.
#'   @param lty Vector of integers specifying linetypes for each curve.
#'   @param lwd Vector of numeric values for line widths.
#' @param ... Other arguments passed to \link{plot}
#'
#' @export

plot.psdesign <- function(psdesign, t, summary = "VE", sig.level = .05, n.samps = 500, xlab = "S.1", ylab = summary, col = 1, lty = 1, lwd = 1, ...){


  VE.me <- VE(psdesign, t, sig.level = sig.level, n.samps = n.samps)
  n.curve.base <- ifelse(summary == "risk", 2, 1)
  ncurve <- ifelse("VE.boot.se" %in% colnames(VE.me), n.curve.base * 3, n.curve.base)

  ## some logic taken from plot.survfit


    if (length(lty)==1 && is.numeric(lty))
      lty <- rep(c(lty, lty+1, lty+1), 2)
    else if (length(lty) < ncurve)
      lty <- rep(rep(lty, each=3), length.out=(ncurve*3))
    else lty <- rep(lty, length.out= ncurve*3)

    if (length(col) <= ncurve) col <- rep(rep(col, each=3), length.out=3*ncurve)
    else col <- rep(col, length.out=3*ncurve)

    if (length(lwd) <= ncurve) lwd <- rep(rep(lwd, each=3), length.out=3*ncurve)
    else lwd <- rep(lwd, length.out=3*ncurve)


  mainme <- switch(summary, VE = parse(text = "VE"),
                   RR = parse(text = "1 - VE"),
                   logRR = parse(text = "log(1 - VE)"),
                   risk = parse(text = "list(R1, R0)"),
                   riskdiff = parse(text = "R1 - R0"))

  if(is.null(mainme)) stop(paste("Plots of summary", summary, "are not supported."))

  if(is.factor(VE.me[, 1])){
    envir <- unique(VE.me)
  } else envir <- VE.me

  if(summary == "risk"){

    rlist <- eval(mainme, envir = envir)
    plot(rlist[[1]] ~ envir[, 1], type = 'l', col = col[1], lty = lty[1], lwd = lwd[1], ylab = ylab, xlab = xlab)
    if(is.factor(VE.me[, 1])){

      segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(rlist[[2]], 2),
               x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[4], lty = lty[4], lwd = lwd[4])

    } else {

      lines(rlist[[2]] ~ envir[, 1], type = 'l', col = col[4], lty = lty[4], lwd = lwd[4])

    }

  } else {

    plot(eval(mainme, envir = envir) ~ envir[, 1], col = col[1], lty = lty[1], lwd = lwd[1], type = 'l', ylab = ylab, xlab = xlab)

  }


  if("VE.boot.se" %in% colnames(VE.me)){

    lnme <- switch(summary, VE = parse(text = grep("VE.lower.", colnames(VE.me), fixed = TRUE, value = TRUE)),
                   RR =  parse(text = paste("1 - ", grep("VE.lower.", colnames(VE.me), fixed = TRUE, value = TRUE))),
                   logRR =  parse(text = paste("log(1 - ", grep("VE.lower.", colnames(VE.me), fixed = TRUE, value = TRUE), ")")),
                   risk = parse(text = paste("list(", paste(sapply(c("R1.lower.", "R0.lower."),
                                                                   function(x) grep(x, colnames(VE.me), fixed = TRUE, value = TRUE)), collapse = ", "), ")")),
                   riskdiff = parse(text = paste(sapply(c("R1.lower.", "R0.lower."), function(x) grep(x, colnames(VE.me), fixed = TRUE, value = TRUE)), collapse = " - ")))

    unme <- switch(summary, VE = parse(text = grep("VE.upper.", colnames(VE.me), fixed = TRUE, value = TRUE)),
                   RR =  parse(text = paste("1 - ", grep("VE.upper.", colnames(VE.me), fixed = TRUE, value = TRUE))),
                   logRR =  parse(text = paste("log(1 - ", grep("VE.upper.", colnames(VE.me), fixed = TRUE, value = TRUE), ")")),
                   risk = parse(text = paste("list(", paste(sapply(c("R1.upper.", "R0.upper."),
                                                                   function(x) grep(x, colnames(VE.me), fixed = TRUE, value = TRUE)), collapse = ", "), ")")),
                   riskdiff = parse(text = paste(sapply(c("R1.upper.", "R0.upper."), function(x) grep(x, colnames(VE.me), fixed = TRUE, value = TRUE)), collapse = " - ")))

    if(summary == "risk"){

      rlistl <- eval(lnme, envir = envir)
      rlistu <- eval(unme, envir = envir)

      if(is.factor(VE.me[, 1])){

        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(rlistu[[1]], 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[2], lty = lty[2], lwd = lwd[2], ...)
        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(rlistu[[2]], 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[5], lty = lty[5], lwd = lwd[5], ...)
        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(rlistl[[1]], 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[3], lty = lty[3], lwd = lwd[3], ...)
        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(rlistl[[2]], 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[6], lty = lty[6], lwd = lwd[6], ...)

      } else {

        lines(rlistu[[1]] ~ envir[, 1], type = 'l', col = col[2], lty = lty[2], lwd = lwd[2])
        lines(rlistu[[2]] ~ envir[, 1], type = 'l', col = col[5], lty = lty[5], lwd = lwd[5])
        lines(rlistl[[1]] ~ envir[, 1], type = 'l', col = col[3], lty = lty[3], lwd = lwd[3])
        lines(rlistl[[2]] ~ envir[, 1], type = 'l', col = col[6], lty = lty[6], lwd = lwd[6])

      }

    } else {

      if(is.factor(VE.me[, 1])){

        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(eval(unme, envir = envir), 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[2], lty = lty[2], lwd = lwd[2], ...)
        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(eval(lnme, envir = envir), 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[3], lty = lty[3], lwd = lwd[3], ...)

      } else {

        lines(eval(unme, envir = envir) ~ envir[, 1], col = col[2], lty = lty[2], lwd = lwd[2], type = 'l')
        lines(eval(lnme, envir = envir) ~ envir[, 1], col = col[3], lty = lty[3], lwd = lwd[3], type = 'l')
      }
    }


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
    marg.VE <- mean(VE(marg.est, bootstraps = FALSE)[, 2])
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
