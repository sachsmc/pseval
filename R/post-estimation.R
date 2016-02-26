#' Calculate the Standardized total gain
#'
#' Computes the standardized total gain for the risk difference. Optionally produces bootstrap standard errors, and permutation test statistic. The standardized total gain is the area between the risk difference curve and the horizontal line at the marginal risk difference.
#'
#' @param psdesign A psdesign object. It must contain a risk model, an
#'   integration model, and estimated parameters. Bootstrapped parameters are
#'   optional
#' @param t For time to event outcomes, a fixed time \code{t} may be provided to
#'   compute the cumulative distribution function. If not, the restricted mean
#'   survival time is used. Omit for binary outcomes.
#' @param sig.level Significance level for bootstrap confidence intervals
#' @param n.samps The number of samples to take over the range of S.1 at which
#'   the VE is calculated
#' @param bootstraps If true, and bootstrapped estimates are present, will
#'   calculate bootstrap standard errors and confidence bands.
#'
#' @export
#'
calc_STG <- function(psdesign, t, sig.level = .05, n.samps = 5000, bootstraps = TRUE){

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

  risks <- riskcalc(psdesign$risk.function, psdesign$augdata$Y, psdesign$estimates$par, t, dat0, dat1)

  obsdelta <- data.frame(S.1 = Splot, delta = risks$R1 - risks$R0)
  obstheta <- mean(obsdelta$delta)

  obsSTG <- stg(risks$R1, risks$R0, Splot)

  if(bootstraps && "bootstraps" %in% names(psdesign)){

    bsests <- psdesign$bootstraps
    bootSTGs <- matrix(NA, nrow = nrow(bsests), ncol = 2)

    for(i in 1:nrow(bsests)){

      thispar <- as.numeric(bsests[i, -ncol(bsests)])

      thisrisk <- riskcalc(psdesign$risk.function, psdesign$augdata$Y,
                           thispar, t, dat0, dat1)

      bootSTGs[i, ] <- c(stg(thisrisk$R1, thisrisk$R0, Splot), bsests[i, "convergence"])

    }

    bootSTGs <- as.data.frame(bootSTGs)
    colnames(bootSTGs)[ncol(bootSTGs)] <- "convergence"
    A1 <- as.data.frame(summarize_bs(bootSTGs, obsSTG, sig.level = sig.level, CI.type = "pointwise")$table)

    colnames(A1) <- gsub("%", "", paste("STG", colnames(A1), sep = "."), fixed = TRUE)

    obsSTG <- cbind(obsSTG, A1)

  }

  obsSTG

}


#' Compute the standardized total gain
#' @param R1 Risk in the treatment group
#' @param R0 Risk in the control group
#' @param S Values of S at which risk is calculated

stg <- function(R1, R0, S){

  delt <- R0 - R1
  theta <- mean(delt)
  mean(abs(delt - theta)) / (2 * theta * (1 - theta))

}