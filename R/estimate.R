#' Estimate parameters from a specified model using estimated maximum likelihood
#'
#' @param psdesign A psdesign object with a risk model and imputation model specified
#' @param start Vector of starting values, if NULL, will come up with starting values
#' @param control List of control parameters for passed to \link{optim}
#' @param ... Arguments passed to \link{optim}
#'
ps_estimate <- function(psdesign, start = NULL, control = list(fnscale = -1), ...){

  stopifnot("impute.model" %in% names(psdesign) && "risk.model" %in% names(psdesign))

  if(is.null(start)){

    start <- rep(0, psdesign$nparam)

  }

  est1 <- optim(start, fn = psdesign$likelihood, control = control, ...)

  est1

}


#' Estimate parameters from a specified model using bootstrap resampling and estimated maximum likelihood
#'
#' @param psdesign A psdesign object with a risk model and imputation model specified
#' @param n.boots Number of bootstrap replicates
#' @param progress.bar Logical, if true will display a progress bar in the console
#' @param start Vector of starting values, if NULL, will come up with starting values
#' @param control List of control parameters for passed to \link{optim}
#' @param ... Arguments passed to \link{optim}
#


ps_bootstrap <- function(psdesign, n.boots = 200, progress.bar = TRUE, start = NULL, control = list(fnscale = -1), ...){

  stopifnot("impute.model" %in% names(psdesign) && "risk.model" %in% names(psdesign))

  bootpar <- vector(mode = "list", length = n.boots)

  if(progress.bar){
     pb <- txtProgressBar(min = 1, max = n.boots, title = paste("Bootstrapping", n.boots, "replicates"))
  }

  psdesign.0 <- psdesign
  for(i in 1:n.boots){
    # resample augdata

    sampdex <- sample(1:nrow(psdesign$augdata), nrow(psdesign$augdata), replace = TRUE)
    psdesign.0$augdata <- psdesign$augdata[sampdex, ]

    ## re-call imputation model

    psdesign2 <- psdesign.0 + do.call(as.character(psdesign$impute.model$args[[1]]), psdesign$impute.model$args[-1])

    ## re-call risk model

    psdesign3 <- psdesign2 + do.call(as.character(psdesign$risk.model$args[[1]]), psdesign$risk.model$args[-1])

    # estimate

    bpar <- ps_estimate(psdesign3, start = start, control = control, ...)

    bootpar[[i]] <- c(bpar$par, convergence = bpar$convergence)

    if(progress.bar){
      setTxtProgressBar(pb, value = i)
      flush.console()
    }

  }

  if(progress.bar) close(pb)

  as.data.frame(do.call(rbind, bootpar))

}
