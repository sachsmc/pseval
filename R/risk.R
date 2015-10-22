#' Risk model for binary outcome
#'
#' @param model Formula specifying the risk model
#' @param D number of samples for the simulated annealing integration
#' @param risk Function for transforming a linear predictor into a probability.
#'   E.g., risk.expit for the logistic model, risk.probit for the probit model


risk_binary <- function(model = Y ~ S.1 * Z, D = 5000, risk = risk.expit, ...){

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    # refactor to include possible CPV and/or BSM
    # break of into separate function called expand_augdata

    trtmat <- model.matrix(model, psdesign$augdata[psdesign$augdata$Z == 1, ])
    Y.trt <- psdesign$augdata[psdesign$augdata$Z == 1, ]$Y

    untrtsamp <- c(psdesign$icdf_sbarw(runif(D)))
    dex <- (1:nrow(psdesign$augdata))[psdesign$augdata$Z == 0]
    untrtobs <- psdesign$augdata[rep(dex, D), ]
    untrtobs$S.1 <- untrtsamp

    untrt.expand <- model.matrix(model, untrtobs)
    Y.untrt <- untrtobs$Y

    likelihood <- function(beta){

      trted <- risk(trtmat %*% beta)^Y.trt *
        (1 - risk(trtmat %*% beta))^(1 - Y.trt)

      # for each W in untrted, sample a bunch from cdf_sbarw and take the mean

      untrted <- matrix(risk(untrt.expand %*% beta)^Y.untrt *
        (1 - risk(untrt.expand %*% beta))^(1 - Y.untrt), nrow = D, byrow = TRUE)

      sum(log(trted)) + sum(log(colMeans(untrted)))

    }

    psdesign$likelihood <- likelihood
    psdesign$risk.model <- list(model = "logistic", args = arglist)
    psdesign$nparam <- ncol(trtmat)

    psdesign

  }
  ## return a likelihood closure

  class(rval) <- c("ps", "riskmodel")
  rval
}


#' Weibull risk model for time to event outcome
#'
#' @param model Formula specifying the risk model. The outcome should be a \link{Surv} object specifying right censoring
#' @param D number of samples for simulated annealing
#' @param ...
#'
risk_weibull <- function(model = Y ~ S.1 * Z, D = 5000, ... ){

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    trtmat <- model.matrix(model, psdesign$augdata[psdesign$augdata$Z == 1, ])
    Y.trt <- psdesign$augdata[psdesign$augdata$Z == 1, ]$Y[,1]
    delt.trt <- psdesign$augdata[psdesign$augdata$Z == 1, ]$Y[, 2]

    untrtsamp <- c(psdesign$icdf_sbarw(runif(D)))
    dex <- (1:nrow(psdesign$augdata))[psdesign$augdata$Z == 0]
    untrtobs <- psdesign$augdata[rep(dex, D), ]
    untrtobs[ , attr(terms.formula(model), "term.labels")[1]] <- untrtsamp

    untrt.expand <- model.matrix(model, untrtobs)
    Y.untrt <- untrtobs$Y[, 1]
    delt.untrt <- untrtobs[, 2]

    likelihood <- function(beta){


      gamma0<-beta[1]
      beta0<-beta[-1]
      shapepram<-exp(gamma0)
      scalepram<-exp(trtmat %*% beta0)

      trtlike <- ((shapepram/scalepram)*(Y.trt/scalepram)^(shapepram-1))^delt.trt*exp(-(Y.trt/scalepram)^shapepram)

      # for each W in untrted, sample a bunch from cdf_sbarw and take the mean

      scale.untrt <- exp(untrt.expand %*% beta0)

      untrted <- matrix(((shapepram/scale.untrt)*(Y.untrt/scale.untrt)^(shapepram-1))^delt.untrt *
                          exp(-(Y.untrt/scale.untrt)^shapepram), nrow = D, byrow = TRUE)

      sum(log(trtlike)) + sum(log(colMeans(untrted)))

    }

    psdesign$likelihood <- likelihood
    psdesign$risk.model <- list(model = "weibull", args = arglist )
    psdesign$nparam <- ncol(trtmat) + 1
    psdesign$censored <- cens

    psdesign

  }
  ## return a likelihood closure

  class(rval) <- c("ps", "riskmodel")
  rval

}


risk.expit <- function(x) {

  exp(x)/(1 + exp(x))

}


risk.probit <- function(x) {

  pnorm(x)

}

