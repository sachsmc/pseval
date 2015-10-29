#' Risk model for binary outcome
#'
#' @param model Formula specifying the risk model
#' @param D number of samples for the simulated annealing integration
#' @param risk Function for transforming a linear predictor into a probability.
#'   E.g., risk.expit for the logistic model, risk.probit for the probit model
#'   @export


risk_binary <- function(model = Y ~ S.1 * Z, D = 5000, risk = risk.expit, ...){

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    expanded <- expand_augdata(model, psdesign, D = D)

    trtmat <- expanded$noimp
    Y.trt <- expanded$noimp.Y

    untrt.expand <- expanded$imp
    Y.untrt <- expanded$imp.Y

    likelihood <- function(beta){

      trted <- risk(trtmat %*% beta)^Y.trt *
        (1 - risk(trtmat %*% beta))^(1 - Y.trt)

      # for each W in untrted, sample a bunch from cdf_sbarw and take the mean

      untrted <- matrix(risk(untrt.expand %*% beta)^Y.untrt *
        (1 - risk(untrt.expand %*% beta))^(1 - Y.untrt), nrow = D, byrow = TRUE)

      -1 * (sum(log(trted)) + sum(log(colMeans(untrted))))

    }

    psdesign$risk.function <- function(data, beta){  ## P(D = 1 | S, Z)

      risk(as.vector(model.matrix(model[-2], data) %*% beta))

    }

    psdesign$likelihood <- likelihood
    psdesign$risk.model <- list(model = "binary", args = arglist)
    psdesign$nparam <- ncol(trtmat)
    psdesign$param.names <- colnames(trtmat)

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
#'@export

risk_weibull <- function(model = Y ~ S.1 * Z, D = 5000, ... ){

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    expanded <- expand_augdata(model, psdesign, D = D)

    if(!inherits(expanded$noimp.Y, "Surv")) stop("Requires a survival outcome specified with Surv(time, event)")

    trtmat <- expanded$noimp
    Y.trt <- expanded$noimp.Y[, 1]
    delt.trt <- expanded$noimp.Y[, 2]

    untrt.expand <- expanded$imp
    Y.untrt <- expanded$imp.Y[, 1]
    delt.untrt <- expanded$imp.Y[, 2]

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

      -1 * (sum(log(trtlike)) + sum(log(colMeans(untrted))))

    }

    psdesign$risk.function <- function(data, beta, t){  #P(Y < t | S, Z)

      mat <- model.matrix(model[-2], data)
      shape <- exp(beta[1])
      beta0 <- beta[-1]

      scale <- as.vector(exp(mat %*% beta0))

      1 - exp(-(t / scale)^shape)

    }

    psdesign$likelihood <- likelihood
    psdesign$risk.model <- list(model = "weibull", args = arglist )
    psdesign$nparam <- ncol(trtmat) + 1
    psdesign$param.names <- c("shape", colnames(trtmat))

    psdesign

  }
  ## return a likelihood closure

  class(rval) <- c("ps", "riskmodel")
  rval

}



#' Exponential risk model for time to event outcome
#'
#' @param model Formula specifying the risk model. The outcome should be a \link{Surv} object specifying right censoring
#' @param D number of samples for simulated annealing
#' @param ...
#'
#'@export

risk_exponential <- function(model = Y ~ S.1 * Z, D = 5000, ... ){

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    expanded <- expand_augdata(model, psdesign, D = D)

    if(!inherits(expanded$noimp.Y, "Surv")) stop("Requires a survival outcome specified with Surv(time, event)")

    trtmat <- expanded$noimp
    Y.trt <- expanded$noimp.Y[, 1]
    delt.trt <- expanded$noimp.Y[, 2]

    untrt.expand <- expanded$imp
    Y.untrt <- expanded$imp.Y[, 1]
    delt.untrt <- expanded$imp.Y[, 2]

    likelihood <- function(beta){

      scalepram <- 1/exp(trtmat %*% beta)

      trtlike <- scalepram^delt.trt * exp(-scalepram * Y.trt)

      # for each W in untrted, sample a bunch from cdf_sbarw and take the mean

      scale.untrt <- 1/exp(untrt.expand %*% beta)

      untrted <- matrix(scale.untrt^delt.untrt * exp(-scale.untrt * Y.untrt), nrow = D, byrow = TRUE)

      -1 * (sum(log(trtlike)) + sum(log(colMeans(untrted))))

    }

    psdesign$risk.function <- function(data, beta, t){ # P(T < t | S, Z)

      scale <- as.vector(exp(model.matrix(model[-2], data) %*% beta))

      1 - exp(-scale * t)

    }

    psdesign$likelihood <- likelihood
    psdesign$risk.model <- list(model = "exponential", args = arglist )
    psdesign$nparam <- ncol(trtmat)
    psdesign$param.names <- colnames(trtmat)

    psdesign

  }
  ## return a likelihood closure

  class(rval) <- c("ps", "riskmodel")
  rval

}


#' @export

risk.expit <- function(x) {

  exp(x)/(1 + exp(x))

}

#' @export
risk.probit <- function(x) {

  pnorm(x)

}


expand_augdata <- function(model, psdesign, D = 500){


  vars <- paste(attr(terms(model), "variables"))[-1]
  stopifnot("S.1" %in% vars)
  noimpdex <- !is.na(psdesign$augdata["S.1"])

  trtmat <- model.matrix(model, psdesign$augdata[noimpdex, ])
  Y.trt <- psdesign$augdata[noimpdex, ]$Y

  untrtsamp <- c(psdesign$imputation.models$S.1$icdf_sbarw(runif(D)))
  dex <- (1:nrow(psdesign$augdata))[!noimpdex]
  untrtobs <- psdesign$augdata[rep(dex, D), ]
  untrtobs$S.1 <- untrtsamp

  untrt.expand <- model.matrix(model, untrtobs)
  Y.untrt <- untrtobs$Y

  list(noimp = trtmat, noimp.Y = Y.trt, imp = untrt.expand, imp.Y = Y.untrt)

}

