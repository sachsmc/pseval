#' Logistic risk model for binary outcome
#'
#' @param model Formula specifying the risk model
#' @param D number of samples for the simulated annealing integration


risk_logistic <- function(model = Y ~ S.1 * Z, D = 5000, ...){

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    risk <- function(beta, mat) {

      exp(mat %*% beta)/(1 + exp(mat %*% beta))

    }

    trtmat <- model.matrix(model, psdesign$augdata[psdesign$augdata$Z == 1, ])
    Y.trt <- psdesign$augdata[psdesign$augdata$Z == 1, ]$Y

    untrtsamp <- c((psdesign$icdf_sbarw(runif(D))))
    dex <- (1:nrow(psdesign$augdata))[psdesign$augdata$Z == 0]
    untrtobs <- psdesign$augdata[rep(dex, D), ]
    untrtobs$S.1 <- untrtsamp

    untrt.expand <- model.matrix(model, untrtobs)
    Y.untrt <- untrtobs$Y

    likelihood <- function(beta){

      trted <- risk(beta, trtmat)^Y.trt *
        (1 - risk(beta, trtmat))^(1 - Y.trt)

      # for each W in untrted, sample a bunch from cdf_sbarw and take the mean

      untrted <- matrix(risk(beta, untrt.expand)^Y.untrt *
        (1 - risk(beta, untrt.expand))^(1 - Y.untrt), nrow = D, byrow = TRUE)

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
#' @param model Formula specifying the risk model
#' @param cens Symbol identifying the censoring indicator
#' @param D number of samples for simulated annealing
#' @param ...
#'
risk_weibull <- function(model = Y ~ S.1 * Z, cens = delt, D = 5000, ... ){


  arglist <- as.list(match.call())
  rval <- function(psdesign){

    trtmat <- model.matrix(model, psdesign$augdata[psdesign$augdata$Z == 1, ])
    Y.trt <- psdesign$augdata[psdesign$augdata$Z == 1, ]$Y
    delt.trt <- psdesign$augdata[psdesign$augdata$Z == 1, cens]

    untrtsamp <- c((psdesign$icdf_sbarw(runif(D))))
    dex <- (1:nrow(psdesign$augdata))[psdesign$augdata$Z == 0]
    untrtobs <- psdesign$augdata[rep(dex, D), ]
    untrtobs[ , attr(terms.formula(model), "term.labels")[1]] <- untrtsamp

    untrt.expand <- model.matrix(model, untrtobs)
    Y.untrt <- untrtobs$Y
    delt.untrt <- untrtobs[, cens]

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
