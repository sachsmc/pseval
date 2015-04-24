#' Logistic risk model for binary outcome
#'
#' @param model Formula specifying the risk model
#' @param D number of samples for the simulated annealing integration


risk_logistic <- function(model = Y ~ S.1 * Z, D = 5000, ...){

  rval <- function(psdesign){

    risk <- function(beta, mat) {

      exp(mat %*% beta)/(1 + exp(mat %*% beta))

    }

    trtmat <- model.matrix(model, psdesign$augdata[psdesign$augdata$Z == 1, ])
    Y.trt <- psdesign$augdata[psdesign$augdata$Z == 1, ]$Y

    untrtsamp <- c(t(psdesign$icdf_sbarw(runif(D))))
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
    psdesign$riskmodel <- "logistic"
    psdesign$nparam <- ncol(trtmat)
    psdesign$riskformula <- model

    psdesign

  }
  ## return a likelihood closure

  class(rval) <- c("ps", "riskmodel")
  rval
}