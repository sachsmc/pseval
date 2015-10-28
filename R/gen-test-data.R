
trunc01 <- function(x){

  pmax(pmin(x, 1), 0)

}

expit <- function(x) exp(x)/(1 + exp(x))

#' Generate sample data used for testing
#'
#' @export

generate_gh_data <- function(n){

  Z <- rbinom(n, 1, .5)
  X <- rnorm(n)  ## pre-treatment covariate

  ## generate S conditional on Z

  S.0 <- X + rnorm(n, sd = .1)
  S.1 <- 2 + X + rnorm(n, sd = .1)

  risk.obs <- expit(-1 + 0 * Z + 2 * S.1 - 1 * S.1 * Z)
  risk.0 <- expit(-1 + 2 * S.1 )
  risk.1 <- expit(-1 + 2 * S.1 - 1 * S.1 * 1)

  Y.0 <- rbinom(n, 1, trunc01(risk.0))
  Y.1 <- rbinom(n, 1, trunc01(risk.1))
  Y.obs <- ifelse(Z == 1, Y.1, Y.0)

  ## CPV measure noisy S.1 at the end of the study for placebo subjects non-event

  CPV <- S.1 + rnorm(n, sd = .1)
  CPV[Z == 1 | Y.obs == 1] <- NA

  ## BSM measure noisy S.0 at start

  BSM <- S.0 + rnorm(n, sd = .1)
  S.1[Z == 0] <- NA
  S.0[Z == 1] <- NA

  ## make categorical variables
  qwantz <- c(-Inf, quantile(c(S.0, S.1), c(.25, .5, .75), na.rm = TRUE), Inf)
  S.1.cat <- cut(S.1, qwantz)
  S.0.cat <- cut(S.0, qwantz)

  BIP.cat <- cut(X, c(-Inf, quantile(X, c(.25, .5, .75)), Inf))

  data.frame(Z, X, CPV, BSM, S.0, S.1, risk.obs, risk.0, risk.1, Y.obs, Y.0, Y.1, S.1.cat, S.0.cat, BIP.cat)

}


test <- generate_gh_data(1000)




