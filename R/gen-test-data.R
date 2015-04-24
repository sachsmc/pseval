
trunc01 <- function(x){

  pmax(pmin(x, 1), 0)

}

expit <- function(x) exp(x)/(1 + exp(x))

#' Generate sample data used for testing
#'

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

  data.frame(Z, X, S.0, S.1, risk.obs, risk.0, risk.1, Y.obs, Y.0, Y.1)

}


test <- generate_gh_data(1000)




