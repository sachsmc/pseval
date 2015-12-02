library(pseval)
library(survival)
library(splines)

test_that("Gilbert Hudgens estimates work", {

  fakedata <- generate_example_data(n = 500)
  binary.est <- psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP) +
    integrate_parametric(S.1 ~ BIP) +
    risk_binary(model = Y ~ S.1 * Z, D = 1000, risk = risk.logit) +
    ps_estimate(method = "BFGS") + ps_bootstrap(method = "BFGS", n.boots = 10)
  summary(binary.est)

})