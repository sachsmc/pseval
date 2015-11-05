library(pseval)
library(survival)
library(splines)

test_that("Gilbert Hudgens estimates work", {

  estpar <- NULL
  for(i in 1:100){
    ghdat <- generate_example_data(200)
    ghdes <- psdesign(ghdat, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP)

    check <- ghdes + integrate_parametric(S.1 ~ BIP) + risk_binary(D = 500) + ps_estimate(method = "BFGS")
    estpar <- rbind(estpar, check$estimates$par)

  }

    ghdes.surv <- psdesign(ghdat, Z = Z, Y = Surv(time.obs, event.obs), S = S.obs, BIP = BIP)
    ghdes.cat <- psdesign(ghdat, Z = Z, Y = Y.obs, S = S.obs.cat, BIP = BIP.cat)

    check <- ghdes.surv + integrate_parametric(S.1 ~ BIP) + risk_exponential(D = 1000) + ps_estimate()

    truepar <- c(-1, 2, 0, -1)

  # classes

  expect_true(inherits(ghdes, c("ps", "psdesign")))
  expect_true(inherits(ghdes2, c("ps", "psdesign")))
  expect_true(inherits(ghdes3, c("ps", "psdesign")))

})