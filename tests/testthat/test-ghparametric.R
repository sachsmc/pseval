library(pseval)
library(splines)

test_that("Gilbert Hudgens estimates work", {

    ghdat <- generate_gh_data(200)
    ghdes <- psdesign(ghdat, Z = Z, Y = Y.obs, S = S.1, BIP = X, CPV = CPV, BSM = BSM)

    ghdes.cat <- psdesign(ghdat, Z = Z, Y = Y.obs, S = S.1.cat, BIP = BIP.cat)
    ghdes.cat2 <- ghdes.cat + impute_nonparametric(S.1 ~ BIP)
    test.cat <- ps_estimate(ghdes.cat2 + risk_binary(D = 1000, risk = risk.expit), method = "BFGS")

    ghdes2 <- ghdes + impute_parametric(S.1 ~ BIP)
    ghdes2b <- ghdes2 + impute_semiparametric(S.1 ~ BIP, S.1 ~ 1)

    ghdes3 <- ghdes2 + risk_binary(D = 1000, risk = risk.expit)
    ghdes3b <- ghdes2b + risk_binary(D = 1000, risk = risk.expit)


    truepar <- c(-1, 2, 0, -1)

    test2 <- ps_estimate(ghdes3, start = rep(0, 4), method = "BFGS")
    test2.boot <- ps_bootstrap(test2, n.boots = 100, start = test2$estimates$par, method = "BFGS")

    test2b <- ps_estimate(ghdes3b, start = truepar, method = "BFGS")
    #test <- ps_bootstrap(ghdes3, n.boots = 50, start = truepar, method = "BFGS")


  # classes

  expect_true(inherits(ghdes, c("ps", "psdesign")))
  expect_true(inherits(ghdes2, c("ps", "psdesign")))
  expect_true(inherits(ghdes3, c("ps", "psdesign")))

})