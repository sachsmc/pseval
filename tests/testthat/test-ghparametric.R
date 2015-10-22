library(pseval)

test_that("Gilbert Hudgens estimates work", {

    ghdat <- generate_gh_data(400)
    ghdes <- psdesign(ghdat, Z = Z, Y = Y.obs, S = S.1, BIP = X, CPV = CPV, BSM = BSM)
    ghdes2 <- ghdes + impute_parametric(S.1 ~ BIP)
    ghdes3 <- ghdes2 + risk_binary(D = 500, risk = risk.expit)


    truepar <- c(-1, 2, 0, -1)

    test2 <- ps_estimate(ghdes3, start = rep(0, 4), method = "BFGS")

    #test <- ps_bootstrap(ghdes3, n.boots = 50, start = truepar, method = "BFGS")


    estpar <- optim(truepar, fn = ghdes3$likelihood, method = "Nelder",
                    control = list(fnscale = -1, maxit = 5000))




  # classes

  expect_true(inherits(ghdes, c("ps", "psdesign")))
  expect_true(inherits(ghdes2, c("ps", "psdesign")))
  expect_true(inherits(ghdes3, c("ps", "psdesign")))

})