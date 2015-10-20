library(pseval)

test_that("Gilbert Hudgens estimates work", {

    ghdat <- generate_gh_data(400)
    ghdes <- psdesign(ghdat, Z, Y.obs, surrogate = S.1, BIP = X)
    ghdes2 <- ghdes + impute_parametric(S.1 ~ X)
    ghdes3 <- ghdes2 + risk_logistic(D = 500)


    truepar <- c(-1, 2, 0, -1)

    test2 <- ps_estimate(ghdes3, start = truepar, method = "BFGS")

    test <- ps_bootstrap(ghdes3, n.boots = 50, start = truepar, method = "BFGS")


    estpar <- optim(truepar, fn = ghdes3$likelihood, method = "Nelder",
                    control = list(fnscale = -1, maxit = 5000))




  # classes

  expect_true(inherits(ghdes, c("ps", "psdesign")))
  expect_true(inherits(ghdes2, c("ps", "psdesign")))
  expect_true(inherits(ghdes3, c("ps", "psdesign")))

})