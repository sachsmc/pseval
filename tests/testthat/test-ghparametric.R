library(pseval)

test_that("Gilbert Hudgens estimates work", {


  ghdat <- generate_gh_data(1000)
  ghdes <- psdesign(ghdat, Z, Y.obs, surrogate = S.1, BIP = X)
  ghdes2 <- ghdes + impute_parametric(S.1 ~ X)
  ghdes3 <- ghdes2 + risk_logistic(D = 10000)

  truepar <- c(0, 2, 0, -1)

  estpar <- optim(truepar, fn = ghdes3$likelihood, method = "Nelder",
                  control = list(fnscale = -1, reltol = 1e-5))



  # classes

  expect_true(inherits(ghdes, c("ps", "psdesign")))
  expect_true(inherits(ghdes2, c("ps", "psdesign")))
  expect_true(inherits(ghdes3, c("ps", "psdesign")))

})