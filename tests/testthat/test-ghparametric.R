library(pseval)

test_that("Gilbert Hudgens estimates work", {

  ghdat <- generate_gh_data(1000)
  ghdes <- psdesign(ghdat, Z, Y.obs, BIP = X)
  ghdes2 <- ghdes + impute_parametric(S.1 ~ X)

  # classes

  expect_true(inherits(ghdes, c("ps", "psdesign")))
  expect_true(inherits(ghdes2, c("ps", "psdesign")))

  # risk model
  ghdes <- add_risk_model(ghdes, risk_logistic())

})