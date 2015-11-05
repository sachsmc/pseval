## ----knitup, include = FALSE---------------------------------------------
library(knitr)
library(printr)
set.seed(205)

## ----eval = FALSE--------------------------------------------------------
#  devtools::install_github("sachsmc/pseval")

## ----setup---------------------------------------------------------------
library(pseval)
library(survival)

## ------------------------------------------------------------------------
fakedata <- generate_example_data(n = 500)
head(fakedata)

## ----psdesign------------------------------------------------------------
binary.ps <- psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP)
binary.ps

## ----helpimp-------------------------------------------------------------
?add_integration

## ----impp----------------------------------------------------------------
binary.ps <- binary.ps + integrate_parametric(S.1 ~ BIP)
binary.ps

## ----imp2----------------------------------------------------------------
binary.ps + integrate_parametric(S.0 ~ BIP)

## ----imp3----------------------------------------------------------------
library(splines)
binary.ps + integrate_parametric(S.1 ~ BIP^2)
binary.ps + integrate_parametric(S.1 ~ bs(BIP, df = 3))

## ----riskhelp------------------------------------------------------------
?add_riskmodel

## ----riskbin-------------------------------------------------------------
binary.ps <- binary.ps + risk_binary(model = Y ~ S.1 * Z, D = 2000, risk = risk.logit)
binary.ps

