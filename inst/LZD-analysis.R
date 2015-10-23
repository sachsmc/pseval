
dat.out <- read.csv("inst/LZD_analysisreadyData_Oct7th15.csv", stringsAsFactors = FALSE)
dat.out$SImp <- with(dat.out, ifelse(is.na(S_I1), S_D1, S_I1))

dat <- dat.out

library(survival)

sr <- survreg(Surv(Time_AE, event_AE) ~ SImp * Z, data= dat.out)

trtmat <- model.matrix(~ SImp * Z, data = dat.out)
Y.trt <- dat.out$Time_AE
delt.trt <- dat.out$event_AE

mu <- sapply(dat[, c("S_I1", "S_D1")], mean, na.rm = TRUE)
sd <- sapply(dat[, c("S_I1", "S_D1")], sd, na.rm = TRUE)
rho.init <- cor(dat[, c("M1", "M2")], use = "pairwise")[1, 2]

ghdes <- psdesign(dat, Z, Time_AE, surrogate = S_I1, BIP = SImp)
ghdes2 <- ghdes + impute_bivnorm(mu = mu, sd = sd, rho = rho.init)

dat$SImp[dat$Z == 0] <- rowMeans(ghdes2$icdf_sbarw(runif(500)))


likelihood <- function(beta){

  gamma0<-beta[1]
  beta0<-beta[-1]
  shapepram<-exp(gamma0)
  scalepram<-exp(trtmat %*% beta0)

  trtlike <- ((shapepram/scalepram)*(Y.trt/scalepram)^(shapepram-1))^delt.trt*exp(-(Y.trt/scalepram)^shapepram)

  -1 * sum(log(trtlike))

}

start <- c(.856, 6.376, -.326, -.528, .134)
optim(start, likelihood)$par


swank.one <- function(dat){
    mu <- sapply(dat[, c("S_I1", "S_D1")], mean, na.rm = TRUE)
    sd <- sapply(dat[, c("S_I1", "S_D1")], sd, na.rm = TRUE)
    rho.init <- cor(dat[, c("M1", "M2")], use = "pairwise")[1, 2]

    ghdes <- psdesign(dat, Z = Z, Y = Surv(Time_AE, event_AE), S = S_I1, BIP = S_D1)
    ghdes2 <- ghdes + impute_bivnorm(x = S.1, mu = mu, sd = sd, rho = rho.init)
    ghdes3 <- ghdes2 + risk_weibull(Y ~ S.1 * Z)
    ghdes3b <- ghdes2 + risk_exponential(Y ~ S.1 * Z)



    start <- c(.856, 6.376, -.326, -.528, .134)

    est1 <- ps_estimate(ghdes3, start = rep(0, 5), method = "BFGS")
    est1b <- ps_estimate(ghdes3b, start = rep(0, 4), method = "BFGS")

    ghdes2.lower <- ghdes + impute_bivnorm(mu = mu, sd = sd, rho = .2)
    ghdes3.lower <- ghdes2.lower + risk_weibull(Time_AE ~ S_I1 * Z, cens = "event_AE")

    est1.lower <- optim(start, fn = ghdes3.lower$likelihood, method = "BFGS",
                  control = list(fnscale = -1, maxit = 2000))

    list(est.rho = est1, rho.2 = est1.lower)

}


system.time(orig.fit <- swank.one(dat.out))
rownames(dat.out) <- paste(dat.out$PID)

resample <- function(){

  trted <- sample(dat.out[dat.out$Z == 1, "PID"], size = sum(dat.out$Z == 1), replace = TRUE)
  untrted <- sample(dat.out[dat.out$Z == 0, "PID"], size = sum(dat.out$Z == 0), replace = TRUE)

  dat.out[paste(c(trted, untrted)), ]

}


boots <- vector(mode = "list", length = 200)
for(i in 1:length(boots)){

  boots[[i]] <- swank.one(resample())

}

pnames <- c("log(shape)", "intercept", "S1", "Z", "S1:Z")

swag <- function(l){

  if(l$est.rho$convergence != 0){
    b1 <- rep(NA, length(pnames))
  } else{
    b1 <- l$est.rho$par
  }
  if(l$rho.2$convergence != 0){
    b2 <- rep(NA, length(pnames))
  } else{
    b2 <- l$est.rho$par
  }

  res <- c(b1, b2)
  names(res) <- c(paste(pnames, "est.rho", sep = "."), paste(pnames, "rho.2", sep = "."))
  res
}

mat.res <- t(sapply(boots, swag))

write.csv(mat.res, file = paste0("inst/LZD-bootstrap-results", Sys.Date(), ".csv"))


## Ogawa

swank.one <- function(dat){
  mu <- sapply(dat[, c("S_I1", "S_D1")], mean, na.rm = TRUE)
  sd <- sapply(dat[, c("S_I1", "S_D1")], sd, na.rm = TRUE)
  rho.init <- cor(dat[, c("M1", "M2")], use = "pairwise")[1, 2]

  ghdes <- psdesign(dat, Z, Ogawa_time, surrogate = S_I1, BIP = S_D1)
  ghdes2 <- ghdes + impute_bivnorm(mu = mu, sd = sd, rho = rho.init)
  ghdes3 <- ghdes2 + risk_weibull(Time_AE ~ S_I1 * Z, cens = "Ogawa_convert")

  start <- c(.8892, 5.16, -.085, -.9288, .3473)
  est1 <- optim(start, fn = ghdes3$likelihood, method = "BFGS",
                control = list(fnscale = -1, maxit = 2000))


  ghdes2.lower <- ghdes + impute_bivnorm(mu = mu, sd = sd, rho = .2)
  ghdes3.lower <- ghdes2.lower + risk_weibull(Ogawa_time ~ S_I1 * Z, cens = "Ogawa_convert")

  est1.lower <- optim(start, fn = ghdes3.lower$likelihood, method = "BFGS",
                      control = list(fnscale = -1, maxit = 2000))

  list(est.rho = est1, rho.2 = est1.lower)

}


system.time(orig.fit <- swank.one(dat.out))
rownames(dat.out) <- paste(dat.out$PID)

resample <- function(){

  trted <- sample(dat.out[dat.out$Z == 1, "PID"], size = sum(dat.out$Z == 1), replace = TRUE)
  untrted <- sample(dat.out[dat.out$Z == 0, "PID"], size = sum(dat.out$Z == 0), replace = TRUE)

  dat.out[paste(c(trted, untrted)), ]

}


boots <- vector(mode = "list", length = 200)
for(i in 1:length(boots)){

  boots[[i]] <- swank.one(resample())

}

pnames <- c("log(shape)", "intercept", "S1", "Z", "S1:Z")

swag <- function(l){

  if(l$est.rho$convergence != 0){
    b1 <- rep(NA, length(pnames))
  } else{
    b1 <- l$est.rho$par
  }
  if(l$rho.2$convergence != 0){
    b2 <- rep(NA, length(pnames))
  } else{
    b2 <- l$est.rho$par
  }

  res <- c(b1, b2)
  names(res) <- c(paste(pnames, "est.rho", sep = "."), paste(pnames, "rho.2", sep = "."))
  res
}

mat.res <- t(sapply(boots, swag))

write.csv(mat.res, file = paste0("inst/LZD-bootstrap-results-ogawa", Sys.Date(), ".csv"))


