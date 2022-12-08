rm(list=ls())
source('EM.R')
dataDir <- paste0(getwd(), '/Simulations/data')

# Define fitting functions ------------------------------------------------
emfit <- function(data, k, family){
  long.formulas <- vector('list', k)
  for(kk in 1:k) long.formulas[[kk]] <- as.formula(paste0('Y.', kk, '~ time + cont + bin + (1 + time|id)'))
  surv.formula <- Surv(survtime, status) ~ bin
  families <- as.list(rep(family, k))
  fit <- tryCatch(suppressMessages(
    EM(long.formulas, surv.formula, data = data, family = families,
       control = list(hessian = 'manual'))),
    error = function(e) NULL
  )
  fit
}

joineRMLfit <- function(data, k){
  long <- rand <- vector('list', k)
  for(kk in 1:k){
    long[[kk]] <- as.formula(paste0('Y.', kk, '~time + cont + bin'))
    rand[[kk]] <- ~ 1 + time | id
  }
  fit <- suppressMessages(joineRML::mjoint(
    formLongFixed = long,
    formLongRandom = rand,
    formSurv = Surv(survtime, status) ~ bin,
    timeVar = 'time', data = data, 
    control = list(
      convCrit = 'rel',
      tol.em = 5e-3,
      type = 'sobol', tol2 = 1e-2
    )
  ))
  summary(fit)
}

JMbayes2fit <- function(data, k, family){
  M <- vector('list', k)
  if(family == 'gaussian'){
    for(kk in 1:k){
      long <- as.formula(paste0('Y.', kk, '~time + cont + bin'))
      random <- ~time|id
      M[[kk]] <- lme(fixed = long, random = random, 
                     data = data, method = 'ML', control = lmeControl(opt='optim'))
    }
  }else{
    for(kk in 1:k){
      long <- as.formula(paste0('Y.', kk, '~time + cont + bin'))
      random <- ~time|id
      TMBfit <- glmmTMB(as.formula(paste0('Y.',kk,'~time+cont+bin+(1+time|id)')),
                        data = data, family = family)
      inits <- list(
        betas = glmmTMB::fixef(TMBfit)$cond,
        D = matrix(glmmTMB::VarCorr(TMBfit)$cond$id, 2, 2)
      )
      M[[kk]] <- mixed_model(fixed = long, random = random,
                             data = data, family = family,
                             initial_values = inits)
    }
  }
  
  survdata <- data[!duplicated(data[, 'id']), ]
  ph <- coxph(Surv(survtime, status) ~ bin, survdata)
  
  fit <- tryCatch(
    jm(ph, M, time_var = 'time', id_var = 'id', data_Surv = survdata,),
    error=function(e) NULL)
  if(!is.null(fit)){
    s <- summary(fit)
    rtn <- list()
    rtn$Outcome1 <- s$Outcome1
    rtn$Outcome2 <- s$Outcome2
    rtn$Outcome3 <- s$Outcome3
    rtn$survival <- s$Survival
    rtn$comp.time <- s$time
    return(rtn)
  }else{ 
    return(NA)
  }
}

INLAfit <- function(data, k, family){
  long.formulas <- vector('list', k)
  for(kk in 1:k) long.formulas[[kk]] <- as.formula(paste0('Y.', kk, '~ time + cont + bin + (1 + time|id)'))
  survdata <- data[!duplicated(data[, 'id']), ]
  fail <<- inla.surv(time = c(survdata$survtime), event = c(survdata$status)) # send to global env or joint() doesn't work.
  assocs <- as.list(rep('SRE', k))
  JMINLA <- joint(
    formLong = long.formulas,
    formSurv = fail ~ bin,
    dataLong = data, dataSurv = survdata,
    id = 'id', timeVar = 'time', corLong = T,
    family = rep(family, k), basRisk = 'rw1',
    assoc = assocs,
    control = list(int.strategy='eb',
                   priorRandom = list(r = 2*k, R = 1),
                   priorAssoc = list(mean = 0, prec = 0.16),
                   priorFixed = list(mean = 0, prec = 0.16, 
                                     mean.intercept = 0,
                                     prec.intercept = 0.16))
  )
  s <- summary(JMINLA)
  longi <- do.call(rbind, s$FixedEff)
  survi <- s$SurvEff[[1]]
  gammas <- s$AssocLS
  time <- s$cpu.used   # Just running time?
  
  list(
    fixed = longi,
    survival = survi,
    gamma = gammas,
    comp.time = time
  )
}

.loader <- function(file){
  assign('data', get(load(file)))
  lapply(data, function(x) x$data)
}

# K=1 ---------------------------------------------------------------------
data1 <- .loader('Simulations/data/gaussianK-1.RData')
fit1 <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data1[[i]]
  fit1[[i]] <- emfit(d, 1, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}
save(fit1, file =  'Simulations/fits/gaussianK-1.RData')

# K=2 ---------------------------------------------------------------------
data2 <- .loader('Simulations/data/gaussianK-2.RData')
fit2 <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data2[[i]]
  fit2[[i]] <- emfit(d, 2, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}
save(fit2, file = 'Simulations/fits/gaussianK-2.RData')

# K=3 ---------------------------------------------------------------------
data3 <- .loader('Simulations/data/gaussianK-3.RData')
fit3 <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data3[[i]]
  fit3[[i]] <- emfit(d, 3, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}
save(fit3, file = 'Simulations/fits/gaussianK-3.RData')


# joineRML fits -----------------------------------------------------------
# K=1 joineRML ------------------------------------------------------------
data1 <- .loader('Simulations/data/gaussianK-1.RData')
fit1.jML <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data1[[i]]
  fit1.jML[[i]] <- joineRMLfit(d, 1)
  utils::setTxtProgressBar(pb, i)
}
save(fit1.jML, file =  'Simulations/fits/gaussianK-1_joineRML.RData')

# K=2 joineRML ------------------------------------------------------------
data2 <- .loader('Simulations/data/gaussianK-2.RData')
fit2.jML <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data2[[i]]
  fit2.jML[[i]] <- joineRMLfit(d, 2)
  utils::setTxtProgressBar(pb, i)
}
save(fit2.jML, file =  'Simulations/fits/gaussianK-2_joineRML.RData')

# K=3 joineRML ------------------------------------------------------------
data3 <- .loader('Simulations/data/gaussianK-3.RData')
fit3.jML <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data3[[i]]
  fit3.jML[[i]] <- joineRMLfit(d, 3)
  utils::setTxtProgressBar(pb, i)
}
save(fit3.jML, file =  'Simulations/fits/gaussianK-3_joineRML.RData')


# JMbayes2 ----------------------------------------------------------------
# K=1 JMbayes 2------------------------------------------------------------
library(JMbayes2)
data2 <- .loader('Simulations/data/gaussianK-2.RData')
fit2.JMb <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data2[[i]]
  fit2.JMb[[i]] <- JMbayes2fit(d, 2, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}
save(fit2.JMb, file = 'Simulations/fits/gaussianK-2_JMbayes2-2.RData')

# K=2 JMbayes 2------------------------------------------------------------
  
# END GAUSSIAN --------
# Poisson -----------------------------------------------------------------
# rm(data1, data2, data3) already done!
data1 <- .loader('Simulations/data/poissonK-1.RData')
data2 <- .loader('Simulations/data/poissonK-2.RData')
data3 <- .loader('Simulations/data/poissonK-3.RData')
# fit1 <- fit2 <- fit3 <- fit1.JMb <- fit2.JMb <- fit3.JMb <- vector('list', 100)
# pb <- utils::txtProgressBar(max = 100, style = 3)
# for(i in 1:100){
#   d1 <- data1[[i]]; d2 <- data2[[i]]; d3 <- data3[[i]]
#   fit1[[i]] <- emfit(d1, 1, 'poisson')
#   fit2[[i]] <- emfit(d2, 2, 'poisson')
#   fit3[[i]] <- emfit(d3, 3, 'poisson')
#   utils::setTxtProgressBar(pb, i)
# }
# save(fit1, file = 'Simulations/fits/poissonK-1.RData')
# save(fit2, file = 'Simulations/fits/poissonK-2.RData')
# save(fit3, file = 'Simulations/fits/poissonK-3.RData')

# JMbayes2
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){ # 21/7/22 do k=2 poisson jmb2
  d1 <- data1[[i]]; d2 <- data2[[i]]; d3 <- data3[[i]]
  # fit1.JMb[[i]] <- JMbayes2fit(d1, 1, 'poisson')
  fit2.JMb[[i]] <- JMbayes2fit(d2, 2, 'poisson')
  # fit3.JMb[[i]] <- JMbayes2fit(d3, 3, 'poisson')
  utils::setTxtProgressBar(pb, i)
}
# save(fit1.JMb, file = 'Simulations/fits/poissonK-1_JMbayes2-2.RData')
save(fit2.JMb, file = 'Simulations/fits/poissonK-2_JMbayes2-2.RData')
# save(fit3.JMb, file = 'Simulations/fits/poissonK-3_JMbayes2-2.RData')


# INLAjoint ---------------------------------------------------------------
# INLA
library(INLA)
library(INLAjoint)
inla.setOption(inla.mode="experimental")
data1 <- .loader('Simulations/data/gaussianK-1.RData')
data2 <- .loader('Simulations/data/gaussianK-2.RData')
data3 <- .loader('Simulations/data/gaussianK-3.RData')

pb <- utils::txtProgressBar(max = 100, style = 3)
fit1.INLA <- fit2.INLA <- fit3.INLA <- vector('list', 100)
for(i in 1:100){
  d1 <- data1[[i]]; d2 <- data2[[i]];# d3 <- data3[[i]]
  #fit1.INLA[[i]] <- tryCatch(suppressMessages(INLAfit(d1, 1, 'gaussian')), error = function(e) NULL)
  fit2.INLA[[i]] <- tryCatch(suppressMessages(INLAfit(d2, 2, 'gaussian')), error = function(e) NULL)
  # fit3.INLA[[i]] <- INLAfit(d3, 3, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}
#save(fit1.INLA, file = 'Simulations/fits/gaussianK-1_INLA.RData')
save(fit2.INLA, file = 'Simulations/fits/gaussianK-2_INLA.RData')
# save(fit3.INLA, file = 'Simulations/fits/gaussianK-3_INLA.RData')

# Poisson
data1 <- .loader('Simulations/data/poissonK-1.RData')
data2 <- .loader('Simulations/data/poissonK-2.RData')
data3 <- .loader('Simulations/data/poissonK-3.RData')

pb <- utils::txtProgressBar(max = 100, style = 3)
fit1.INLA <- fit2.INLA <- fit3.INLA <- vector('list', 100)
for(i in 1:100){
  d1 <- data1[[i]]; 
  d2 <- data2[[i]];
 # d3 <- data3[[i]]
  fit1.INLA[[i]] <- INLAfit(d1, 1, 'poisson')
  fit2.INLA[[i]] <- INLAfit(d2, 2, 'poisson')
  #fit3.INLA[[i]] <- INLAfit(d3, 3, 'poisson')
  utils::setTxtProgressBar(pb, i)
}
save(fit1.INLA, file = 'Simulations/fits/poissonK-1_INLA.RData')
save(fit2.INLA, file = 'Simulations/fits/poissonK-2_INLA.RData')
save(fit3.INLA, file = 'Simulations/fits/poissonK-3_INLA.RData')

# END COUNTS --------
# Binary ------------------------------------------------------------------
data1 <- .loader('Simulations/data/binomialK-1.RData')
data2 <- .loader('Simulations/data/binomialK-2.RData')
data3 <- .loader('Simulations/data/binomialK-3.RData')

fit1 <- fit2 <-fit3 <-  vector('list', 100)
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d1 <- data1[[i]]; d2 <- data2[[i]]; d3 <- data3[[i]]
  #fit1[[i]] <- emfit(d1, 1, 'binomial')
  # fit2[[i]] <- emfit(d2, 2, 'binomial')
  fit3[[i]] <- emfit(d3, 3, 'binomial')
  utils::setTxtProgressBar(pb, i)
}
save(fit3, file = 'Simulations/fits/binomialK-3.RData')

library(INLA)
library(INLAjoint)
fit1.INLA <- fit2.INLA <- fit3.INLA <- vector('list', 100)
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d1 <- data1[[i]]; d2 <- data2[[i]]; d3 <- data3[[i]]
  # fit1.INLA[[i]] <- INLAfit(d1, 1, 'binomial')
  fit2.INLA[[i]] <- INLAfit(d2, 2, 'binomial')
  # fit3[[i]] <- emfit(d3, 3, 'binomial')
  utils::setTxtProgressBar(pb, i)
}
save(fit2.INLA, file = 'Simulations/fits/binomialK-2_INLA.RData')
